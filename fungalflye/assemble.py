import subprocess
from pathlib import Path
import shutil
import shutil as sh


def run(cmd):
    print(f"\n[fungalflye] Running: {cmd}\n")
    subprocess.run(cmd, shell=True, check=True)


def check_dependencies():

    tools = ["flye", "minimap2", "racon", "seqkit"]

    missing = []

    for t in tools:
        if sh.which(t) is None:
            missing.append(t)

    if missing:
        print("\n❌ Missing required tools:")
        for m in missing:
            print(f"  - {m}")

        print("\nInstall with:")
        print("conda install -c bioconda flye minimap2 racon seqkit filtlong\n")

        raise SystemExit(1)


def prune_contained_contigs(fasta, out_fasta, threads=8):

    print("\n[fungalflye] Pruning redundant contigs (>95% contained)\n")

    fasta = Path(fasta)
    out_fasta = Path(out_fasta)

    paf = out_fasta.with_suffix(".self.paf")

    run(f"minimap2 -x asm10 -t {threads} {fasta} {fasta} > {paf}")

    remove = set()

    with open(paf) as f:
        for line in f:

            cols = line.strip().split()

            q = cols[0]
            t = cols[5]

            if q == t:
                continue

            qlen = int(cols[1])
            tlen = int(cols[6])

            matches = int(cols[9])
            aln_len = int(cols[10])

            identity = matches / aln_len if aln_len > 0 else 0
            coverage = aln_len / qlen if qlen > 0 else 0

            if identity > 0.95 and coverage > 0.95:
                if qlen < tlen:
                    remove.add(q)

    with open(out_fasta, "w") as out:
        keep = True
        with open(fasta) as f:
            for line in f:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    keep = name not in remove
                if keep:
                    out.write(line)

    paf.unlink(missing_ok=True)

    print(f"[fungalflye] Removed {len(remove)} redundant contigs")


def run_assembly(
    reads,
    genome_size,
    outdir,
    threads=8,
    min_read_len=0,
    downsample_cov=0,
):

    check_dependencies()

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    reads = Path(reads)

    flye_dir = outdir / "flye"
    filtered_reads = outdir / "reads.filtered.fastq"
    paf = outdir / "reads.paf"

    racon_fasta = outdir / "racon.fasta"
    pruned_fasta = outdir / "pruned.fasta"
    final_fasta = outdir / "final.fasta"

    reads_used = reads

    # filtering
    if min_read_len > 0:
        run(f"seqkit seq -m {min_read_len} {reads} > {filtered_reads}")
        reads_used = filtered_reads

    # Flye
    run(
        f"flye --nano-hq {reads_used} "
        f"--genome-size {genome_size} "
        f"--threads {threads} "
        f"--iterations 3 "
        f"--asm-coverage 60 "
        f"--keep-haplotypes "
        f"-o {flye_dir}"
    )

    assembly = flye_dir / "assembly.fasta"

    # mapping
    run(
        f"minimap2 -x map-ont -t {threads} "
        f"{assembly} {reads_used} > {paf}"
    )

    # racon
    run(f"racon {reads_used} {paf} {assembly} > {racon_fasta}")

    prune_contained_contigs(racon_fasta, pruned_fasta, threads)

    shutil.copy(pruned_fasta, final_fasta)

    print("\n🧬 Assembly Complete")
    print(f"Final assembly: {final_fasta}\n")

    return str(final_fasta)