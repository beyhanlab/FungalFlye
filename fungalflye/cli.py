import typer
import subprocess
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

from .qc import run_qc

app = typer.Typer(help="FungalFlye Long-read fungal genome assembly pipeline")


# -------------------------
# helper runner
# -------------------------
def run(cmd):
    print(f"\n[fungalflye] Running: {cmd}\n")
    subprocess.run(cmd, shell=True, check=True)


# -------------------------
# read length analysis
# -------------------------
def analyze_reads(reads, outdir):

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    lengths_file = outdir / "read_lengths.tsv"

    run(f"seqkit fx2tab -n -l {reads} > {lengths_file}")

    df = pd.read_csv(lengths_file, sep="\t", header=None)
    lengths = df[1]

    total_reads = len(lengths)
    total_bases = lengths.sum()

    # N50
    sorted_lengths = sorted(lengths, reverse=True)
    cumsum = 0
    n50 = 0

    for L in sorted_lengths:
        cumsum += L
        if cumsum >= total_bases / 2:
            n50 = L
            break

    # plot
    plt.figure(figsize=(6, 4))
    plt.hist(lengths, bins=100)
    plt.xlabel("Read length (bp)")
    plt.ylabel("Count")
    plt.title("Read length distribution")
    plt.tight_layout()
    plt.savefig(outdir / "read_length_histogram.png")

    return total_reads, total_bases, n50, lengths


# -------------------------
# filtering preview
# -------------------------
def preview_filter(lengths, cutoff):

    kept = lengths[lengths >= cutoff]
    removed = lengths[lengths < cutoff]

    return len(kept), len(removed)


# -------------------------
# assembly core
# -------------------------
@app.command()
def assemble(
    reads: str,
    gsize: str,
    outdir: str = "fungalflye_out",
    threads: int = 8,
):

    Path(outdir).mkdir(exist_ok=True)

    typer.echo("\n🧬 Running Flye assembly\n")

    run(
        f"flye --nano-hq {reads} "
        f"--genome-size {gsize} "
        f"--threads {threads} "
        f"--iterations 3 "
        f"--asm-coverage 60 "
        f"--keep-haplotypes "
        f"-o {outdir}/flye"
    )

    typer.echo("\n🧬 Mapping reads for polishing\n")

    run(
        f"minimap2 -x map-ont "
        f"{outdir}/flye/assembly.fasta {reads} "
        f"> {outdir}/reads.paf"
    )

    typer.echo("\n🧬 Running Racon polishing\n")

    run(
        f"racon {reads} "
        f"{outdir}/reads.paf "
        f"{outdir}/flye/assembly.fasta "
        f"> {outdir}/racon.fasta"
    )

    typer.echo("\n🧬 Sorting contigs by length\n")

    run(
        f"seqkit fx2tab -n -l {outdir}/racon.fasta "
        f"| sort -k2,2nr > {outdir}/contig_sizes_sorted.tsv"
    )

    typer.echo(f"\n✅ Final polished assembly: {outdir}/racon.fasta")
    typer.echo(f"📊 Contig sizes: {outdir}/contig_sizes_sorted.tsv\n")


# -------------------------
# QC
# -------------------------
@app.command()
def qc(
    fasta: str,
    telomere: str | None = None
):
    run_qc(fasta, telomere)


# -------------------------
# INTERACTIVE MODE
# -------------------------
@app.command()
def interactive():

    typer.echo("\n🧬 Welcome to FungalFlye")
    typer.echo("Long-read fungal genome assembly assistant\n")

    reads = typer.prompt("Enter path to Nanopore reads (FASTQ/FASTQ.GZ)")
    gsize_input = typer.prompt("Estimated genome size (e.g. 40m or 40000000)")
    threads = typer.prompt("Threads", default=8)
    outdir = typer.prompt("Output directory", default="fungalflye_out")

    Path(outdir).mkdir(exist_ok=True)

    typer.echo("\n🔎 Analyzing read lengths...\n")

    total_reads, total_bases, read_n50, lengths = analyze_reads(reads, outdir)

    # convert genome size
    gsize_numeric = gsize_input.lower().replace("m", "000000")
    gsize_numeric = int(gsize_numeric)

    coverage = total_bases / gsize_numeric

    typer.echo(f"Total reads: {total_reads:,}")
    typer.echo(f"Total bases: {total_bases:,}")
    typer.echo(f"Read N50: {read_n50:,} bp")
    typer.echo(f"Estimated coverage: {coverage:.1f}×")

    suggested_cutoff = int(read_n50 * 0.7)

    typer.echo(f"\nSuggested minimum read length: {suggested_cutoff} bp")

    apply_filter = typer.confirm("Apply read filtering?", default=True)

    filtered_reads = reads

    if apply_filter:

        cutoff = typer.prompt(
            "Minimum read length",
            default=suggested_cutoff
        )

        kept, removed = preview_filter(lengths, cutoff)

        typer.echo(f"\nReads kept: {kept:,}")
        typer.echo(f"Reads removed: {removed:,}")

        confirm = typer.confirm("Continue with filtering?", default=True)

        if confirm:

            filtered_reads = f"{outdir}/filtered.fastq"

            run(
                f"seqkit seq -m {cutoff} {reads} > {filtered_reads}"
            )

    typer.echo("\n🚀 Starting assembly pipeline\n")

    assemble(
        reads=filtered_reads,
        gsize=gsize_input,
        outdir=outdir,
        threads=threads
    )

    typer.echo("\n🎉 Pipeline complete!\n")


if __name__ == "__main__":
    app()

