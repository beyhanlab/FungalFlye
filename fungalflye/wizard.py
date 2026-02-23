from pathlib import Path
import os
import typer

from .assemble import run_assembly
from .qc import run_qc, discover_telomere_motif
from .cli import analyze_reads, preview_filter   # <-- bring in the good stuff

app = typer.Typer()

BANNER = r"""
🧬🐉  FungalFlye  🐉🧬
Long-read fungal genome assembly pipeline

Modes:
  1) Full assembly + QC
  2) QC only (use existing FASTA)

(Press Ctrl+C anytime to exit.)
"""


# -------------------------
# helpers
# -------------------------

def pause(msg="Press Enter to continue..."):
    typer.echo("")
    typer.prompt(msg, default="", show_default=False)
    typer.echo("")


def path_exists(p: str) -> str:
    p2 = os.path.expanduser(p.strip())
    if not Path(p2).exists():
        raise typer.BadParameter(f"Path not found: {p2}")
    return p2


def normalize_gsize(g: str) -> str:
    g = g.strip().lower()
    if g.isdigit():
        val = int(g)
        if val < 1000:
            return f"{val}m"
        return str(val)
    return g


# -------------------------
# telomere setup
# -------------------------

def get_telomere_setup():

    run_telomeres = typer.confirm("Run telomere analysis?", default=True)

    if not run_telomeres:
        return False, None, False

    typer.echo("\nDo you know the telomere motif?")
    typer.echo("  1) Yes — I will enter it")
    typer.echo("  2) No — discover automatically")

    choice = typer.prompt("Enter 1 or 2", default="2")

    if choice.strip() == "1":
        motif = typer.prompt(
            "Enter telomere motif sequence (e.g., TTAGGG)"
        ).strip().upper()
        return True, motif, False

    return True, None, True


# -------------------------
# wizard command
# -------------------------

@app.command()
def wizard():

    typer.echo(BANNER)

    typer.echo("\nSelect mode:")
    typer.echo("  1) Full assembly pipeline")
    typer.echo("  2) QC only (use existing FASTA)")

    mode = typer.prompt("Enter 1 or 2", default="1")

    # Telomere setup
    run_telomeres, tel_motif, auto_tel = get_telomere_setup()

    # =====================================================
    # MODE 2 — QC ONLY
    # =====================================================

    if mode.strip() == "2":

        fasta = typer.prompt(
            "Path to assembly FASTA",
            value_proc=path_exists
        )

        typer.echo("\n🧾 QC Plan:")
        typer.echo(f"  FASTA: {fasta}")

        if run_telomeres:
            if tel_motif:
                typer.echo(f"  Telomeres: ON ({tel_motif})")
            else:
                typer.echo("  Telomeres: ON (auto discovery)")
        else:
            typer.echo("  Telomeres: OFF")

        pause("Press Enter to start QC")

        if run_telomeres and auto_tel:
            tel_motif = discover_telomere_motif(fasta)

        run_qc(
            fasta,
            telomere=tel_motif,
            run_telomeres=run_telomeres
        )

        typer.echo("\n✅ QC complete.")
        raise typer.Exit()

    # =====================================================
    # MODE 1 — FULL PIPELINE
    # =====================================================

    reads = typer.prompt(
        "Path to raw reads (FASTQ/FASTQ.GZ)",
        value_proc=path_exists
    )

    gsize = typer.prompt("Genome size (e.g., 40m)", default="40m")
    gsize = normalize_gsize(gsize)

    outdir = typer.prompt("Output folder", default="fungalflye_out")
    threads = typer.prompt("Threads", default=8, type=int)

    Path(outdir).mkdir(exist_ok=True)

    # -------------------------
    # READ ANALYSIS
    # -------------------------

    typer.echo("\n🔎 Analyzing read lengths...\n")

    total_reads, total_bases, read_n50, lengths = analyze_reads(reads, outdir)

    # convert genome size numeric
    gsize_numeric = gsize.lower().replace("m", "000000")
    gsize_numeric = int(gsize_numeric)

    coverage = total_bases / gsize_numeric

    typer.echo(f"Total reads: {total_reads:,}")
    typer.echo(f"Total bases: {total_bases:,}")
    typer.echo(f"Read N50: {read_n50:,} bp")
    typer.echo(f"Estimated coverage: {coverage:.1f}×")

    suggested_cutoff = int(read_n50 * 0.7)

    typer.echo(f"\nSuggested minimum read length: {suggested_cutoff} bp")

    # -------------------------
    # FILTERING
    # -------------------------

    apply_filter = typer.confirm("Apply read filtering?", default=True)

    filtered_reads = reads
    min_read_len = suggested_cutoff

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

            typer.echo("\n🧬 Filtering reads\n")

            os.system(
                f"seqkit seq -m {cutoff} {reads} > {filtered_reads}"
            )

            min_read_len = cutoff

    # -------------------------
    # DOWNSAMPLE OPTION
    # -------------------------

    downsample_cov = typer.prompt(
        "Downsample coverage? (0 = no downsample)",
        default=0,
        type=int
    )

    typer.echo("\n🧾 Plan:")
    typer.echo(f"  Reads: {filtered_reads}")
    typer.echo(f"  Genome size: {gsize}")
    typer.echo(f"  Outdir: {outdir}")
    typer.echo(f"  Threads: {threads}")
    typer.echo(f"  Min read len: {min_read_len}")

    pause("Ready. Press Enter to start assembly")

    # -------------------------
    # RUN ASSEMBLY
    # -------------------------

    run_assembly(
        reads=filtered_reads,
        genome_size=gsize,
        outdir=outdir,
        threads=threads,
        min_read_len=min_read_len,
        downsample_cov=downsample_cov,
    )

    pause("Assembly finished. Press Enter to run QC")

    final_fasta = Path(outdir) / "final.fasta"

    if not final_fasta.exists():
        typer.echo(f"❌ ERROR: Final assembly not found: {final_fasta}")
        raise typer.Exit(code=1)

    if run_telomeres and auto_tel:
        tel_motif = discover_telomere_motif(str(final_fasta))

    run_qc(
        str(final_fasta),
        telomere=tel_motif,
        run_telomeres=run_telomeres
    )

    typer.echo("\n✅ All done.")
    typer.echo(f"Final assembly: {final_fasta}\n")


if __name__ == "__main__":
    app()