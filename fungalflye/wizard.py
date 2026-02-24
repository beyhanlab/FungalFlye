from pathlib import Path
import os
import typer

from .assemble import run_assembly
from .qc import run_qc, discover_telomere_motif
from .cli import analyze_reads, preview_filter

app = typer.Typer()

BANNER = r"""
🧬🐉  FungalFlye  🐉🧬
Long-read fungal genome assembly pipeline

Interactive genome assembly for fungi
"""


# ------------------------------------------------
# helpers
# ------------------------------------------------

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


# ------------------------------------------------
# telomere setup
# ------------------------------------------------

def get_telomere_setup():

    run_telomeres = typer.confirm("Run telomere analysis?", default=True)

    if not run_telomeres:
        return False, None, False

    typer.echo("\nTelomere motif:")
    typer.echo("  1) I know the motif")
    typer.echo("  2) Auto discover")

    choice = typer.prompt("Enter 1 or 2", default="2")

    if choice.strip() == "1":
        motif = typer.prompt("Enter motif sequence").strip().upper()
        return True, motif, False

    return True, None, True


# ------------------------------------------------
# resume detection
# ------------------------------------------------

def detect_existing_run(outdir):

    outdir = Path(outdir)

    final = outdir / "final.fasta"

    if final.exists():
        return "complete"

    flye = outdir / "flye" / "assembly.fasta"
    racon = outdir / "racon.fasta"

    if racon.exists():
        return "racon"

    if flye.exists():
        return "flye"

    return None


# ------------------------------------------------
# MAIN WIZARD
# ------------------------------------------------

@app.command()
def wizard():

    typer.echo(BANNER)

    while True:

        typer.echo("\nSelect mode:")
        typer.echo("  1) Full pipeline")
        typer.echo("  2) Assembly only")
        typer.echo("  3) QC only")
        typer.echo("  4) Exit")

        mode = typer.prompt("Enter choice", default="1")

        if mode == "4":
            raise typer.Exit()

        # -------------------------
        # QC ONLY
        # -------------------------

        if mode == "3":

            fasta = typer.prompt(
                "Path to assembly FASTA",
                value_proc=path_exists
            )

            run_telomeres, tel_motif, auto_tel = get_telomere_setup()

            if run_telomeres and auto_tel:
                tel_motif = discover_telomere_motif(fasta)

            run_qc(
                fasta,
                telomere=tel_motif,
                run_telomeres=run_telomeres
            )

            raise typer.Exit()

        # -------------------------
        # ASSEMBLY INPUT LOOP
        # -------------------------

        while True:

            reads = typer.prompt(
                "Path to raw reads",
                value_proc=path_exists
            )

            gsize = typer.prompt("Genome size (e.g., 40m)", default="40m")
            gsize = normalize_gsize(gsize)

            outdir = typer.prompt("Output folder", default="fungalflye_out")
            threads = typer.prompt("Threads", default=8, type=int)

            typer.echo("\nPlan:")
            typer.echo(f"Reads: {reads}")
            typer.echo(f"Genome size: {gsize}")
            typer.echo(f"Outdir: {outdir}")
            typer.echo(f"Threads: {threads}")

            typer.echo("\nProceed?")
            typer.echo("  1) Yes")
            typer.echo("  2) Edit")
            typer.echo("  3) Cancel")

            choice = typer.prompt("Enter choice", default="1")

            if choice == "1":
                break
            elif choice == "3":
                raise typer.Exit()

        Path(outdir).mkdir(exist_ok=True)

        # ------------------------------------------------
        # RESUME DETECTION
        # ------------------------------------------------

        existing = detect_existing_run(outdir)

        if existing:
            typer.echo(f"\n⚠ Existing run detected: {existing}")

            resume = typer.confirm("Resume previous run?", default=True)

            if resume and existing == "complete":
                typer.echo("Assembly already complete.")
                final_fasta = Path(outdir) / "final.fasta"
            else:
                final_fasta = None
        else:
            final_fasta = None

        # ------------------------------------------------
        # READ ANALYSIS
        # ------------------------------------------------

        typer.echo("\n🔎 Analyzing reads...\n")

        total_reads, total_bases, read_n50, lengths = analyze_reads(reads, outdir)

        typer.echo(f"Total reads: {total_reads:,}")
        typer.echo(f"Total bases: {total_bases:,}")
        typer.echo(f"Read N50: {read_n50:,} bp")

        suggested_cutoff = int(read_n50 * 0.7)

        typer.echo(f"Suggested cutoff: {suggested_cutoff} bp")

        # ------------------------------------------------
        # FILTER
        # ------------------------------------------------

        min_read_len = 0

        if typer.confirm("Apply read filtering?", default=True):

            cutoff = typer.prompt(
                "Minimum read length",
                default=suggested_cutoff
            )

            kept, removed = preview_filter(lengths, cutoff)

            typer.echo(f"Reads kept: {kept:,}")
            typer.echo(f"Reads removed: {removed:,}")

            if typer.confirm("Continue?", default=True):
                min_read_len = cutoff

        # ------------------------------------------------
        # DOWNSAMPLE
        # ------------------------------------------------

        downsample_cov = typer.prompt(
            "Downsample coverage? (0 = none)",
            default=0,
            type=int
        )

        # ------------------------------------------------
        # TELEMERES
        # ------------------------------------------------

        run_telomeres, tel_motif, auto_tel = get_telomere_setup()

        pause("Ready. Press Enter to start assembly")

        # ------------------------------------------------
        # RUN ASSEMBLY
        # ------------------------------------------------

        if final_fasta is None:

            final_fasta = run_assembly(
                reads=reads,
                genome_size=gsize,
                outdir=outdir,
                threads=threads,
                min_read_len=min_read_len,
                downsample_cov=downsample_cov,
            )

        # ------------------------------------------------
        # QC
        # ------------------------------------------------

        if mode == "1":

            pause("Assembly finished. Press Enter to run QC")

            if run_telomeres and auto_tel:
                tel_motif = discover_telomere_motif(final_fasta)

            run_qc(
                final_fasta,
                telomere=tel_motif,
                run_telomeres=run_telomeres
            )

        typer.echo("\n🎉 Pipeline complete.")
        typer.echo(f"Final assembly: {final_fasta}\n")

        raise typer.Exit()


if __name__ == "__main__":
    app()