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
"""


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


def get_telomere_setup():

    run_telomeres = typer.confirm("Run telomere analysis?", default=True)

    if not run_telomeres:
        return False, None, False

    typer.echo("\nDo you know the telomere motif?")
    typer.echo("  1) Yes")
    typer.echo("  2) Auto discover")

    choice = typer.prompt("Enter 1 or 2", default="2")

    if choice.strip() == "1":
        motif = typer.prompt(
            "Enter telomere motif sequence"
        ).strip().upper()
        return True, motif, False

    return True, None, True


@app.command()
def wizard():

    typer.echo(BANNER)

    # -------------------------
    # INPUT LOOP (editable)
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

        typer.echo("\n🧾 Plan:")
        typer.echo(f"Reads: {reads}")
        typer.echo(f"Genome size: {gsize}")
        typer.echo(f"Outdir: {outdir}")
        typer.echo(f"Threads: {threads}")

        typer.echo("\nProceed?")
        typer.echo("  1) Yes")
        typer.echo("  2) Edit parameters")
        typer.echo("  3) Cancel")

        choice = typer.prompt("Enter 1/2/3", default="1")

        if choice == "1":
            break
        elif choice == "3":
            raise typer.Exit()

    Path(outdir).mkdir(exist_ok=True)

    # -------------------------
    # READ ANALYSIS
    # -------------------------

    typer.echo("\n🔎 Analyzing read lengths...\n")

    total_reads, total_bases, read_n50, lengths = analyze_reads(reads, outdir)

    typer.echo(f"Total reads: {total_reads:,}")
    typer.echo(f"Total bases: {total_bases:,}")
    typer.echo(f"Read N50: {read_n50:,} bp")

    suggested_cutoff = int(read_n50 * 0.7)

    typer.echo(f"\nSuggested minimum read length: {suggested_cutoff} bp")

    # -------------------------
    # FILTER OPTION
    # -------------------------

    apply_filter = typer.confirm("Apply read filtering?", default=True)

    min_read_len = 0

    if apply_filter:

        cutoff = typer.prompt(
            "Minimum read length",
            default=suggested_cutoff
        )

        kept, removed = preview_filter(lengths, cutoff)

        typer.echo(f"\nReads kept: {kept:,}")
        typer.echo(f"Reads removed: {removed:,}")

        confirm = typer.confirm("Continue?", default=True)

        if confirm:
            min_read_len = cutoff

    # -------------------------
    # DOWNSAMPLE
    # -------------------------

    downsample_cov = typer.prompt(
        "Downsample coverage? (0 = none)",
        default=0,
        type=int
    )

    run_telomeres, tel_motif, auto_tel = get_telomere_setup()

    pause("Ready. Press Enter to start assembly")

    final_fasta = run_assembly(
        reads=reads,
        genome_size=gsize,
        outdir=outdir,
        threads=threads,
        min_read_len=min_read_len,
        downsample_cov=downsample_cov,
    )

    pause("Assembly finished. Press Enter to run QC")

    if run_telomeres and auto_tel:
        tel_motif = discover_telomere_motif(final_fasta)

    run_qc(
        final_fasta,
        telomere=tel_motif,
        run_telomeres=run_telomeres
    )

    typer.echo("\n✅ All done.")


if __name__ == "__main__":
    app()