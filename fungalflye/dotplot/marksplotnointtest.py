import os
import matplotlib.pyplot as plt
from FastaGenome import FastaGenome
from Genome import MemGeneSet, spanningLocus
from GenomeCoord import GenomeCoord
from Locus import Locus
from MUMmerTools import NucmerMap, ClickCoords

# Change directory to the full path of "Markscripts"
os.chdir("/Users/jdurazo/Markscripts")

# Define paths for genomes
genome1_path = "/Users/jdurazo/Desktop/Projects/Sequencing/CAassemblies/Allfastas/example2/CA25.fasta"
genome2_path = "/Users/jdurazo/Desktop/Projects/Sequencing/CAassemblies/Allfastas/example2/CA25test.fasta"

# Extract filenames without extensions
genome1_name = os.path.splitext(os.path.basename(genome1_path))[0]
genome2_name = os.path.splitext(os.path.basename(genome2_path))[0]

# Define paths for delta and coordinate files based on genome filenames
delta_file = f"CA25_CA25test.delta"
coords_file = f"CA25_CA25test.coords"
plots_dir = "/Users/jdurazo/Desktop/Dotplots/Newcomparisons/"

# Load genomes
genome1 = FastaGenome(open(genome1_path, "rt"), genome1_name)
genome2 = FastaGenome(open(genome2_path, "rt"), genome2_name)

# Align genomes using MUMmer
alignment_output = f"{plots_dir}/{genome1_name}_{genome2_name}"
os.system(f"nucmer -p {alignment_output} {genome1_path} {genome2_path}")

# Generate coordinate file
os.system(f"show-coords -r -c -l -d {delta_file} > {coords_file}")

# Create NucmerMap
comp = NucmerMap.from_coords(coords_file, genome1, genome2, GenomeCoord(genome1), GenomeCoord(genome2))

# Copy genome1 file to plots directory
os.system(f"cp {genome1_path} {plots_dir}")

# Plot
fig = comp.plot()
click_coords = ClickCoords(comp)
fig.canvas.mpl_connect('button_press_event', click_coords)
plt.title(f"{genome1_name}_{genome2_name}")
plt.savefig(f"{plots_dir}/{genome1_name}_{genome2_name}_dotplot.pdf")
plt.show()
