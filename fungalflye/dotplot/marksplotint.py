import os
import matplotlib.pyplot as plt
from FastaGenome import FastaGenome
from Genome import MemGeneSet, spanningLocus
from GenomeCoord import GenomeCoord
from Locus import Locus
from MUMmerTools import NucmerMap, ClickCoords

# Change directory to the full path of "Markscripts"
os.chdir("/Users/jdurazo/Markscripts")

# Enable matplotlib interactive mode
plt.ion()

# Define paths
genome1_path = "/Users/jdurazo/Desktop/Projects/Sequencing/CAassemblies/Allfastas/CA56.fasta"
genome2_path = "/Users/jdurazo/Desktop/Projects/Sequencing/CAassemblies/Allfastas/CA58.fasta"
delta_file = "CA56_CA58.delta"
coords_file = "CA56_CA58.coords"
plots_dir = "/Users/jdurazo/Desktop/Projects/Sequencing/CAassemblies/Allfastas/Mummerplots/"

# Load genomes
genome1 = FastaGenome(open(genome1_path, "rt"), "CA25")
genome2 = FastaGenome(open(genome2_path, "rt"), "CA34")

# Run show-coords command
os.system(f"show-coords -r -c -l -d {delta_file} > {coords_file}")

# Create NucmerMap
comp = NucmerMap.from_coords(coords_file, genome1, genome2, GenomeCoord(genome1), GenomeCoord(genome2))

# Copy delta file
os.system(f"cp {genome1_path} {plots_dir}")

# Display head of coords file
with open(coords_file, "r") as f:
    print(f.read())

# Plot
fig = comp.plot()
click_coords = ClickCoords(comp)
fig.canvas.mpl_connect('button_press_event', click_coords)
plt.title("CA25_CA34")
plt.savefig("mydotplot.pdf")
plt.show()
