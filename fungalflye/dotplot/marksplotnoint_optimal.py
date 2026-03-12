import os
import csv
from Bio import SeqIO
from FastaGenome import FastaGenome
from GenomeCoord import GenomeCoord
from MUMmerTools import NucmerMap

# Define paths
genome1_path = "/Users/jdurazo/Library/CloudStorage/OneDrive-J.CraigVenterInstitute/Desktop/Projects/Sequencing/PolishingOutput/CA_polished_genomes_multiround/CA103/Round_4/CA103_polished_round_4.fasta"
genome2_path = "/Users/jdurazo/Library/CloudStorage/OneDrive-J.CraigVenterInstitute/Desktop/Projects/Sequencing/PolishingOutput/CA_polished_genomes_multiround/CA35/Round_4/CA35_polished_round_4.fasta"

delta_file = f"/Users/jdurazo/Markscripts/CA103_polished_round_4_CA35_polished_round_4.delta"
coords_file = f"/Users/jdurazo/Markscripts/CA103_polished_round_4_CA35_polished_round_4.coords"
output_csv = f"/Users/jdurazo/Markscripts/CA103_CA35_alignment_data.csv"

sorted_genome1_path = "/Users/jdurazo/Markscripts/CA103_sorted.fasta"
sorted_genome2_path = "/Users/jdurazo/Markscripts/CA35_sorted.fasta"

# Helper function to compute contig lengths and sort
def sort_fasta_by_length(fasta_path, output_path):
    contig_lengths = []
    with open(fasta_path, "rt") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            contig_lengths.append((record.id, len(record.seq)))
    
    # Sort contigs by length
    sorted_contigs = sorted(contig_lengths, key=lambda x: x[1], reverse=True)
    
    # Create a sorted FASTA file
    with open(fasta_path, "rt") as fasta_file, open(output_path, "wt") as sorted_file:
        records = list(SeqIO.parse(fasta_file, "fasta"))
        record_dict = {record.id: record for record in records}
        
        for contig_id, _ in sorted_contigs:
            if contig_id in record_dict:
                SeqIO.write(record_dict[contig_id], sorted_file, "fasta")

# Sort genomes
sort_fasta_by_length(genome1_path, sorted_genome1_path)
sort_fasta_by_length(genome2_path, sorted_genome2_path)

# Load sorted genomes
genome1 = FastaGenome(open(sorted_genome1_path, "rt"), "genome1")
genome2 = FastaGenome(open(sorted_genome2_path, "rt"), "genome2")

# Create NucmerMap
comp = NucmerMap.from_coords(coords_file, genome1, genome2, GenomeCoord(genome1), GenomeCoord(genome2))

# Debugging output
print(f"Number of pairs: {len(comp.pairs)}")
for i, pair in enumerate(comp.pairs[:10]):  # Print the first 10 pairs
    print(f"Pair {i}: {pair}")

# Export alignment data to CSV
with open(output_csv, 'w', newline='') as csvfile:
    fieldnames = ['contigA', 'startA', 'endA', 'contigB', 'startB', 'endB']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for pair in comp.pairs:
        # Ensure there are enough elements in the pair before accessing them
        if len(pair) >= 6:
            writer.writerow({
                'contigA': pair[0],
                'startA': pair[1],
                'endA': pair[2],
                'contigB': pair[3],
                'startB': pair[4],
                'endB': pair[5]
            })
        else:
            print(f"Warning: Skipping pair with insufficient elements: {pair}")











