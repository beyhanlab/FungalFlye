# Load necessary libraries
library(GenomicAlignments)
library(karyoploteR)

# Read BAM file
bam_file <- "/Users/jdurazo/Desktop/Projects/Sequencing/CAassemblies/Allfastas/alignmentresults/alignment_results.sorted.bam"
CDC_8Osort <- readGAlignments(bam_file)

# Generate contig labels
seqinfo <- seqinfo(CDC_8Osort)
contig_labels <- paste("Contig", seqinfo@seqnames)

# Sort contigs by length
sorted_seqnames <- seqinfo@seqnames[order(-seqinfo@seqlengths)]
CDC_8Ogenome <- GRanges(seqnames = sorted_seqnames, 
                        ranges = IRanges(start = 1, end = seqinfo@seqlengths[order(-seqinfo@seqlengths)]))

# Define plot parameters
pp <- getDefaultPlotParams(plot.type = 2)
pp$leftmargin <- 0.15

# Plot karyotype
kp <- plotKaryotype(genome = CDC_8Ogenome, plot.type = 2, plot.params = pp, chromosomes = "all", labels=contig_labels)

# Calculate coverage
CDC_8Ocoverage <- coverage(CDC_8Osort)

# Plot coverage
kpPlotCoverage(kp, data = CDC_8Ocoverage, col = "#f60909", ylab="Coverage")

# Save plot image
png("~/Candida_auris/AurisMappingstoReference/AurisWorkspace.png")
dev.off()

