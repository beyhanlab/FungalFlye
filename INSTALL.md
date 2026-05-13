# FunCAT Installation Guide — COMPLETE

**FunCAT — Fungal Chromosome Assembly Tool**

Developed by Jacob Durazo, Beyhan Lab, J. Craig Venter Institute

---

## One-Command Installation (RECOMMENDED)

```bash
git clone https://github.com/beyhanlab/FunCAT.git
cd FunCAT
conda env create -f environment.yml
conda activate funcat
funcat
```

**That's it.** Everything installs automatically, including:
- ✅ All assembly tools (Flye, minimap2, seqkit, filtlong)
- ✅ All polishing tools (Medaka, Racon, Polypolish, Pilon)
- ✅ All analysis tools (purge_dups, MUMmer4, BWA)
- ✅ All Python dependencies (pandas, matplotlib, biopython, typer)

---

## What Gets Installed

### Core Assembly (Required for all runs)
- **Flye** — de novo long-read assembly engine
- **minimap2** — universal read-to-contig mapping
- **seqkit** — read/contig statistics and filtering
- **filtlong** — coverage-based read downsampling
- **samtools ≥1.17** — BAM handling (MUST be ≥1.17 for Medaka)
- **bcftools** — VCF processing (Medaka requirement)
- **htslib** — bgzip/tabix VCF compression (Medaka requirement)

### Polishing & Enhancement (All included for flexibility)
- **Medaka** — ONT Nanopore consensus polishing
- **Racon** — PacBio HiFi read polishing
- **Polypolish** — Illumina hybrid polishing (recommended)
- **Pilon** — Alternative Illumina polisher
- **BWA** — Illumina read alignment

### Assembly Optimization & Analysis (All included)
- **purge_dups** — Diploid haplotig and repeat detection
- **mummer4** — Whole-genome SNP detection and dotplots

### Python Libraries (All included)
- **Biopython** — FASTA/FASTQ parsing
- **pandas** — Data handling and TSV processing
- **matplotlib** — QC plots and visualizations
- **typer** — CLI framework for all commands
- **pyabpoa** — Medaka performance optimization

---

## Manual Installation (if needed)

### Step 1: Create Conda environment
```bash
conda create -n funcat python=3.10
conda activate funcat
```

### Step 2: Install core tools from bioconda
```bash
conda install -c bioconda \
  flye minimap2 seqkit filtlong \
  'samtools>=1.17' bcftools htslib \
  bwa racon polypolish pilon \
  purge_dups mummer4
```

### Step 3: Install Python packages
```bash
conda install -c conda-forge \
  pandas>=1.5 matplotlib>=3.6 \
  biopython>=1.81 decorator>=4.0 \
  typer>=0.9
```

### Step 4: Install Medaka (AFTER bcftools/htslib on PATH)
```bash
pip install medaka pyabpoa
```

### Step 5: Install FunCAT
```bash
pip install git+https://github.com/beyhanlab/FunCAT.git
```

### Verify installation
```bash
funcat --help
```

---

## Common Issues & Fixes

### samtools version too old
If you see `samtools: invalid option -- 'a'` or similar:
```bash
conda install -c bioconda 'samtools>=1.17'
which samtools    # Should point to conda env, not system
```

### Medaka complains about missing bcftools/bgzip/tabix
These must be on PATH BEFORE you `pip install medaka`:
```bash
conda install -c bioconda bcftools htslib
pip install medaka
```

### polypolish_insert_filter not found
It's included in the polypolish conda package. Reinstall:
```bash
conda install -c bioconda polypolish
```

### Environment already exists, want to rebuild
```bash
conda env remove -n funcat
conda env create -f environment.yml
conda activate funcat
```

---

## Quick Start

```bash
# Full interactive pipeline
funcat

# Or individual commands
funcat assemble reads.fastq.gz 40m --threads 8
funcat qc assembly.fasta
funcat report assembly.fasta
funcat snps reference.fasta query.fasta
funcat dotplot reference.fasta query.fasta
funcat telo-scaffold assembly.fasta reads.fastq.gz
```

---

## System Requirements

- **OS**: macOS or Linux (Windows via WSL2)
- **CPU**: 4+ cores recommended (8+ for optimal speed)
- **RAM**: 16 GB minimum, 32+ GB recommended
- **Disk**: 200+ GB free for typical fungal genomes (50-100 Mb assemblies + reads)
- **Python**: 3.9+

---

## Citation

Durazo J, et al. FunCAT — Fungal Chromosome Assembly Tool (2025) [Manuscript in preparation]

---

## Support

For issues: https://github.com/beyhanlab/FunCAT/issues
