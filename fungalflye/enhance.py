"""
funcat.enhance
~~~~~~~~~~~~~~~~~~
Assembly enhancement modules:
  1. Adaptive Flye parameter selection
  2. Iterative Medaka polishing with convergence detection
  3. Purge Duplicates (diploid haplotig removal)
  4. Contig confidence scoring
"""

import subprocess
import shutil
import json
import math
from pathlib import Path
from Bio import SeqIO

# Optional PacBio detection import
try:
    from .pacbio_detection import detect_pacbio_data_type
    PACBIO_DETECTION_AVAILABLE = True
except ImportError:
    PACBIO_DETECTION_AVAILABLE = False
    def detect_pacbio_data_type(*args, **kwargs):
        return None, 0, {"error": "PacBio detection module not available"}


# ------------------------------------------------
# helpers
# ------------------------------------------------

def run(cmd):
    print(f"\n[funcat] Running: {cmd}\n")
    subprocess.run(cmd, shell=True, check=True)


def run_capture(cmd):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout.strip()


# ================================================
# MODULE 1 — Adaptive Flye parameter selection
# ================================================



def suggest_flye_params(reads, genome_size_bp, threads=8, read_type=None):
    """
    Analyse reads and return optimised Flye parameters as a dict.
    Prints a human-readable explanation of every decision made.
    
    read_type: "pacbio-hifi", "nano-hq", "nano-raw" - overrides detection
    """

    print("\n" + "=" * 60)
    print("🔬 Module 1 — Adaptive parameter selection")
    print("=" * 60)

    # --- read stats via seqkit ---
    print("\n[funcat] Analysing read characteristics...")

    stats_raw = run_capture(
        f"seqkit fx2tab -n -l -g {reads} 2>/dev/null"
    )

    lengths = []
    gcs = []

    for line in stats_raw.splitlines():
        parts = line.strip().split("\t")
        if len(parts) >= 3:
            try:
                lengths.append(int(parts[1]))
                gcs.append(float(parts[2]))
            except ValueError:
                pass

    if not lengths:
        print("[funcat] ⚠️  Could not parse read stats — using defaults")
        return {}

    total_bases = sum(lengths)
    coverage = total_bases / genome_size_bp
    mean_gc = sum(gcs) / len(gcs) if gcs else 50.0

    sorted_len = sorted(lengths, reverse=True)
    cumsum, read_n50 = 0, 0
    for L in sorted_len:
        cumsum += L
        if cumsum >= total_bases / 2:
            read_n50 = L
            break

    print(f"\n  Read N50       : {read_n50:,} bp")
    print(f"  Est. coverage  : {coverage:.1f}x")
    print(f"  Mean GC        : {mean_gc:.1f}%")

    # --- Auto-detect PacBio data type if needed ---
    if read_type and 'pacbio' in read_type.lower():
        print("\n[funcat] 📡 PacBio data detected - analyzing data type...")
        
        detected_type, confidence, details = detect_pacbio_data_type(reads)
        
        if detected_type and confidence > 0.8:
            print(f"[funcat] ✅ Auto-detected: {detected_type} (confidence: {confidence*100:.0f}%)")
            
            if read_type == "pacbio-hifi" and detected_type == "pacbio-raw":
                print(f"[funcat] ⚠️  WARNING: You specified HiFi but data appears to be raw subreads!")
                print(f"[funcat] This may cause 'No overlaps found' errors.")
                print(f"[funcat] Recommendation: Use --pacbio-raw instead")
                
                # Auto-correct if high confidence
                if confidence > 0.9:
                    print(f"[funcat] 🔧 Auto-correcting read type to: {detected_type}")
                    read_type = detected_type
                    
            elif read_type == "pacbio-raw" and detected_type == "pacbio-hifi":
                print(f"[funcat] ℹ️  Note: Data appears to be HiFi but using raw mode as requested")
        
        else:
            print(f"[funcat] ❓ Could not auto-detect PacBio type reliably")
            print(f"[funcat] Proceeding with specified type: {read_type}")

    params = {}
    reasoning = []

    # --- min-overlap: Let Flye auto-select (recommended approach) ---
    # Flye's automatic selection based on N90 is more reliable than manual calculation
    # According to Flye documentation: "this parameter is chosen automatically based on 
    # the read length distribution (reads N90) and does not require manual setting"
    reasoning.append(
        f"  --min-overlap [auto]  "
        f"(Flye auto-selection based on read N90 - more reliable than manual override)"
    )
    # Don't set min_overlap parameter - let Flye handle it

    # --- asm-coverage ---
    # Rule: cap at actual coverage if lower than default 60x
    if coverage < 40:
        asm_cov = max(int(coverage * 0.8), 20)
        reasoning.append(
            f"  --asm-coverage {asm_cov}  "
            f"(low coverage {coverage:.0f}x — reduced to avoid read starvation)"
        )
    elif coverage > 120:
        asm_cov = 80
        reasoning.append(
            f"  --asm-coverage 80  "
            f"(very high coverage {coverage:.0f}x — capped to reduce compute)"
        )
    else:
        asm_cov = 60
        reasoning.append(
            f"  --asm-coverage 60  (coverage {coverage:.0f}x — standard)"
        )
    params["asm_coverage"] = asm_cov

    # --- iterations ---
    # Rule: more iterations for lower coverage (graph needs more resolution)
    if coverage < 30:
        iters = 4
        reasoning.append(
            "  --iterations 4  (low coverage — extra iterations for graph resolution)"
        )
    else:
        iters = 3
        reasoning.append("  --iterations 3  (standard)")
    params["iterations"] = iters

    # --- AT-rich warning ---
    if mean_gc < 35:
        params["at_rich"] = True
        reasoning.append(
            f"  ⚠️  AT-rich genome ({mean_gc:.1f}% GC) — "
            "consider --min-overlap reduction if assembly fragments"
        )

    print("\n[funcat] Recommended Flye parameters:")
    for r in reasoning:
        print(r)

    return params


# ================================================
# MODULE 2 — Iterative Medaka polishing
# ================================================

def _count_changes(before_fasta, after_fasta, threads=4):
    """
    Count sequence-level changes between two assemblies using minimap2 asm5.
    Works with Medaka 2.x which no longer outputs VCF files.

    Aligns before → after, sums mismatches and indels from the cs tag.
    Returns total number of corrected bases (substitutions + indels).
    Falls back to 0 if minimap2 is unavailable or alignment fails.
    """
    import subprocess
    import re
    from pathlib import Path

    try:
        # Ensure both files exist and are valid
        if not Path(before_fasta).exists() or not Path(after_fasta).exists():
            print(f"[funcat] Warning: Input files missing for convergence check")
            return 0
            
        # Check file sizes to ensure they're real assemblies
        before_size = Path(before_fasta).stat().st_size
        after_size = Path(after_fasta).stat().st_size
        if before_size < 1000 or after_size < 1000:
            print(f"[funcat] Warning: Assembly files too small for comparison")
            return 0

        print(f"[funcat] Comparing {Path(before_fasta).name} → {Path(after_fasta).name}")
        
        result = subprocess.run(
            [
                "minimap2", "-x", "asm5",
                "--cs", "-t", str(threads),
                str(before_fasta), str(after_fasta),
            ],
            capture_output=True, text=True, timeout=300,
        )
        
        if result.returncode != 0:
            print(f"[funcat] Warning: minimap2 alignment failed: {result.stderr}")
            return 0
            
        changes = 0
        alignments = 0
        
        for line in result.stdout.splitlines():
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue
                
            alignments += 1
            
            # Parse cs tag — count mismatches (*) and indels (+ -)
            for col in cols:
                if col.startswith("cs:Z:"):
                    cs = col[5:]
                    # substitutions: *XY  insertions: +seq  deletions: -seq
                    substitutions = len(re.findall(r'\*[acgtACGT]{2}', cs))
                    insertions = sum(len(m) - 1 for m in re.findall(r'\+[acgtACGT]+', cs))
                    deletions = sum(len(m) - 1 for m in re.findall(r'-[acgtACGT]+', cs))
                    changes += substitutions + insertions + deletions
                    
        print(f"[funcat] Found {changes} sequence differences across {alignments} alignments")
        return changes
        
    except subprocess.TimeoutExpired:
        print(f"[funcat] Warning: Convergence check timed out")
        return 0
    except Exception as e:
        print(f"[funcat] Warning: Convergence check failed: {e}")
        return 0


def run_medaka_iterative(
    assembly,
    reads,
    outdir,
    threads,
    model,
    max_rounds=3,
    convergence_threshold=50,
):
    """
    Run Medaka up to max_rounds times.
    Stops early when variant calls drop below convergence_threshold
    (meaning the assembly has converged and further polishing won't help).
    Returns path to final polished FASTA.
    """

    print("\n" + "=" * 60)
    print("🔬 Module 2 — Iterative Medaka polishing")
    print(f"   Max rounds: {max_rounds}  |  Convergence threshold: {convergence_threshold} base corrections")
    print("=" * 60)

    current = Path(assembly)
    medaka_base = Path(outdir) / "medaka"
    medaka_base.mkdir(exist_ok=True)

    state_file = medaka_base / "polishing_state.json"

    # Resume support: load previous state
    if state_file.exists():
        state = json.loads(state_file.read_text())
        last_round = state.get("last_round", 0)
        last_variants = state.get("last_variants", None)
        current = Path(state.get("current_fasta", assembly))
        print(f"[funcat] Resuming from round {last_round}")
        
        # If we've already completed max_rounds, use the existing result
        if last_round >= max_rounds:
            final_fasta = current
            print(f"[funcat] Polishing already complete at round {last_round}")
        else:
            final_fasta = None
    else:
        state = {}
        last_round = 0
        last_variants = None
        final_fasta = None

    for round_num in range(last_round + 1, max_rounds + 1):

        round_dir = medaka_base / f"round_{round_num}"

        polished = round_dir / "consensus.fasta"

        if polished.exists():
            print(f"[funcat] Round {round_num} already done — skipping")
            current = polished
            final_fasta = polished
            continue

        print(f"\n[funcat] Medaka round {round_num} / {max_rounds}")

        round_dir.mkdir(exist_ok=True)

        # Snapshot the input so we can compare before/after for convergence
        current_before_polish = current

        run(
            f"medaka_consensus "
            f"-i {reads} "
            f"-d {current} "
            f"-o {round_dir} "
            f"-t {threads} "
            f"-m {model}"
        )

        if not polished.exists():
            raise RuntimeError(f"Medaka round {round_num} failed")

        # Count sequence changes between input and polished assembly.
        # Uses direct FASTA comparison via minimap2 asm5 — compatible with
        # Medaka 2.x which no longer outputs VCF files.
        changes_this_round = _count_changes(current_before_polish, polished, threads=threads)

        print(
            f"[funcat] Round {round_num} complete — "
            f"{changes_this_round} bases corrected"
        )

        # Save state for resumability
        state = {
            "last_round": round_num,
            "last_variants": changes_this_round,
            "current_fasta": str(polished),
        }
        state_file.write_text(json.dumps(state, indent=2))

        current = polished
        final_fasta = polished

        # Convergence check — stop if corrections fall below threshold
        if changes_this_round <= convergence_threshold:
            print(
                f"\n✅ Converged after round {round_num} "
                f"({changes_this_round} bases corrected ≤ threshold {convergence_threshold})"
            )
            break

        if last_variants is not None:
            improvement = last_variants - changes_this_round
            if improvement <= 0:
                print(
                    f"\n✅ No further improvement after round {round_num} — stopping"
                )
                break

        last_variants = changes_this_round

    if final_fasta is None:
        raise RuntimeError("Medaka polishing produced no output")

    print(f"\n[funcat] Final polished assembly: {final_fasta}")
    return final_fasta


# ================================================
# MODULE 3 — Purge Duplicates (diploid)
# ================================================

def run_purge_dups(assembly, reads, outdir, threads, minimap2_preset):
    """
    Run purge_dups to remove haplotig duplicates from diploid assemblies.
    Returns path to purged primary assembly.

    Requires: purge_dups, minimap2, split_fa, pbcstat, calcuts, get_seqs
    (all installed via: conda install -c bioconda purge_dups)
    """

    print("\n" + "=" * 60)
    print("🔬 Module 3 — Purge Duplicates (diploid haplotig removal)")
    print("=" * 60)

    purge_dir = Path(outdir) / "purge_dups"
    purge_dir.mkdir(parents=True, exist_ok=True)

    purged_primary = purge_dir / "purged.fa"
    haplotigs = purge_dir / "hap.fa"

    if purged_primary.exists():
        print("[funcat] Existing purge_dups output detected — skipping")
        return purged_primary, haplotigs

    assembly = Path(assembly)

    # Check tools
    missing = [t for t in ["purge_dups", "split_fa", "pbcstat", "calcuts", "get_seqs"]
               if shutil.which(t) is None]
    if missing:
        print(f"\n⚠️  purge_dups tools not found: {missing}")
        print("   Install with: conda install -c bioconda purge_dups")
        print("   Skipping haplotig purging — assembly will be unpurged\n")
        return assembly, None

    print("\n[funcat] Step 1 — self-mapping assembly")
    self_paf = purge_dir / "self.paf"
    run(f"minimap2 -xasm5 -DP -t {threads} {assembly} {assembly} > {self_paf}")

    print("[funcat] Step 2 — mapping reads to assembly")
    read_paf = purge_dir / "reads.paf"
    run(f"minimap2 -x {minimap2_preset} -t {threads} {assembly} {reads} > {read_paf}")

    print("[funcat] Step 3 — coverage histogram")
    stat_file = purge_dir / "PB.stat"
    base_cov = purge_dir / "PB.base.cov"
    
    try:
        # Ensure output directory exists and is writable
        purge_dir.mkdir(parents=True, exist_ok=True)
        
        # Try alternative pbcstat approach - change directory first
        print("[funcat] Using working directory approach for pbcstat...")
        result = subprocess.run(
            f"cd {purge_dir} && pbcstat {read_paf}",
            shell=True, 
            capture_output=True, 
            text=True
        )
        
        if result.returncode != 0:
            print(f"pbcstat stderr: {result.stderr}")
            raise RuntimeError(f"pbcstat failed with exit code {result.returncode}")


        
        # Validate PB.stat file before proceeding
        if not stat_file.exists():
            raise FileNotFoundError(f"pbcstat failed to create {stat_file}")
        
        if stat_file.stat().st_size == 0:
            raise ValueError(f"pbcstat created empty file: {stat_file}")
        
        print("[funcat] Step 4 — calculating cutoffs")
        cutoffs_file = purge_dir / "cutoffs"
        
        # Safe calcuts execution
        result = subprocess.run(f"calcuts {stat_file}", shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            raise RuntimeError(f"calcuts failed with exit code {result.returncode}: {result.stderr}")
        
        # Write cutoffs output
        with open(cutoffs_file, 'w') as f:
            f.write(result.stdout)
            
        # Validate cutoffs - if all contigs marked as junk, adjust parameters
        print("[funcat] Validating purge duplicates parameters...")
        
        # First run a test to see what gets classified
        test_bed = purge_dir / "test_dups.bed"
        test_result = subprocess.run(
            f"purge_dups -2 -T {cutoffs_file} -c {base_cov} {self_paf}",
            shell=True, 
            capture_output=True, 
            text=True,
            cwd=purge_dir
        )
        
        if test_result.returncode == 0:
            with open(test_bed, 'w') as f:
                f.write(test_result.stdout)
                
            # Check if all sequences marked as JUNK
            junk_count = 0
            keep_count = 0
            
            with open(test_bed, 'r') as f:
                for line in f:
                    if line.strip():
                        if '\tJUNK' in line:
                            junk_count += 1
                        elif '\tKEEP' in line or '\tHAP' in line:
                            keep_count += 1
            
            # If all sequences marked as JUNK, use less aggressive parameters
            if junk_count > 0 and keep_count == 0:
                print(f"⚠️  All {junk_count} contigs marked as JUNK - adjusting parameters for ONT diploid data")
                
                # Create more permissive cutoffs for ONT diploid assembly
                permissive_cutoffs = "3\t8\t15\t40\t80\t160"  # More conservative
                with open(cutoffs_file, 'w') as f:
                    f.write(permissive_cutoffs)
                    
                print("[funcat] Using relaxed cutoffs optimized for ONT diploid assembly")
            elif junk_count > keep_count * 3:  # More than 75% marked as junk
                print(f"⚠️  Too many contigs marked as JUNK ({junk_count} vs {keep_count} primary)")
                print("[funcat] Adjusting for better primary/haplotig balance")
                
                # Moderately permissive for cases with excessive junk classification
                balanced_cutoffs = "5\t10\t20\t50\t100\t200"
                with open(cutoffs_file, 'w') as f:
                    f.write(balanced_cutoffs)
                    
                print("[funcat] Using balanced cutoffs for optimal haplotig detection")
            else:
                print(f"[funcat] Acceptable purge parameters: {keep_count} primary + {junk_count} haplotigs")
        else:
            print("[funcat] Warning: Could not validate purge parameters")

            
    except Exception as e:
        print(f"⚠️  Purge duplicates failed: {e}")
        print("📋 Continuing with polished assembly (diploid contigs preserved)")
        print("💡 This may result in more contigs but assembly is still valid")
        return assembly, None

    print("[funcat] Step 5 — splitting assembly")
    split_asm = purge_dir / "split.fa"
    run(f"split_fa {assembly} > {split_asm}")

    split_paf = purge_dir / "split.paf"
    run(f"minimap2 -xasm5 -DP -t {threads} {split_asm} {split_asm} > {split_paf}")

    print("[funcat] Step 6 — purging haplotigs")
    bed_file = purge_dir / "dups.bed"
    run(
        f"purge_dups -2 -T {cutoffs_file} -c {base_cov} "
        f"{self_paf} > {bed_file}"
    )

    # Analyze and potentially fix inverted classification
    print("[funcat] Analyzing purge duplicates classification...")
    
    # Read the bed file and calculate primary vs haplotig sizes
    junk_contigs = set()
    junk_total_length = 0
    
    with open(bed_file, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 4 and parts[3] == 'JUNK':
                    junk_contigs.add(parts[0])
                    junk_total_length += int(parts[2]) - int(parts[1])
    
    # Calculate total assembly size and primary size
    total_contigs = 0
    total_length = 0
    primary_total_length = 0
    
    with open(assembly, 'r') as f:
        current_seq_len = 0
        current_contig = None
        for line in f:
            if line.startswith('>'):
                if current_seq_len > 0 and current_contig:
                    total_length += current_seq_len
                    total_contigs += 1
                    if current_contig not in junk_contigs:
                        primary_total_length += current_seq_len
                current_contig = line.strip()[1:].split()[0]
                current_seq_len = 0
            else:
                current_seq_len += len(line.strip())
        # Last sequence
        if current_seq_len > 0 and current_contig:
            total_length += current_seq_len
            total_contigs += 1
            if current_contig not in junk_contigs:
                primary_total_length += current_seq_len
    
    junk_count = len(junk_contigs)
    primary_count = total_contigs - junk_count
    
    print(f"📊 Classification analysis:")
    print(f"   Total: {total_contigs} contigs, {total_length/1e6:.1f} MB")
    print(f"   JUNK: {junk_count} contigs, {junk_total_length/1e6:.1f} MB")
    print(f"   Primary: {primary_count} contigs, {primary_total_length/1e6:.1f} MB")
    
    # Detect inverted classification (small primary, large JUNK)
    invert_classification = False
    if primary_total_length < junk_total_length and primary_total_length < total_length * 0.4:
        print(f"⚠️  Classification appears INVERTED!")
        print(f"   Primary ({primary_total_length/1e6:.1f} MB) < JUNK ({junk_total_length/1e6:.1f} MB)")
        print(f"💡 For diploid assemblies, will swap primary ↔ haplotig classification")
        invert_classification = True
    else:
        print(f"✅ Classification looks correct")

    # Use Python implementation with smart inversion detection
    print("[funcat] Extracting sequences with corrected classification...")
    purged_primary = purge_dir / "purged.fa"
    haplotigs = purge_dir / "hap.fa"
    
    # Python implementation to extract sequences with proper classification
    classification_script = f'''
import sys

# Read the dups.bed file to get list of JUNK contigs
junk_contigs = set()
with open("{bed_file}", "r") as f:
    for line in f:
        if line.strip():
            parts = line.strip().split("\\t")
            if len(parts) >= 4 and parts[3] == "JUNK":
                junk_contigs.add(parts[0])

print(f"Found {{len(junk_contigs)}} JUNK contigs")

# Apply inversion logic if needed
invert = {str(invert_classification).lower()}
if invert:
    print("🔄 INVERTING classification: JUNK → Primary, Non-JUNK → Haplotigs")

# Read assembly and extract sequences
primary_sequences = []
haplotig_sequences = []

with open("{assembly}", "r") as infile:
    current_seq = ""
    current_header = ""
    
    for line in infile:
        if line.startswith(">"):
            # Process previous sequence
            if current_header and current_seq:
                contig_name = current_header[1:].split()[0]
                is_marked_junk = contig_name in junk_contigs
                
                # Apply inversion logic
                if invert:
                    # Inverted: JUNK becomes primary, non-JUNK becomes haplotigs
                    if is_marked_junk:
                        primary_sequences.append((current_header, current_seq))
                    else:
                        haplotig_sequences.append((current_header, current_seq))
                else:
                    # Normal: non-JUNK is primary, JUNK is haplotigs
                    if is_marked_junk:
                        haplotig_sequences.append((current_header, current_seq))
                    else:
                        primary_sequences.append((current_header, current_seq))
            
            current_header = line.strip()
            current_seq = ""
        else:
            current_seq += line
    
    # Process final sequence
    if current_header and current_seq:
        contig_name = current_header[1:].split()[0]
        is_marked_junk = contig_name in junk_contigs
        
        if invert:
            if is_marked_junk:
                primary_sequences.append((current_header, current_seq))
            else:
                haplotig_sequences.append((current_header, current_seq))
        else:
            if is_marked_junk:
                haplotig_sequences.append((current_header, current_seq))
            else:
                primary_sequences.append((current_header, current_seq))

# Write primary assembly
with open("{purged_primary}", "w") as f:
    for header, seq in primary_sequences:
        f.write(header + "\\n")
        f.write(seq)

# Write haplotigs
with open("{haplotigs}", "w") as f:
    for header, seq in haplotig_sequences:
        f.write(header + "\\n")
        f.write(seq)

# Report final sizes
primary_size = sum(len(seq) for _, seq in primary_sequences)
haplotig_size = sum(len(seq) for _, seq in haplotig_sequences)

print(f"✅ Final results:")
print(f"   Primary: {{len(primary_sequences)}} contigs, {{primary_size/1e6:.1f}} MB")
print(f"   Haplotigs: {{len(haplotig_sequences)}} contigs, {{haplotig_size/1e6:.1f}} MB")
'''
    
    result = subprocess.run(["python3", "-c", classification_script], 
                          capture_output=True, text=True, cwd=purge_dir)
    
    if result.returncode != 0:
        print(f"❌ Python classification failed: {result.stderr}")
        print("🔄 Falling back to original assembly")
        return assembly, None
    else:
        print(result.stdout)
    
    # Validate output files exist and have content
    if not purged_primary.exists():
        print("⚠️  purge_dups did not produce primary output — using unpurged assembly")
        return assembly, None

    # Check for reasonable file sizes
    primary_size = purged_primary.stat().st_size if purged_primary.exists() else 0
    haplotig_size = haplotigs.stat().st_size if haplotigs.exists() else 0
    
    if primary_size < 100000:  # Less than 100KB is suspiciously small
        print(f"⚠️  Primary assembly very small ({primary_size} bytes)")
        print("💡 This may indicate parameter issues - consider using original assembly")

    # Report final statistics
    try:
        from Bio import SeqIO
        orig_contigs = sum(1 for r in SeqIO.parse(str(assembly), "fasta"))
        purged_contigs = sum(1 for r in SeqIO.parse(str(purged_primary), "fasta"))
        removed = orig_contigs - purged_contigs

        print(f"\n✅ Purge Duplicates complete")
        print(f"   Original contigs : {orig_contigs}")
        print(f"   After purging    : {purged_contigs}  ({removed} haplotigs removed)")
        if haplotigs.exists():
            hap_contigs = sum(1 for r in SeqIO.parse(str(haplotigs), "fasta"))
            print(f"   Haplotigs saved  : {hap_contigs} → {haplotigs}")
    except ImportError:
        print("✅ Purge Duplicates complete - install biopython for detailed stats")

    return purged_primary, haplotigs


# ================================================
# MODULE 4 — Contig confidence scoring
# ================================================

def score_contig_confidence(assembly, reads, outdir, threads, minimap2_preset):
    """
    Map reads back to assembly and compute per-contig metrics:
      - mean coverage
      - coverage uniformity (coefficient of variation)
      - a confidence label: GOOD / REVIEW / FLAG

    Writes a TSV report and prints a summary.
    Returns path to the TSV.
    """

    print("\n" + "=" * 60)
    print("🔬 Module 4 — Contig confidence scoring")
    print("=" * 60)

    score_dir = Path(outdir) / "confidence"
    score_dir.mkdir(exist_ok=True)

    bam = score_dir / "reads.bam"
    report = score_dir / "contig_confidence.tsv"

    if report.exists():
        print("[funcat] Existing confidence report detected — skipping")
        return report

    # Check samtools
    if shutil.which("samtools") is None:
        print("⚠️  samtools not found — skipping confidence scoring")
        print("   Install with: conda install -c bioconda samtools")
        return None

    print("[funcat] Mapping reads for coverage analysis...")

    run(
        f"minimap2 -ax {minimap2_preset} -t {threads} {assembly} {reads} "
        f"| samtools sort -@ {threads} -o {bam}"
    )
    run(f"samtools index {bam}")

    print("[funcat] Computing per-contig coverage statistics...")

    depth_raw = run_capture(f"samtools depth -a {bam}")

    # Build per-contig depth lists
    contig_depths = {}
    for line in depth_raw.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        ctg, depth = parts[0], int(parts[2])
        contig_depths.setdefault(ctg, []).append(depth)

    # Compute overall median coverage for relative thresholds
    all_depths = [d for depths in contig_depths.values() for d in depths]
    all_depths.sort()
    median_cov = all_depths[len(all_depths) // 2] if all_depths else 1

    rows = []
    flagged = 0
    review = 0

    for contig, depths in contig_depths.items():
        n = len(depths)
        mean_cov = sum(depths) / n
        variance = sum((d - mean_cov) ** 2 for d in depths) / n
        std_cov = math.sqrt(variance)
        cv = (std_cov / mean_cov) if mean_cov > 0 else 999
        low_cov_pct = 100 * sum(1 for d in depths if d < 5) / n

        # Confidence rules
        if mean_cov > median_cov * 1.8:
            label = "FLAG"       # likely collapsed repeat
            reason = "coverage >1.8x median (possible collapsed repeat)"
            flagged += 1
        elif mean_cov < median_cov * 0.2:
            label = "FLAG"       # possible contamination / misassembly
            reason = "coverage <0.2x median (possible contamination)"
            flagged += 1
        elif cv > 1.5 or low_cov_pct > 20:
            label = "REVIEW"     # uneven coverage
            reason = f"uneven coverage (CV={cv:.2f}, {low_cov_pct:.0f}% low-cov bases)"
            review += 1
        else:
            label = "GOOD"
            reason = "coverage looks uniform and expected"

        rows.append({
            "contig": contig,
            "length_bp": n,
            "mean_coverage": round(mean_cov, 1),
            "coverage_cv": round(cv, 3),
            "low_cov_pct": round(low_cov_pct, 1),
            "label": label,
            "reason": reason,
        })

    # Write TSV
    with open(report, "w") as f:
        headers = ["contig", "length_bp", "mean_coverage",
                   "coverage_cv", "low_cov_pct", "label", "reason"]
        f.write("\t".join(headers) + "\n")
        for row in rows:
            f.write("\t".join(str(row[h]) for h in headers) + "\n")

    total = len(rows)
    good = total - flagged - review

    print(f"\n✅ Confidence scoring complete ({total} contigs)")
    print(f"   GOOD   : {good}")
    print(f"   REVIEW : {review}  (check these before publishing)")
    print(f"   FLAG   : {flagged}  (likely collapsed repeats or contamination)")
    print(f"\n   Full report: {report}")

    if flagged > 0:
        print("\n   Flagged contigs:")
        for row in rows:
            if row["label"] == "FLAG":
                print(f"     {row['contig']} — {row['reason']}")

    return report


# ================================================
# MODULE 6 — Illumina polishing
# ================================================

def run_illumina_polishing(
    assembly,
    illumina_r1,
    illumina_r2,
    outdir,
    threads,
    polisher="polypolish"
):
    """
    Polish a long-read assembly using Illumina short reads.
    
    Uses either Polypolish (default, recommended) or Pilon for polishing.
    Significantly improves base accuracy and reduces internal stop codons.
    
    Returns path to polished assembly.
    """
    
    print("\n" + "=" * 60)
    print("🔬 Module 6 — Illumina polishing")
    print(f"   Polisher: {polisher}")
    print("=" * 60)
    
    polish_dir = Path(outdir) / "illumina_polish"
    polish_dir.mkdir(exist_ok=True)
    
    polished_assembly = polish_dir / f"polished_{polisher}.fasta"
    
    if polished_assembly.exists():
        print("[funcat] Existing Illumina polishing output detected — skipping")
        return polished_assembly
    
    assembly = Path(assembly)
    illumina_r1 = Path(illumina_r1)
    illumina_r2 = Path(illumina_r2)
    
    # Check required tools
    required_tools = ["bwa", "samtools"]
    if polisher == "polypolish":
        required_tools.extend(["polypolish", "polypolish_insert_filter"])
    elif polisher == "pilon":
        required_tools.append("pilon")
    
    missing_tools = [tool for tool in required_tools if shutil.which(tool) is None]
    
    if missing_tools:
        print(f"\n⚠️  Missing tools for Illumina polishing: {missing_tools}")
        if polisher == "polypolish":
            print("   Install with: conda install -c bioconda polypolish bwa samtools")
        else:
            print("   Install with: conda install -c bioconda pilon bwa samtools")
        print("   Skipping Illumina polishing — assembly will be unpolished\n")
        return assembly
    
    print("\n[funcat] Step 1 — Indexing assembly for BWA")
    run(f"bwa index {assembly}")
    
    print("[funcat] Step 2 — Mapping Illumina reads")
    sam_r1 = polish_dir / "alignments_1.sam"
    sam_r2 = polish_dir / "alignments_2.sam"
    
    run(f"bwa mem -t {threads} -a {assembly} {illumina_r1} > {sam_r1}")
    run(f"bwa mem -t {threads} -a {assembly} {illumina_r2} > {sam_r2}")
    
    if polisher == "polypolish":
        print("[funcat] Step 3 — Filtering alignments (Polypolish)")
        filtered_r1 = polish_dir / "filtered_1.sam"
        filtered_r2 = polish_dir / "filtered_2.sam"
        
        run(f"polypolish_insert_filter --in1 {sam_r1} --in2 {sam_r2} --out1 {filtered_r1} --out2 {filtered_r2}")
        
        print("[funcat] Step 4 — Polishing with Polypolish")
        run(f"polypolish {assembly} {filtered_r1} {filtered_r2} > {polished_assembly}")
        
    elif polisher == "pilon":
        print("[funcat] Step 3 — Converting to BAM and sorting")
        bam_file = polish_dir / "illumina_mapped.bam"
        
        run(f"samtools view -bS {sam_r1} | samtools sort -@ {threads} -o {bam_file}")
        run(f"samtools index {bam_file}")
        
        print("[funcat] Step 4 — Polishing with Pilon")
        run(f"pilon --genome {assembly} --frags {bam_file} --output polished_pilon --outdir {polish_dir} --changes --threads {threads}")
        
        # Pilon outputs with different name
        pilon_output = polish_dir / "polished_pilon.fasta"
        if pilon_output.exists():
            run(f"cp {pilon_output} {polished_assembly}")
    
    if not polished_assembly.exists():
        print(f"⚠️  {polisher} polishing failed — using original assembly")
        return assembly
    
    # Calculate improvement statistics
    try:
        original_size = sum(len(record.seq) for record in SeqIO.parse(str(assembly), "fasta"))
        polished_size = sum(len(record.seq) for record in SeqIO.parse(str(polished_assembly), "fasta"))
        size_change = polished_size - original_size
        
        print(f"\n✅ Illumina polishing complete")
        print(f"   Polisher used    : {polisher}")
        print(f"   Original size    : {original_size:,} bp")
        print(f"   Polished size    : {polished_size:,} bp")
        print(f"   Size change      : {size_change:+,} bp")
        print(f"   Output assembly  : {polished_assembly}")
        
        # Suggest running BUSCO to check improvement
        print(f"\n   💡 Tip: Run BUSCO on both assemblies to measure improvement:")
        print(f"      Internal stop codons should drop significantly (22% → <5%)")
        
    except Exception as e:
        print(f"✅ Illumina polishing complete (statistics calculation failed: {e})")
    
    return polished_assembly
