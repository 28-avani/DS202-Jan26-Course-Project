import os
import re
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

# --- CONFIGURATION ---
references = {
    "P_aeruginosa": "sequences/P_aeruginosa.fasta",
    "M_tuberculosis": "sequences/M_tuberculosis.fasta",
    "S_aureus": "sequences/S_aureus.fasta",
    "L_acidophilus": "sequences/L_acidophilus.fasta"
}
runs = ["results/assembly_basic", "results/assembly_mercy", "results/assembly_pruned", "results/assembly_full"]

# 1. Calculate True Total Reference Size for NG50
total_ref_size = 0
for path in references.values():
    if os.path.exists(path):
        total_ref_size += sum(len(r) for r in SeqIO.parse(path, "fasta"))

def get_fasta_stats(filepath, expected_genome_size):
    lengths = []
    gcs = []
    for rec in SeqIO.parse(filepath, "fasta"):
        lengths.append(len(rec))
        gcs.append(gc_fraction(rec.seq) * 100)
    
    if not lengths: return 0, 0, 0, 0, 0, 0, [], []
    
    sorted_lengths = sorted(lengths, reverse=True)
    total_bp = sum(sorted_lengths)
    max_contig = sorted_lengths[0]
    
    # Standard N50 / L50
    csum = np.cumsum(sorted_lengths)
    n50_idx = np.where(csum >= total_bp / 2)[0][0]
    n50 = sorted_lengths[n50_idx]
    l50 = n50_idx + 1
    
    # NG50 / LG50 (Using true genome size)
    try:
        ng50_idx = np.where(csum >= expected_genome_size / 2)[0][0]
        ng50 = sorted_lengths[ng50_idx]
        lg50 = ng50_idx + 1
    except IndexError:
        ng50, lg50 = 0, 0 # Assembly is less than half the expected size!

    # auN (Area under the Nx curve)
    aun = sum(l**2 for l in sorted_lengths) / total_bp if total_bp > 0 else 0
    
    return total_bp, max_contig, n50, l50, ng50, lg50, aun, sorted_lengths, gcs

def parse_megahit_log(log_path):
    """Extracts RAM (bytes) and Time (seconds) from the log."""
    time_s, max_ram_b = 0.0, 0
    if not os.path.exists(log_path): return time_s, max_ram_b
    
    with open(log_path, 'r') as f:
        content = f.read()
        # Regex to catch: "Real: 12.34s, max rss: 123456"
        time_match = re.search(r'Real: ([\d\.]+)s', content)
        ram_match = re.search(r'max rss: (\d+)', content)
        
        if time_match: time_s = float(time_match.group(1))
        if ram_match: max_ram_b = int(ram_match.group(1))
            
    return time_s, max_ram_b

def calculate_recovery_and_chimeras(ref_dict, query_fasta):
    """Maps to all refs and calculates recovery + rough chimera count."""
    if not os.path.exists(query_fasta):
        return {k: 0.0 for k in ref_dict}, 0
        
    recoveries = {}
    aligned_contigs = {} # Map contig_id to a set of species it aligned to
    
    for species, ref_path in ref_dict.items():
        if not os.path.exists(ref_path): continue
        cmd = f"minimap2 -x asm5 {ref_path} {query_fasta} 2>/dev/null"
        try:
            output = subprocess.check_output(cmd, shell=True).decode().strip()
            aligned_bases = 0
            for line in output.split('\n'):
                if not line: continue
                parts = line.split('\t')
                q_name, mapped_len = parts[0], int(parts[10])
                aligned_bases += mapped_len
                
                if mapped_len > 500: # Only count significant alignments
                    if q_name not in aligned_contigs: aligned_contigs[q_name] = set()
                    aligned_contigs[q_name].add(species)
                    
            total_ref_len = sum(len(r) for r in SeqIO.parse(ref_path, "fasta"))
            recoveries[species] = min(100.0, (aligned_bases / total_ref_len) * 100)
        except:
            recoveries[species] = 0.0

    # Chimeras: Contigs that map significantly to >1 species
    chimeras = sum(1 for hits in aligned_contigs.values() if len(hits) > 1)
    return recoveries, chimeras

data = []
all_lens = {}

runs_paths = [f"results/{r}" if not os.path.exists(r) else r for r in runs]

for run, path in zip(runs, runs_paths):
    fasta = os.path.join(path, "final.contigs.fa")
    log = os.path.join(path, "megahit.log")
    if not os.path.exists(fasta): continue

    tot_bp, max_c, n50, l50, ng50, lg50, aun, lens, gcs = get_fasta_stats(fasta, total_ref_size)
    time_s, ram_b = parse_megahit_log(log)
    recs, chimeras = calculate_recovery_and_chimeras(references, fasta)
    
    all_lens[run] = lens
    
    row = {
        "Phase": run.split("_")[1].capitalize(),
        "Total_MB": round(tot_bp / 1000000, 2),
        "NG50": ng50,
        "auN": int(aun),
        "Chimeras": chimeras,
        "RAM_MB": round(ram_b / 1048576, 1) if ram_b else 0,
        "Time_s": time_s
    }
    for species, rec in recs.items():
        row[f"{species[:3]}_Rec_%"] = round(rec, 2)
        
    data.append(row)

df = pd.DataFrame(data)

print("\n" + "="*110)
print("ALGORITHMIC BIOLOGY PROJECT: ADVANCED PIPELINE METRICS")
print("="*110)
print(df.to_string(index=False))
print("="*110)

# --- PLOTTING (Advanced 6-Panel) ---
fig, axes = plt.subplots(2, 3, figsize=(20, 10))

# 1. Recovery
colors = ["#2ca02c", "#9467bd", "#1f77b4", "#d62728"]
rec_cols = ["P_a_Rec_%", "M_t_Rec_%", "S_a_Rec_%", "L_a_Rec_%"]
df.plot(x="Phase", y=rec_cols, kind="bar", ax=axes[0,0], color=colors)
axes[0,0].set_title("1. Sensitivity (Mercy Rescue)")
axes[0,0].set_ylabel("Genome Recovery (%)")
axes[0,0].legend(["P. aeruginosa", "M. tuberculosis", "S. aureus", "L. acidophilus"], fontsize='small')

# 2. NG50 vs N50 (True Contiguity)
df.plot(x="Phase", y=["NG50", "auN"], kind="line", marker='o', ax=axes[0,1])
axes[0,1].set_title("2. True Contiguity (NG50 & auN)")
axes[0,1].set_yscale('log')
axes[0,1].set_ylabel("Base Pairs (Log Scale)")


# 3. Cumulative Assembly Curve (NGx)
for run_name, lengths in all_lens.items():
    if not lengths: continue
    phase_name = run_name.split("_")[-1].capitalize()
    y = np.cumsum(lengths) / 1000000
    x = np.arange(1, len(lengths) + 1)
    axes[0,2].plot(x, y, label=phase_name)
    
axes[0,2].set_title("3. NGx Cumulative Assembly Curve")
axes[0,2].set_xlabel("Number of Contigs")
axes[0,2].set_ylabel("Cumulative Size (MB)")
axes[0,2].set_xlim(0, max(len(l) for l in all_lens.values()) * 0.5) # Zoom in
axes[0,2].legend()

plt.tight_layout()
plt.savefig("advanced_ablation_metrics.png", dpi=300)
print("\nHigh-res plot saved to 'advanced_ablation_metrics.png'")