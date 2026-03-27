import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

# --- CONFIGURATION ---
references = {
    "P_aeruginosa": "P_aeruginosa.fasta",
    "M_tuberculosis": "M_tuberculosis.fasta",
    "S_aureus": "S_aureus.fasta",
    "L_acidophilus": "L_acidophilus.fasta"
}
runs = ["results/assembly_basic", "results/assembly_mercy", "results/assembly_pruned", "results/assembly_full"]

def get_fasta_stats(filepath):
    lengths = []
    gcs = []
    for rec in SeqIO.parse(filepath, "fasta"):
        lengths.append(len(rec))
        gcs.append(gc_fraction(rec.seq) * 100)
    
    if not lengths: return 0, 0, [], []
    return sum(lengths), max(lengths), sorted(lengths, reverse=True), gcs

def get_n50_l50(lengths):
    if not lengths: return 0, 0
    total_len = sum(lengths)
    csum = np.cumsum(lengths)
    n50_idx = np.where(csum >= total_len / 2)[0][0]
    return lengths[n50_idx], n50_idx + 1

def calculate_recovery(ref_fasta, query_fasta):
    if not os.path.exists(ref_fasta) or not os.path.exists(query_fasta):
        return 0.0
    
    cmd = f"minimap2 -x asm5 {ref_fasta} {query_fasta} 2>/dev/null | awk '{{sum += $10}} END {{print sum}}'"
    try:
        output = subprocess.check_output(cmd, shell=True).decode().strip()
        aligned_bases = int(output) if output else 0
    except:
        aligned_bases = 0
        
    total_ref_len = sum(len(r) for r in SeqIO.parse(ref_fasta, "fasta"))
    return min(100.0, (aligned_bases / total_ref_len) * 100)

data = []
all_gcs = {}

# Make sure we look in the 'results' folder if that's where your runs are
# If your folders are just 'assembly_basic', remove 'results/' from the runs list.
runs_paths = [f"results/{r}" if not os.path.exists(r) else r for r in runs]

for run, path in zip(runs, runs_paths):
    fasta = os.path.join(path, "final.contigs.fa")
    if not os.path.exists(fasta): 
        continue

    total_bp, max_contig, lengths, gcs = get_fasta_stats(fasta)
    n50, l50 = get_n50_l50(lengths)
    all_gcs[run] = gcs
    
    row = {
        "Phase": run.split("_")[1].capitalize(),
        "Total_BP": total_bp,
        "Max_Contig": max_contig,
        "N50": n50,
        "L50": l50
    }
    
    for species, ref_path in references.items():
        row[f"{species[:3]}_Rec_%"] = calculate_recovery(ref_path, fasta)
        
    data.append(row)

df = pd.DataFrame(data)

print("\n" + "="*100)
print("ALGORITHMIC BIOLOGY PROJECT: 4-SPECIES ABLATION STUDY RESULTS")
print("="*100)
print(df.to_string(index=False))
print("="*100)

# --- PLOTTING ---
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Plot 1: Recovery for ALL FOUR Species
recovery_cols = ["P_a_Rec_%", "M_t_Rec_%", "S_a_Rec_%", "L_a_Rec_%"]
colors = ["#2ca02c", "#9467bd", "#1f77b4", "#d62728"] # Green, Purple, Blue, Red

df.plot(x="Phase", y=recovery_cols, kind="bar", ax=axes[0,0], color=colors)
axes[0,0].set_title("Exp 2 & 4: Genome Recovery Across Pipeline")
axes[0,0].set_ylabel("Genome Recovery (%)")
axes[0,0].legend(["P. aeruginosa", "M. tuberculosis", "S. aureus", "L. acidophilus"], title="Species")

# Plot 2: Fragmentation (L50)
df.plot(x="Phase", y="L50", kind="line", marker='v', color='darkorange', ax=axes[0,1])
axes[0,1].set_title("Exp 3: Assembly Fragmentation (L50 - Lower is Better)")
axes[0,1].set_ylabel("L50 (Number of Contigs to reach 50%)")
axes[0,1].invert_yaxis()

# Plot 3: Spanning (Max Contig & N50)
df.plot(x="Phase", y=["Max_Contig", "N50"], kind="line", marker='o', ax=axes[1,0])
axes[1,0].set_title("Exp 4: Repeat Spanning (Iterative k-mer)")
axes[1,0].set_ylabel("Base Pairs (bp)")

# Plot 4: GC Content Histogram (Using the 'Full' assembly)
if "assembly_full" in all_gcs or "results/assembly_full" in all_gcs:
    key = "assembly_full" if "assembly_full" in all_gcs else "results/assembly_full"
    axes[1,1].hist(all_gcs[key], bins=50, color="purple", alpha=0.7)
    axes[1,1].set_title("Taxonomic Resolution (GC Content in Full Assembly)")
    axes[1,1].set_xlabel("GC Percentage (%)")
    axes[1,1].set_ylabel("Number of Contigs")
    axes[1,1].axvline(x=34, color='blue', linestyle='--', label='Low GC Peak (~34%)')
    axes[1,1].axvline(x=65, color='red', linestyle='--', label='High GC Peak (~65%)')
    axes[1,1].legend()

plt.tight_layout()
plt.savefig("ablation_results_4species.png")
print("\nGraphs saved as 'ablation_results_4species.png'")