import os
import math
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from collections import Counter

# --- CONFIGURATION ---
points = ["pt1", "pt2", "pt3", "pt4", "pt5", "pt6", "pt7", "pt8", "pt9", "pt10"]
BASE_DIR = "results_entropy"

def calc_shannon_entropy(seq_str):
    if not seq_str: return 0
    counts = Counter(seq_str)
    length = len(seq_str)
    return -sum((count/length) * math.log2(count/length) for count in counts.values())

def get_ram_from_megahit_log(log_path):
    """Extracts the ABSOLUTE MAXIMUM maxrss from the log."""
    if not os.path.exists(log_path): 
        return 0
    with open(log_path, 'r') as f:
        content = f.read()
        # Find ALL maxrss mentions and take the highest one
        matches = re.findall(r'maxrss:\s+(\d+)', content)
        if matches:
            max_kb = max([int(m) for m in matches])
            return max_kb / 1024 # Convert KB to MB
    return 0

def get_stats(run_id):
    # Matches your exact VS Code file structure
    folder_name = f"ref_{run_id}"
    fasta_path = os.path.join(BASE_DIR, folder_name, "final.contigs.fa")
    log_path = os.path.join(BASE_DIR, folder_name, "log")
    
    if not os.path.exists(fasta_path):
        print(f"⚠️ Skipping {run_id}: Could not find {fasta_path}")
        return None

    lengths = [len(rec) for rec in SeqIO.parse(fasta_path, "fasta")]
    if not lengths: 
        print(f"⚠️ Skipping {run_id}: Fasta file is empty.")
        return None
    
    full_seq = "".join(str(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta"))
    entropy = calc_shannon_entropy(full_seq)
    
    total_bp = sum(lengths)
    lengths.sort(reverse=True)
    n50 = lengths[np.where(np.cumsum(lengths) >= total_bp/2)[0][0]]
    ram = get_ram_from_megahit_log(log_path)
    
    efficiency = ram / (total_bp / 1e6) if total_bp > 0 else 0 

    return {
        "Point": run_id,
        "Entropy": round(entropy, 3),
        "RAM_MB": round(ram, 2),
        "N50": n50,
        "Total_MB": round(total_bp / 1e6, 2),
        "MB_RAM_per_Mbp": round(efficiency, 2)
    }

print("Analyzing MEGAHIT assemblies...")
results = []
for p in points:
    stat = get_stats(p)
    if stat:
        results.append(stat)

if not results:
    print("\nCRITICAL ERROR: No results were collected. Check your folder structure.")
else:
    df = pd.DataFrame(results).sort_values("Entropy")
    
    # --- PLOTTING ---
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    
    # Plot 1: SdBG Memory vs. Entropy
    axes[0].scatter(df["Entropy"], df["RAM_MB"], color='red', s=100, edgecolors='black')
    axes[0].plot(df["Entropy"], df["RAM_MB"], color='red', linestyle='--', alpha=0.6)
    axes[0].set_title("1. SdBG Memory vs. Information Entropy")
    axes[0].set_xlabel("Shannon Entropy (bits/base)")
    axes[0].set_ylabel("Peak RAM (MB)")
    axes[0].grid(True, alpha=0.3)
    
    # Plot 2: N50 (Connectivity)
    axes[1].plot(df["Entropy"], df["N50"], marker='o', color='teal', linewidth=2)
    axes[1].set_title("2. Graph Connectivity (N50)")
    axes[1].set_xlabel("Shannon Entropy (bits/base)")
    axes[1].set_ylabel("N50 (bp)")
    axes[1].grid(True, alpha=0.3)
    
    # Plot 3: Effective Genome Coverage
    axes[2].bar(df["Point"], df["Total_MB"] / 2.0, color='skyblue', alpha=0.8)
    axes[2].set_title("3. Effective Genome Coverage vs. Entropy")
    axes[2].set_xlabel("Entropy Point")
    axes[2].set_ylabel("Assembled Fold-Coverage (x)")
    axes[2].set_xticklabels(df["Point"], rotation=45)

    plt.tight_layout()
    plt.savefig("entropy_10point_analysis.png", dpi=300)
    
    print("\n" + "="*80)
    print("EXPERIMENT 6: 10-POINT ENTROPY ANALYSIS")
    print("="*80)
    print(df.to_string(index=False))
    print("="*80)
    print("\n✅ Plot saved to entropy_10point_analysis.png")