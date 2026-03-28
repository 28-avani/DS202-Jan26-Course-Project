import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO

# --- CONFIGURATION ---
experiments = {
    "Low (<70% ANI)": {
        "folder": "results_taxonomic_distance/run_low_sim",
        "ref1": "S_aureus_8325.fasta",
        "ref2": "L_monocytogenes.fasta"
    },
    "Med (~80% ANI)": {
        "folder": "results_taxonomic_distance/run_med_sim",
        "ref1": "S_aureus_8325.fasta",
        "ref2": "S_epidermidis.fasta"
    },
    "High (>99% ANI)": {
        "folder": "results_taxonomic_distance/run_high_sim",
        "ref1": "S_aureus_8325.fasta",
        "ref2": "S_aureus_USA300.fasta"
    }
}

def get_mapped_contigs(ref_fasta, query_fasta):
    mapped = {}
    if not os.path.exists(ref_fasta):
        print(f"  ❌ Missing reference FASTA: {ref_fasta}")
        return mapped
    if not os.path.exists(query_fasta):
        print(f"  ❌ Missing assembly file: {query_fasta}")
        return mapped
    
    cmd = f"minimap2 -x asm5 {ref_fasta} {query_fasta} 2>/dev/null"
    try:
        output = subprocess.check_output(cmd, shell=True).decode().strip()
        for line in output.split('\n'):
            if not line: continue
            parts = line.split('\t')
            # Ensure line has enough columns
            if len(parts) > 10:
                contig_name = parts[0]
                mapped_len = int(parts[10]) # Number of matching bases
                
                # We only care about substantial alignments to avoid false noise
                if mapped_len > 500: 
                    if contig_name not in mapped:
                        mapped[contig_name] = 0
                    mapped[contig_name] += mapped_len
    except Exception as e:
        print(f"  ❌ Mapping failed for {ref_fasta}: {e}")
        
    return mapped

def get_n50(fasta_path):
    if not os.path.exists(fasta_path): return 0
    lengths = sorted([len(rec) for rec in SeqIO.parse(fasta_path, "fasta")], reverse=True)
    if not lengths: return 0
    total_bp = sum(lengths)
    return lengths[np.where(np.cumsum(lengths) >= total_bp/2)[0][0]]

results = []
print("Running Taxonomic Resolution Analysis (This may take a minute...)...\n")

for label, data in experiments.items():
    print(f"Checking {label}...")
    assembly_path = os.path.join(data["folder"], "final.contigs.fa")
    
    if not os.path.exists(assembly_path):
        print(f"  ❌ Assembly not found: {assembly_path}")
        results.append({"Similarity": label, "N50 (bp)": 0, "Chimeras": 0})
        continue
        
    n50 = get_n50(assembly_path)
    map1 = get_mapped_contigs(data["ref1"], assembly_path)
    map2 = get_mapped_contigs(data["ref2"], assembly_path)
    
    chimeras = sum(1 for contig in map1 if contig in map2)
            
    results.append({
        "Similarity": label,
        "N50 (bp)": n50,
        "Chimeras": chimeras
    })
    print(f"  ✅ Analyzed successfully.")

df = pd.DataFrame(results)

print("\n" + "="*60)
print("EXPERIMENT 5: TAXONOMIC RESOLUTION RESULTS")
print("="*60)
print(df.to_string(index=False))
print("="*60)

# --- PLOTTING ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].bar(df["Similarity"], df["Chimeras"], color='darkred', alpha=0.8)
axes[0].set_title("1. Inter-Species Chimeras vs. Genome Similarity")
axes[0].set_xlabel("Genomic Similarity (ANI)")
axes[0].set_ylabel("Number of Chimeric Contigs")

axes[1].plot(df["Similarity"], df["N50 (bp)"], marker='o', color='teal', linewidth=2, markersize=8)
axes[1].set_title("2. Assembly Fragmentation (N50) vs. Similarity")
axes[1].set_xlabel("Genomic Similarity (ANI)")
axes[1].set_ylabel("N50 (bp)")

plt.tight_layout()
plt.savefig("taxonomic_resolution_exp5.png", dpi=300)
print("\n✅ Plot saved to taxonomic_resolution_exp5.png")