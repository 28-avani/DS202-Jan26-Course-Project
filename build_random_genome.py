"""
build_random_genome.py
──────────────────────
Constructs a random ACGT genome of the same length as the L. acidophilus
genome (2,022,041 bp) and writes it to a FASTA file.

Output: random_genome.fasta
"""

import numpy as np

# ── parameters ────────────────────────────────────────────────────────────────

L_ACIDOPHILUS_GENOME_LEN = 1_993_560
OUTPUT_FASTA             = "random_genome.fasta"
LINE_WIDTH               = 60          # bases per FASTA line
SEED                     = 42          # set to None for a different genome each run

# ── build ─────────────────────────────────────────────────────────────────────

print(f"Building random genome ({L_ACIDOPHILUS_GENOME_LEN:,} bp, seed={SEED}) …")

rng   = np.random.default_rng(SEED)
idx   = rng.integers(0, 4, size=L_ACIDOPHILUS_GENOME_LEN, dtype=np.uint8)
bases = np.array(['A', 'C', 'G', 'T'], dtype='U1')
seq   = "".join(bases[idx])

# ── write FASTA ───────────────────────────────────────────────────────────────

with open(OUTPUT_FASTA, "w") as fh:
    fh.write(">Random_genome_L_acidophilus_size\n")
    for i in range(0, len(seq), LINE_WIDTH):
        fh.write(seq[i:i + LINE_WIDTH] + "\n")

print(f"Written to {OUTPUT_FASTA}  ({len(seq):,} bp)")
