"""
make_random_reads.py
────────────────────
Reads random_genome.fasta (produced by build_random_genome.py) and generates
a FASTQ file of simulated reads using the same strategy as make_reads_ablation.c:
  • pick a random start position in the genome
  • slice out a fixed-length substring
  • write a FASTQ record with a perfect Phred+33 quality string ('I')

Output: mock_random_exp1.fq

Usage:
    python make_random_reads.py [--num-reads N] [--read-length L] [--seed S]
"""

import argparse
import random
import sys

# ── defaults (mirror make_reads_ablation.c) ───────────────────────────────────

DEFAULT_NUM_READS   = 20_000
DEFAULT_READ_LENGTH = 300
DEFAULT_SEED        = None       # None → non-deterministic (time-based)
INPUT_FASTA         = "random_genome.fasta"
OUTPUT_FASTQ        = "mock_random_exp1.fq"

# ── CLI ───────────────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser(description="Generate simulated reads from random_genome.fasta")
parser.add_argument("--num-reads",   type=int, default=DEFAULT_NUM_READS,
                    help=f"Number of reads to generate (default: {DEFAULT_NUM_READS})")
parser.add_argument("--read-length", type=int, default=DEFAULT_READ_LENGTH,
                    help=f"Read length in bp (default: {DEFAULT_READ_LENGTH})")
parser.add_argument("--seed",        type=int, default=DEFAULT_SEED,
                    help="Random seed (default: non-deterministic)")
args = parser.parse_args()

NUM_READS   = args.num_reads
READ_LENGTH = args.read_length
SEED        = args.seed

# ── load genome ───────────────────────────────────────────────────────────────

print(f"Loading genome from {INPUT_FASTA} …")
seq_parts = []
with open(INPUT_FASTA) as fh:
    for line in fh:
        if line.startswith(">"):
            continue
        seq_parts.append(line.strip().upper())
genome = "".join(seq_parts)
print(f"  Genome length: {len(genome):,} bp")

if READ_LENGTH >= len(genome):
    sys.exit(f"Error: read length ({READ_LENGTH}) must be shorter than genome ({len(genome):,} bp)")

# ── generate reads ────────────────────────────────────────────────────────────

if SEED is not None:
    random.seed(SEED)

qual = "I" * READ_LENGTH     # perfect Phred+33 quality (same as the C code)
max_start = len(genome) - READ_LENGTH

print(f"Generating {NUM_READS:,} reads of length {READ_LENGTH} bp …")
with open(OUTPUT_FASTQ, "w") as out:
    for i in range(NUM_READS):
        start = random.randint(0, max_start)
        read  = genome[start:start + READ_LENGTH]
        out.write(f"@Read_{i}_Source_random\n{read}\n+\n{qual}\n")

print(f"Written to {OUTPUT_FASTQ}  ({NUM_READS:,} reads × {READ_LENGTH} bp)")
