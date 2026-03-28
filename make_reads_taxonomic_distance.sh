#!/bin/bash
echo "Simulating reads for all genomes..."
wgsim -N 300000 -1 100 -2 100 S_aureus_8325.fasta sa_anchor.fq /dev/null
wgsim -N 300000 -1 100 -2 100 L_monocytogenes.fasta l_mono.fq /dev/null
wgsim -N 300000 -1 100 -2 100 S_epidermidis.fasta s_epi.fq /dev/null
wgsim -N 300000 -1 100 -2 100 S_aureus_USA300.fasta sa_usa300.fq /dev/null

echo "Mixing Pair 1: Low Similarity (S. aureus + L. mono)"
cat sa_anchor.fq l_mono.fq > mix_low.fq

echo "Mixing Pair 2: Medium Similarity (S. aureus + S. epidermidis)"
cat sa_anchor.fq s_epi.fq > mix_med.fq

echo "Mixing Pair 3: High Similarity (S. aureus 8325 + S. aureus USA300)"
cat sa_anchor.fq sa_usa300.fq > mix_high.fq

echo "Done! 3 mixed FASTQ files created."