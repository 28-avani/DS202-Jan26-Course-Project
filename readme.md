# 1 Problem Description

Assembling metagenomic data from scratch takes a huge amount of memory if we use standard de Bruijn graphs. Real-world samples contain a mixture of microbes with varying abundances, requiring a lot of RAM to process via traditional methods. Essentially, the algorithmic problem this project tackles is: how can we represent and navigate this massive graph to generate contigs without running out of memory, while still being able to distinguish real biological variations and sequencing errors in rare species.

# 2 Importance and Applications of this Problem

To study microbial communities like those found in the human/animal gut or in environmental/soil samples, we don’t have existing reference genomes to map our reads against, because most of these microbes are unknown. De novo metagenomic assembly is therefore the first step to ID the microbes present in the sample. It is computationally expensive to use standard procedures to sequence these genomes. Thus, this memory-optimization problem is essential for practical research.
Solving this enables many applications like:
- Large-Scale Environmental Microbiome Profiling
- Clinical Pathogen Discovery in Unculturable Samples
- AMR Gene Detection and Tracking
- Human Gut Microbiome Characterization
- Novel Biocatalyst and Enzyme Discovery
- Wastewater Bio-surveillance and Outbreak Tracking
- Soil Health and Plant-Microbe Interaction Analysis

# 3 Topics and Sub-topics for the Lecture

## 3.1 Introduction
First, we explain why de Bruijn graphs fail on mixed-genome samples, and highlight how in biological samples, highly abundant species drown out the rare species.

## 3.2 Constructing the Succinct de Bruijn Graph
We then explain how the graph is compressed into $O(m)$ bits and the $(k+1)$-mer edges are sorted in reverse lexicographical order to construct a Burrows-Wheeler Transform string, in contrast to the regular de Bruijn Graph. We also demonstrate how to use Last-to-First (LF) mapping and rank/select queries on the BWT string to navigate from one DNA node to its neighbours without explicitly storing the graph in memory. (1)

## 3.3 Addressing Errors and Rare Species
We explain how to discard the low-frequency $k$-mers to remove sequencing errors, and the “mercy $k$-mer” rescue strategy to recover actual low-depth edges from rare species that would otherwise be lost. (3)

## 3.4 Dynamic Graph Pruning
We explain the algorithm to detect and remove dead-end paths caused by sequencing errors at the ends of reads (tip removal) and how the algorithm uses breadth-first search to detect diverging-then-converging paths caused by strain variations or middle-read errors and collapses them by selecting the path with higher sequence coverage (bubble merging) and also progressively flips the valid bit-vector to 0 for edges that have suspiciously low local depth compared to their neighbours (low coverage pruning).

## 3.5 The Iterative Multiple $k$-mer Strategy
Finally, we explain why a single $k$-mer size fails for metagenomics. Small $k$-mers bridge gaps in low-coverage regions, while large $k$-mers are required to resolve repetitive DNA. (3). We then summarise the MEGAHIT algorithm loop: building an SdBG at $k_{min}$, cleaning it, extracting contigs, and then stepping up the $k$ size and showing how the algorithm extracts new $(k+step+1)$-mers from both the original reads and the newly formed contigs to build the next higher-order SdBG, repeating until $k_{max}$ is reached to give the final assembly. (2)

# 4 Planned Experiments
To evaluate the algorithmic design, the experiments systematically enable distinct components of the algorithm on simulated datasets, to quantify the exact problem each component solves.
In addition to the component-wise analysis detailed in Table 1 Rows 1-4, Experiment 5 tests the edge-case scenario to find points where the algorithm fails, by evaluating sequences across a gradient of taxonomic relatedness, we will demonstrate how the graph remains permanently tangled (failing to merge bubbles) if the shared genomic sequence between two strains strictly exceeds the length of $k_{max}$. Similarly, Phase 6 will establish an information-theoretic baseline. By comparing the assembly of a true biological genome against a purely random nucleotide sequence of identical length, we expect to demonstrate that the SdBG’s compression efficiency and graph connectivity fundamentally rely on the inherent redundancy of biological data. The random sequence is hypothesized to yield a poor BWT compression ratio and fail to form contiguous paths.

![](/images/experiments.png)

# 5 Results

## Experiment 1: Memory Scaling

This experiment measures and compares the memory usage of a Standard de Bruijn Graph (DBG) and a Succinct DBG (BOSS representation) as a function of read count and read length, for both random and biological input data.

### Prerequisites

Install dependencies into the virtual environment:

```bash
pip install -r requirements.txt
```

### Step 1 — Generate the random genome

```bash
python build_random_genome.py
```

Constructs a 1,993,560 bp random ACGT genome (matching the *L. acidophilus* genome length) and writes it to `random_genome.fasta`.

### Step 2 — Generate simulated reads from the random genome

```bash
python make_random_reads.py
```

Samples 20,000 reads of 1,000 bp from `random_genome.fasta` using the same random-substring strategy as `make_reads_ablation.c`. Writes output to `mock_random_exp1.fq`.

Optional flags:
```bash
python make_random_reads.py --num-reads 20000 --read-length 1000 --seed 42
```

### Step 3 — Generate biological reads

Compile and run `make_reads_ablation_exp1.c` to produce `mock_metagenome_exp1.fq` from the *L. acidophilus* genome:

```bash
gcc -O2 -o make_reads_ablation make_reads_ablation_exp1.c -lm
./make_reads_ablation
```

### Step 4 — Run the memory scaling benchmark and plot

```bash
python memory_scaling_plots.py
```

Produces a two-panel figure:
- **Left:** Memory usage vs. number of reads (2 – 1,000), at a fixed read length of 150 bp
- **Right:** Memory usage vs. read length (20 – 1,000 bp), at a fixed read count of 1,000

Each panel overlays four series: Standard DBG and Succinct DBG for both random and biological data.

# References
(1) A. Bowe, T. Onodera, K. Sadakane, and T. Shibuya. 2012. Succinct de Bruijn Graphs. Raphael, B., Tang, J. (eds) Algorithms in Bioinformatics. WABI 2012. Lecture Notes in Computer Science(), vol 7534. Springer, Berlin, Heidelberg (2012). doi:10.1007/978-3-642-33122-0_18

(2) D. Li, Liu, C. M., R. Luo, K. Sadakane, and T. W. Lam. 2015. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics, 31(10), 1674-1676 (2015). doi:10.1093/bioinformatics/btv033

(3) Yu Peng, Henry C. M. Leung, S. M. Yiu, and Francis Y. L. Chin. 2012. IDBA-UD: a de novo assembler for single-cell and metagenomic sequencing data with highly uneven depth. Bioinformatics, Volume 28, Issue 11, June 2012, Pages 1420–1428 (2012). doi:10.1093/bioinformatics/bts174