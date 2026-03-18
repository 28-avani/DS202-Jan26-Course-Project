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

# References
(1) A. Bowe, T. Onodera, K. Sadakane, and T. Shibuya. 2012. Succinct de Bruijn Graphs. Raphael, B., Tang, J. (eds) Algorithms in Bioinformatics. WABI 2012. Lecture Notes in Computer Science(), vol 7534. Springer, Berlin, Heidelberg (2012). doi:10.1007/978-3-642-33122-0_18

(2) D. Li, Liu, C. M., R. Luo, K. Sadakane, and T. W. Lam. 2015. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics, 31(10), 1674-1676 (2015). doi:10.1093/bioinformatics/btv033

(3) Yu Peng, Henry C. M. Leung, S. M. Yiu, and Francis Y. L. Chin. 2012. IDBA-UD: a de novo assembler for single-cell and metagenomic sequencing data with highly uneven depth. Bioinformatics, Volume 28, Issue 11, June 2012, Pages 1420–1428 (2012). doi:10.1093/bioinformatics/bts174