import random
import numpy as np
import matplotlib.pyplot as plt
from memory_scaling_compressed import StandardDBG, SuccinctDBG

# ── helpers ──────────────────────────────────────────────────────────────────

# L. acidophilus genome length — used so the random genome is the same size,
# making the two read-sampling distributions directly comparable.
L_ACIDOPHILUS_GENOME_LEN = 1_993_560

def build_random_genome(length):
    """Generate a random ACGT sequence of *length* bp as a plain string."""
    rng = np.random.default_rng()
    idx = rng.integers(0, 4, size=length, dtype=np.uint8)
    bases = np.array(['A', 'C', 'G', 'T'], dtype='U1')
    return "".join(bases[idx])


def sample_reads_from_genome(genome, count, length):
    """
    Mimic make_reads_ablation.c: pick *count* random start positions and
    return fixed-length substrings — the same strategy used to produce the
    biological FASTQ.
    """
    genome_len = len(genome)
    if length >= genome_len:
        raise ValueError(f"Read length {length} >= genome length {genome_len}")
    starts = random.sample(range(genome_len - length), count)
    return [genome[s:s + length] for s in starts]


def load_bio_reads(fq_path):
    """Return all sequence strings from a FASTQ file."""
    reads = []
    with open(fq_path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().strip()
            fh.readline()   # '+'
            fh.readline()   # quality
            reads.append(seq)
    return reads


def measure(reads):
    sd = StandardDBG(K)
    sd.build_from_reads(reads)
    ss = SuccinctDBG(sd)
    ss.convert_to_boss()
    std = sd.nodes.nbytes + sd.edge_labels.nbytes + sd.edge_offsets.nbytes
    suc = ss.W.nbytes + ss.last.nbytes
    return std, suc

# ── benchmarks ────────────────────────────────────────────────────────────────

def bench_vs_count_random(counts, length, genome):
    """Sample *count* reads of fixed *length* from the random genome."""
    std_mem, suc_mem = [], []
    for n in counts:
        reads = sample_reads_from_genome(genome, n, length)
        std, suc = measure(reads)
        std_mem.append(std)
        suc_mem.append(suc)
    return std_mem, suc_mem


def bench_vs_count_bio(counts, all_bio_reads):
    """Sample *count* reads from the biological pool."""
    std_mem, suc_mem = [], []
    for n in counts:
        reads = random.sample(all_bio_reads, n)
        std, suc = measure(reads)
        std_mem.append(std)
        suc_mem.append(suc)
    return std_mem, suc_mem


def bench_vs_length_random(lengths, count, genome):
    """Sample *count* reads at each target length from the random genome."""
    std_mem, suc_mem = [], []
    for length in lengths:
        reads = sample_reads_from_genome(genome, count, length)
        std, suc = measure(reads)
        std_mem.append(std)
        suc_mem.append(suc)
    return std_mem, suc_mem


def bench_vs_length_bio(lengths, all_bio_reads, count):
    """Trim *count* biological reads to each target length."""
    std_mem, suc_mem = [], []
    pool = random.sample(all_bio_reads, count)
    for length in lengths:
        reads = [r[:length] for r in pool]
        std, suc = measure(reads)
        std_mem.append(std)
        suc_mem.append(suc)
    return std_mem, suc_mem

# ── parameters ────────────────────────────────────────────────────────────────

K           = 31
FASTQ       = "mock_metagenome_exp1.fq"
READ_LENGTH = 150   # fixed length for the vs-count benchmark (matches bio reads)
READ_COUNT  = 1000   # fixed count  for the vs-length benchmark

counts      = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1000]
rnd_lengths = [20, 40, 60, 80, 100, 150, 200, 300, 500, 750, 1000]  # random genome is large enough for 1000 bp
bio_lengths = [20, 40, 60, 80, 100, 150, 200, 300, 500, 750, 1000]  # trimmed from bio reads

# ── build random genome (same size as L. acidophilus) ─────────────────────────

print(f"Building random genome ({L_ACIDOPHILUS_GENOME_LEN:,} bp) …")
rnd_genome = build_random_genome(L_ACIDOPHILUS_GENOME_LEN)

# ── load biological data ──────────────────────────────────────────────────────

all_bio = load_bio_reads(FASTQ)
print(f"Loaded {len(all_bio):,} biological reads (length {len(all_bio[0])} bp)")

# ── run benchmarks ────────────────────────────────────────────────────────────

print("Benchmarking memory vs. read count …")
std_cnt_rnd, suc_cnt_rnd = bench_vs_count_random(counts, READ_LENGTH, rnd_genome)
std_cnt_bio, suc_cnt_bio = bench_vs_count_bio(counts, all_bio)

print("Benchmarking memory vs. read length …")
std_len_rnd, suc_len_rnd = bench_vs_length_random(rnd_lengths, READ_COUNT, rnd_genome)
std_len_bio, suc_len_bio = bench_vs_length_bio(bio_lengths, all_bio, READ_COUNT)

# ── plot ──────────────────────────────────────────────────────────────────────

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# --- left: memory vs. read count ---
ax1.plot(counts, std_cnt_rnd, label='Standard DBG (random)',     marker='o', color='steelblue')
ax1.plot(counts, suc_cnt_rnd, label='Succinct DBG (random)',     marker='s', color='steelblue',  linestyle='--')
ax1.plot(counts, std_cnt_bio, label='Standard DBG (biological)', marker='o', color='darkorange')
ax1.plot(counts, suc_cnt_bio, label='Succinct DBG (biological)', marker='s', color='darkorange', linestyle='--')
ax1.set_xlabel(f'Number of Reads (length={READ_LENGTH} bp)')
ax1.set_ylabel('Memory Usage (Bytes)')
ax1.set_title('Memory vs. Read Count')
ax1.legend()
ax1.grid(True, linestyle='--', alpha=0.6)

# --- right: memory vs. read length ---
ax2.plot(rnd_lengths, std_len_rnd, label='Standard DBG (random)',     marker='o', color='steelblue')
ax2.plot(rnd_lengths, suc_len_rnd, label='Succinct DBG (random)',     marker='s', color='steelblue',  linestyle='--')
ax2.plot(bio_lengths, std_len_bio, label='Standard DBG (biological)', marker='o', color='darkorange')
ax2.plot(bio_lengths, suc_len_bio, label='Succinct DBG (biological)', marker='s', color='darkorange', linestyle='--')
ax2.set_xlabel(f'Read Length in bp (count={READ_COUNT})')
ax2.set_ylabel('Memory Usage (Bytes)')
ax2.set_title('Memory vs. Read Length')
ax2.legend()
ax2.grid(True, linestyle='--', alpha=0.6)

fig.suptitle('Memory Scaling: Standard vs. Succinct DBG  |  Random vs. Biological',
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.show()
