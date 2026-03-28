import time
import random
import matplotlib.pyplot as plt
import sys
from memory_scaling_compressed import StandardDBG, SuccinctDBG

def generate_random_reads(count, length):
    alphabet = "ACGT"
    return ["".join(random.choice(alphabet) for _ in range(length)) for _ in range(count)]

K = 13

def measure(reads):
    sd = StandardDBG(K)
    sd.build_from_reads(reads)
    ss = SuccinctDBG(sd)
    ss.convert_to_boss()
    std = sd.nodes.nbytes + sd.edge_labels.nbytes + sd.edge_offsets.nbytes
    suc = ss.W.nbytes + ss.last.nbytes
    return std, suc

def run_benchmark_vs_count(read_counts, read_length):
    std_mem, suc_mem = [], []
    for count in read_counts:
        reads = generate_random_reads(count, read_length)
        std, suc = measure(reads)
        std_mem.append(std)
        suc_mem.append(suc)
    return std_mem, suc_mem

def run_benchmark_vs_length(read_lengths, read_count):
    std_mem, suc_mem = [], []
    for length in read_lengths:
        reads = generate_random_reads(read_count, length)
        std, suc = measure(reads)
        std_mem.append(std)
        suc_mem.append(suc)
    return std_mem, suc_mem

# Parameters
counts  = [2, 4, 8, 16, 32, 64, 128, 256]
lengths = [20, 40, 60, 80, 100, 150, 200, 300]

std_vs_count, suc_vs_count = run_benchmark_vs_count(counts, read_length=50)
std_vs_len,   suc_vs_len   = run_benchmark_vs_length(lengths, read_count=64)

# Plotting
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

ax1.plot(counts, std_vs_count, label='Standard DBG', marker='o')
ax1.plot(counts, suc_vs_count, label='Succinct DBG (BOSS)', marker='s')
ax1.set_xlabel('Number of Reads (Length=50)')
ax1.set_ylabel('Memory Usage (Bytes)')
ax1.set_title('Memory vs. Read Count')
ax1.legend()
ax1.grid(True, linestyle='--')

ax2.plot(lengths, std_vs_len, label='Standard DBG', marker='o')
ax2.plot(lengths, suc_vs_len, label='Succinct DBG (BOSS)', marker='s')
ax2.set_xlabel('Read Length (Count=64)')
ax2.set_ylabel('Memory Usage (Bytes)')
ax2.set_title('Memory vs. Read Length')
ax2.legend()
ax2.grid(True, linestyle='--')

fig.suptitle('Memory Scaling: Standard vs. Succinct DBG', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.show()