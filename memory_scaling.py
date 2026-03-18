# Trying to implement a graph and table output, we can explain succinct dBGs with it!
# This implementation is from the BOSS paper


class StandardDBG:
    def __init__(self, k):
        self.k = k
        self.nodes = set()
        self.edges = {} # k-mer -> bitmask of outgoing bases (forming a k+1 mer)

    def add_kplus1_mer(self, kp1_mer):
        """
        Input is a string of length k+1.
        u (prefix) and v (suffix) are the nodes (k-mers).
        """
        u = kp1_mer[:-1] # First k chars
        v = kp1_mer[1:]  # Last k chars
        
        char_to_bit = {'A':0, 'C':1, 'G':2, 'T':3}
        
        self.nodes.add(u)
        self.nodes.add(v)
        
        # The edge is defined by the character that turns u into v
        # u + new_char = v's full string is not needed, just the last char
        new_char = kp1_mer[-1] 
        
        mask = self.edges.get(u, 0)
        self.edges[u] = mask | (1 << char_to_bit[new_char])

    def build_from_reads(self, reads):
        for read in reads:
            # To get nodes of length k and edges of length k+1,
            # we must slide a window of size k+1.
            for i in range(len(read) - (self.k + 1) + 1):
                kp1_mer = read[i : i + self.k + 1]
                self.add_kplus1_mer(kp1_mer)

    def build_from_reads_with_padding(self, reads):
        for read in reads:
            # Pad the read so the first character is an edge from $$$
            padded_read = ("$" * self.k) + read + "$"
            for i in range(len(padded_read) - (self.k + 1) + 1):
                kp1_mer = padded_read[i : i + self.k + 1]
                self.add_kplus1_mer(kp1_mer)

    def get_outgoing(self, node):
        """Helper to see which bases follow a node"""
        mask = self.edges.get(node, 0)
        bases = "ACGT"
        return [bases[i] for i in range(4) if (mask & (1 << i))]
        
        
class SuccinctDBG:
    def __init__(self, k, W, last, F):
        self.k = k
        self.W = W            # String of edge labels (including A- markers)
        self.last = last      # Bitarray or list of 0s and 1s
        self.F = F            # Dictionary: char -> cumulative count
        self.alphabet = "ACGT"

    # --- Basic Succinct Operations ---
    def rank(self, structure, char, idx):
        """Number of occurrences of char in structure up to idx (inclusive)"""
        return structure[:idx + 1].count(char)

    def select(self, structure, char, count):
        """Position of the n-th occurrence of char"""
        found = 0
        for i, val in enumerate(structure):
            if val == char:
                found += 1
                if found == count:
                    return i
        return -1

    # --- Paper's Navigation Functions ---
    def bwd(self, i):
        """Mapping from node index i to edge index j [cite: 214]"""
        # C(i) is the last character of the node at index i [cite: 201]
        char = self.get_node_last_char(i)
        
        # r = rank_1(last, i) - rank_1(last, F[char]) [cite: 213]
        r = self.rank(self.last, 1, i) - self.rank(self.last, 1, self.F[char])
        
        # j = select_c(W, r) [cite: 213]
        return self.select(self.W, char, r)

    def fwd(self, j):
        """Mapping from edge index j to node index i [cite: 214]"""
        char = self.W[j].replace('-', '') # Handle A- notation [cite: 138]
        r = self.rank(self.W, self.W[j], j)
        
        # i = select_1(last, rank_1(last, F[c]) + r) [cite: 214]
        target_rank = self.rank(self.last, 1, self.F[char]) + r
        return self.select(self.last, 1, target_rank)

    def get_node_last_char(self, i):
        """Finds which character group node i belongs to via F array [cite: 201]"""
        sorted_chars = sorted(self.F.keys(), key=lambda x: self.F[x], reverse=True)
        for char in sorted_chars:
            if i >= self.F[char]:
                return char
        return None

    def out_degree(self, v):
        """Returns number of outgoing edges from node v [cite: 132, 220]"""
        # v is an index where last[v] == 1
        # outdegree = v - pred_1(last, v-1) [cite: 220]
        prev_one = -1
        for i in range(v - 1, -1, -1):
            if self.last[i] == 1:
                prev_one = i
                break
        return v - prev_one

    def display(self):
        """Prints the BOSS structure in the format of Fig 2."""
        rows = []
        for i in range(len(self.W)):
            # Determine the character group for the 'F' column label
            char_group = self.get_node_last_char(i)
            
            # Reconstruction of the 'Node' string is for visualization only
            # In a real succinct structure, we wouldn't do this!
            node_label = self.reconstruct_node_label(i)
            
            rows.append({
                "i": i,
                "F (label)": char_group,
                "last": self.last[i],
                "Node": node_label,
                "W": self.W[i]
            })
        
        # Displaying with Pandas for that clean 'Paper' look
        df = pd.DataFrame(rows)
        print("\n--- Succinct BOSS Representation (Fig 2) ---")
        print(df.to_string(index=False))

    def reconstruct_node_label(self, i):
        """Back-traverses from index i to find the k-mer it represents."""
        curr = i
        label = ""
        for _ in range(self.k):
            char = self.get_node_last_char(curr)
            label = char + label
            curr = self.bwd(curr)
        return label

    

import sys

def compare_memory(standard, succinct):
    # Rough estimation of memory (Python object overhead makes this high for both)
    std_size = sys.getsizeof(standard.nodes) + sys.getsizeof(standard.edges)
    for k, v in standard.edges.items():
        std_size += sys.getsizeof(k) + sys.getsizeof(v)
        
    suc_size = sys.getsizeof(succinct.W) + sys.getsizeof(succinct.last) + sys.getsizeof(succinct.F)
    
    print(f"--- Memory Comparison ---")
    print(f"Standard DBG approx: {std_size} bytes")
    print(f"Succinct DBG approx: {suc_size} bytes")
    print(f"Reduction: {100 - (suc_size/std_size*100):.2f}%")

def convert_to_boss(std_dbg):
    """
    Converts StandardDBG to BOSS-style SuccinctDBG.
    Sorting is done by u[::-1] (co-lexicographical).
    """
    edges_list = []
    # Extract all k+1-mers from the standard graph
    for u, mask in std_dbg.edges.items():
        for i, char in enumerate("ACGT"):
            if mask & (1 << i):
                edges_list.append(u + char)

    # Sort by the reverse of the prefix (u), then by the edge label (c)
    # This is the core requirement for the BOSS navigation properties
    edges_list.sort(key=lambda x: (x[:-1][::-1], x[-1]))

    W = ""
    last = []
    
    for i in range(len(edges_list)):
        current_edge = edges_list[i]
        prefix = current_edge[:-1]
        label = current_edge[-1]
        
        W += label
        
        # If the next edge has the same prefix, this is not the last edge for this node
        if i + 1 < len(edges_list) and edges_list[i+1][:-1] == prefix:
            last.append(0)
        else:
            last.append(1)

    # Build F: maps char to the rank of the first occurrence in the sorted edge-node list
    # For a simple version, we sort all incoming labels of all nodes
    all_incoming = sorted([e[:-1][-1] for e in edges_list if len(e[:-1]) > 0])
    F = {}
    for char in "ACGT":
        try:
            F[char] = all_incoming.index(char)
        except ValueError:
            F[char] = len(all_incoming)

    return SuccinctDBG(std_dbg.k, W, last, F)

def export_to_graphviz(std_dbg, filename="dbg_graph"):
    """
    Generates a DOT file representing the Standard DBG.
    To render: install graphviz and run 'dot -Tpng dbg_graph.dot -o graph.png'
    """
    with open(f"{filename}.dot", "w") as f:
        f.write("digraph StandardDBG {\n")
        # Styling to match Fig. 1
        f.write('  node [shape=box, fontname="Arial", fontsize=10];\n')
        f.write('  edge [fontname="Arial", fontsize=10];\n')
        f.write('  rankdir=LR;\n') # Left to Right layout

        for u, mask in std_dbg.edges.items():
            bases = "ACGT"
            for i in range(4):
                if mask & (1 << i):
                    char = bases[i]
                    # In a DBG, the next node v is the suffix of (u + char)
                    v = u[1:] + char
                    f.write(f'  "{u}" -> "{v}" [label=" {char}"];\n')
        
        f.write("}\n")
    print(f"Successfully exported to {filename}.dot")

# --- Execution ---
if __name__ == '__main__':
    dbg = StandardDBG(k=3)
    # Adding more data to see a larger difference
    reads = ["TACAC", "TACTC", "GACTC"]
    dbg.build_from_reads(reads)

    # Convert
    sdbg = convert_to_boss(dbg)

    print(f"Nodes in Standard: {len(dbg.nodes)}")
    print(f"Edges in Succinct (W): {sdbg.W}")
    print(f"Last vector: {sdbg.last}")
    
    compare_memory(dbg, sdbg)