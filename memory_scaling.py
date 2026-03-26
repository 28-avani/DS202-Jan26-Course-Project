class StandardDBG:
    def __init__(self, k):
        self.k = k
        self.nodes = set()
        self.edges = dict()
        self.concat_reads = '' 

    def add_kplus1_mer(self, kp1_mer):
        """
        Input is a string of length k+1.
        u (prefix) and v (suffix) are the nodes (k-mers).
        """
        u = kp1_mer[:-1] # First k chars
        v = kp1_mer[1:]  # Last k chars
        
        self.nodes.add(u)
        self.nodes.add(v)
        
        # The edge is defined by the character that turns u into v
        # u + new_char = v's full string is not needed, just the last char
        new_char = kp1_mer[-1] 
        
        # Adding the edge
        if u not in self.edges:
            self.edges[u] = set()
        self.edges[u].add(new_char) 

        if v not in self.edges:
            self.edges[v] = set()
  

            
    def concat_string(self, reads):
        # Pad the read so there are sentinel characters in between the read
        # I'm using 0...9 as these special characters
        # I'm assuming for this demostration I will not input >10 reads
        assert len(reads) < 10, "The numbers of reads must not be greater than 10"
        assert len(reads) > 0, "There must be atleast 2 reads"

        # Build the concatenated string
        read_0 = ("0" * self.k) + reads[0]
        concat_string = read_0
        for i in range(1, len(reads)):
            concat_string += str(i) + reads[i]
        concat_string += str(len(reads))

        self.concat_reads = concat_string
        return concat_string

    def build_from_reads(self, reads):
        concat_string = self.concat_string(reads)
        # Build a window of (k+1)-mers that slides across the string
        for i in range(len(concat_string) - (self.k + 1) + 1):
                kp1_mer = concat_string[i : i + self.k + 1]
                self.add_kplus1_mer(kp1_mer)


class SuccinctDBG:
    def __init__(self, dbg):
        self.dbg = dbg
        self.k = self.dbg.k
        self.kp1mers = []
        self.W    = []       # list of edge-label tokens  (one per (k+1)-mer)
        self.Node = []    # list of kmers in co-lex order
        self.last = []    # list of 0/1 — 1 marks the last edge of each node group
        self.F    = {'$': None, 'A': None, 'C': None, 'G': None, 'T': None }     # dict: token → first-row index for that last-char group

    # ── Build ──────────────────────────────────────────────
    def kp1mers_from_concat_string(self, concat_string):
        # Build a window of (k+1)-mers that slides across the string
        for i in range(len(concat_string) - (self.k + 1) + 1):
                kp1_mer = concat_string[i : i + self.k + 1]
                self.kp1mers.append(kp1_mer)

    def convert_to_boss(self):
        """
        Converts a StandardDBG into a BOSS SuccinctDBG.

        Each (k+1)-mer is a tuple of tokens. Edges are sorted co-lexicographically:
        first by the *reversed* source k-mer (primary key), then by edge label (tie-
        break). This is the core invariant that makes BOSS navigation work.

        Alphabet ordering (from std_dbg.alphabet) determines all sort keys, so
        sentinels (0 < 1 < … < 9) naturally sort before A < C < G < T.
        """
        # 1. Collect all (k+1)-mers
        self.kp1mers_from_concat_string(self.dbg.concat_reads)

        # Remove duplicates (BOSS only stores unique edges)
        self.kp1mers = list(set(self.kp1mers))

        # 2. Co-lex sort
        self.kp1mers.sort(key=lambda s: (s[:-1][::-1], s[-1]))

        # 3. Build W (list of tokens), Node (list of kmers) and last (list of 0/1)
        for i in range(len(self.kp1mers)):
            self.W.append(self.kp1mers[i][-1])
            self.Node.append(self.kp1mers[i][:-1])
            if (i+1) <= (len(self.kp1mers) - 1):
                same_next = (self.kp1mers[i][:-1] == self.kp1mers[i+1][:-1])  
            else:
                same_next = 0
            self.last.append(0 if same_next else 1)

        # 4. Build F: for each character c, F[c] = first row index whose source
        #    node ends with c.  Because edges are co-lex sorted (primary key =
        #    last char of source node), each character's rows are contiguous.
        #    The last char of the source node of edge e is e[-2].
        self.F['$'] = 0
        for pos in range(1, len(self.kp1mers)):
            diff_before = (self.kp1mers[pos][-2] != self.kp1mers[pos-1][-2]) and (self.kp1mers[pos][-2] in ("A", "C", "G", "T")) 
            if diff_before:
                self.F[self.kp1mers[pos][-2]] = pos

    # ── Rank / Select primitives ──────────────────────────────────────────────

    def rank(self, structure, char, idx):
        """Count of `char` in structure[0..idx] (inclusive)."""
        return structure[:idx + 1].count(char)

    def select(self, structure, char, count):
        """Index of the count-th occurrence of `char` in structure (1-based count)."""
        found = 0
        for i, val in enumerate(structure):
            if val == char:
                found += 1
                if found == count:
                    return i
        return -1   # not found

    # ── BOSS navigation ───────────────────────────────────────────────────────

    def get_node_last_char(self, i):
        """
        Last character of the k-mer at row i, determined by which F-group row i
        falls into. Checks in reverse alphabet order (largest F value first) to
        find the largest F[char] ≤ i.
        """
        for char in reversed(self.alphabet):
            f = self.F.get(char)
            if f is not None and f != len(self.last) and i >= f:
                return char
        return None

    def bwd(self, i):
        """Edge index j such that edge j points INTO node i. [paper eq. 213-214]"""
        char = self.get_node_last_char(i)
        if char is None:
            return -1
        r = self.rank(self.last, 1, i) - self.rank(self.last, 1, self.F[char])
        return self.select(self.W, char, r)

    def fwd(self, j):
        """Node index i that edge j points TO. [paper eq. 214]"""
        char        = self.W[j]
        r           = self.rank(self.W, char, j)
        target_rank = self.rank(self.last, 1, self.F[char]) + r
        return self.select(self.last, 1, target_rank)

    def out_degree(self, v):
        """Number of outgoing edges from node v. [paper Sec. 3.2]"""
        prev_one = -1
        for i in range(v - 1, -1, -1):
            if self.last[i] == 1:
                prev_one = i
                break
        return v - prev_one

# ── Graph Display ─────────────────────────────────────────────────────────────

import matplotlib.pyplot as plt
import networkx as nx

def display_standard_dbg(dbg, title="Standard De Bruijn Graph"):
    """
    Draws the Standard DBG as a directed graph.
    Nodes are k-mers; edges are labeled with the transition character.
    Sentinel characters (digits) are filtered out.
    """
    G = nx.MultiDiGraph()

    for node in dbg.nodes:
        G.add_node(node)

    edge_labels = {}
    for u, chars in dbg.edges.items():
        for char in chars:
            v = u[1:] + char
            if v in dbg.nodes:
                G.add_edge(u, v, label=char)
                edge_labels[(u, v)] = char

    pos = nx.kamada_kawai_layout(G, scale=2)
    plt.figure(figsize=(14, 10))
    plt.title(title)

    nx.draw_networkx_nodes(G, pos, node_size=2000, node_color='lightblue')
    nx.draw_networkx_labels(G, pos, font_size=9)
    nx.draw_networkx_edges(G, pos, arrows=True, arrowsize=20, min_source_margin=25,
                           min_target_margin=25, connectionstyle='arc3,rad=0.2')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=9,
                                 bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.7))

    plt.axis('off')
    plt.tight_layout()
    plt.show()


def display_succinct_dbg(sdbg, title="Succinct De Bruijn Graph (BOSS)"):
    """
    Displays the BOSS arrays as two side-by-side tables:
      - Left:  F array (character -> first row index)
      - Right: Row-indexed table of last, Node, W
    """
    fig, (ax_f, ax_rows) = plt.subplots(1, 2, figsize=(10, 1.5 + len(sdbg.Node) * 0.4))
    fig.suptitle(title, fontsize=13, fontweight='bold')

    # ── F table (left) ───────────────────────────────────────────────────────
    ax_f.axis('off')
    f_keys = list(sdbg.F.keys())
    f_data = [[str(sdbg.F[k])] for k in f_keys]
    f_table = ax_f.table(
        cellText=f_data,
        colLabels=['index'],
        rowLabels=f_keys,
        loc='center',
        cellLoc='center',
    )
    f_table.auto_set_font_size(False)
    f_table.set_fontsize(10)
    f_table.scale(1, 1.5)
    ax_f.set_title('F', fontsize=10, pad=4)

    # ── Row-indexed table (right): last, Node, W ─────────────────────────────
    ax_rows.axis('off')
    n = len(sdbg.Node)
    rows_data = [[sdbg.last[i], sdbg.Node[i], sdbg.W[i]] for i in range(n)]
    row_labels = [str(i) for i in range(n)]
    rows_table = ax_rows.table(
        cellText=rows_data,
        colLabels=['last', 'Node', 'W'],
        rowLabels=row_labels,
        loc='center',
        cellLoc='center',
    )
    rows_table.auto_set_font_size(False)
    rows_table.set_fontsize(10)
    rows_table.scale(1, 1.5)

    plt.tight_layout()
    plt.show()

# ── Memory comparison ─────────────────────────────────────────────────────────

import sys

def compare_memory(standard, succinct):
    # Rough estimate — Python object overhead inflates both numbers equally
    std_size = sys.getsizeof(standard.nodes) + sys.getsizeof(standard.edges)

    suc_size = (sys.getsizeof(succinct.W)
                + sys.getsizeof(succinct.last)
                + sys.getsizeof(succinct.F))

    print(f'--- Memory Comparison ---')
    print(f'Standard DBG approx : {std_size} bytes')
    print(f'Succinct DBG approx : {suc_size} bytes')
    print(f'Reduction : {100 - (suc_size / std_size * 100):.2f}%')



# --- Execution ---
if __name__ == '__main__':
    dbg = StandardDBG(k=3)
    reads = ["TACAC", "TACTC", "GACTC"]
    dbg.build_from_reads(reads)
    print(f'--- Standard DBG ---')
    print(f"Concatenated String: {dbg.concat_string(reads)}")
    print(f"List of Nodes: {dbg.nodes}")
    print(f"List of Nodes and Edges: {dbg.edges}")


    sdbg = SuccinctDBG(dbg)
    sdbg.convert_to_boss()
    print(f'--- Succinct DBG ---')
    print(f"F array: {sdbg.F}")
    print(f"last array: {sdbg.last}")
    print(f"Node array: {sdbg.Node}")
    print(f"W array: {sdbg.W}")

    compare_memory(dbg, sdbg)

    # --- Plotting Graphs ---
    display_standard_dbg(dbg)
    display_succinct_dbg(sdbg)

    #ToDO
    # figure out the order of two same Nodes (TAC, TAC)
    # plot and compare "reduction" vs read size
  