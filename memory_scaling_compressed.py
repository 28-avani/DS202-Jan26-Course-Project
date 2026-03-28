import numpy as np

# Alphabet encoding: map each character to a compact uint8 integer.
# Uppercase = normal edge labels; lowercase = redundant (A-) edge labels.
CHAR_TO_UINT8 = {'$': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4,
                          'a': 5, 'c': 6, 'g': 7, 't': 8}
UINT8_TO_CHAR = {v: k for k, v in CHAR_TO_UINT8.items()}

class StandardDBG:
    def __init__(self, k):
        self.k = k
        self.nodes        = np.empty((0, k), dtype=np.uint8)  # 2D array: each row is an encoded k-mer
        self._nodes_temp  = set()                             # temporary set used during build
        # CSR edge structure (populated by build_from_reads)
        self.edge_labels  = np.array([], dtype=np.uint8)  # flat array of all outgoing edge labels
        self.edge_offsets = np.array([0], dtype=np.uint32) # edge_offsets[i]:edge_offsets[i+1] → node i's labels
        self._edges_temp  = dict()                         # temporary dict used during build
        self.concat_reads = ''

    def add_kplus1_mer(self, kp1_mer):
        """
        Input is a string of length k+1.
        u (prefix) and v (suffix) are the nodes (k-mers).
        """
        u = kp1_mer[:-1] # First k chars
        v = kp1_mer[1:]  # Last k chars
        
        self._nodes_temp.add(tuple(CHAR_TO_UINT8[c] for c in u))
        self._nodes_temp.add(tuple(CHAR_TO_UINT8[c] for c in v))
        
        # The edge is defined by the character that turns u into v
        # u + new_char = v's full string is not needed, just the last char
        new_char = kp1_mer[-1] 
        
        # Accumulate edge label into temporary dict (label stored as compact uint8 index)
        if u not in self._edges_temp:
            self._edges_temp[u] = set()
        self._edges_temp[u].add(CHAR_TO_UINT8[new_char])

        if v not in self._edges_temp:
            self._edges_temp[v] = set()
  

            
    def concat_string(self, reads):
        # Pad the read so there are sentinel characters in between the read
        # Using '$' as the sentinel character (sorts before A/C/G/T)
        assert len(reads) > 0, "There must be atleast 2 reads"

        # Build the concatenated string
        read_0 = ("$" * self.k) + reads[0]
        concat_string = read_0
        for i in range(1, len(reads)):
            concat_string += "$" + reads[i]
        concat_string += "$"

        self.concat_reads = concat_string
        return concat_string

    def build_from_reads(self, reads):
        concat_string = self.concat_string(reads)
        # Build a window of (k+1)-mers that slides across the string
        for i in range(len(concat_string) - (self.k + 1) + 1):
                kp1_mer = concat_string[i : i + self.k + 1]
                self.add_kplus1_mer(kp1_mer)
        # Finalise nodes: convert unique encoded k-mers to a sorted 2D uint8 array
        self.nodes = np.array(sorted(self._nodes_temp), dtype=np.uint8)
        del self._nodes_temp

        # Finalise edges: build CSR from _edges_temp, ordered by sorted node index
        labels_list = []
        offsets_list = [0]
        for node_row in self.nodes:
            node_str = ''.join(UINT8_TO_CHAR[int(c)] for c in node_row)
            for label in sorted(self._edges_temp.get(node_str, set())):
                labels_list.append(label)
            offsets_list.append(len(labels_list))
        self.edge_labels  = np.array(labels_list, dtype=np.uint8)
        self.edge_offsets = np.array(offsets_list, dtype=np.uint32)
        del self._edges_temp


class SuccinctDBG:
    def __init__(self, dbg):
        self.dbg = dbg
        self.k = self.dbg.k
        self.kp1mers = []
        self.W        = np.array([], dtype=np.uint8)  # uint8 edge-label array (one per (k+1)-mer)
        self.Node     = []    # list of kmers in co-lex order
        self.last     = np.array([], dtype=np.uint8)  # packed bitvector (np.packbits)
        self.last_len = 0                              # number of valid bits in last
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

        # 3. Build W (uint8 array), Node (list of kmers) and last (bit-vector)
        # Track labels used for the CURRENT suffix group
        current_suffix = None
        seen_labels_in_suffix = set()
        w_list    = []
        last_list = []

        for i in range(len(self.kp1mers)):
            curr = self.kp1mers[i]
            edge_label = curr[-1]
            node_suffix = curr[1:-1] # The (k-1) suffix

            # 2. Check if we entered a NEW suffix group
            if node_suffix != current_suffix:
                current_suffix = node_suffix
                seen_labels_in_suffix = set() # Reset for the new destination node group

            # 3. Redundancy Check: Have we seen this edge label for this suffix?
            is_redundant = False
            if edge_label in seen_labels_in_suffix:
                is_redundant = True
            else:
                seen_labels_in_suffix.add(edge_label)

            # 4. Store the marked/unmarked character as compact uint8 index
            # Redundant edges use the lowercase character (A- notation)
            final_w = edge_label.lower() if is_redundant else edge_label
            w_list.append(CHAR_TO_UINT8[final_w])

            self.Node.append(curr[:-1])

            # 5. Build 'last' bitvector (standard BOSS logic)
            # 1 if this is the last outgoing edge for this specific source node
            if (i + 1) < len(self.kp1mers):
                is_last = (curr[:-1] != self.kp1mers[i+1][:-1])
            else:
                is_last = True
            last_list.append(1 if is_last else 0)

        # Convert accumulated arrays: W as compact uint8, last as packed bits
        self.W        = np.array(w_list, dtype=np.uint8)
        last_bool     = np.array(last_list, dtype=np.uint8)
        self.last_len = len(last_bool)
        self.last     = np.packbits(last_bool)

        # 4. Build F: for each character c, F[c] = first row index whose source
        #    node ends with c.  Because edges are co-lex sorted (primary key =
        #    last char of source node), each character's rows are contiguous.
        #    The last char of the source node of edge e is e[-2].
        self.F['$'] = 0
        for pos in range(1, len(self.kp1mers)):
            diff_before = (self.kp1mers[pos][-2] != self.kp1mers[pos-1][-2]) and (self.kp1mers[pos][-2] in ("A", "C", "G", "T"))
            if diff_before:
                self.F[self.kp1mers[pos][-2]] = pos

    @property
    def last_unpacked(self):
        """Unpack the stored bit-vector back to a uint8 array of 0/1."""
        return np.unpackbits(self.last)[:self.last_len]

    # ── Rank / Select primitives ──────────────────────────────────────────────

    def rank(self, structure, char, idx):
        """Count of `char` in structure[0..idx] (inclusive)."""
        if isinstance(structure, np.ndarray):
            return int(np.sum(structure[:idx + 1] == char))
        return structure[:idx + 1].count(char)

    def select(self, structure, char, count):
        """Index of the count-th occurrence of `char` in structure (1-based count)."""
        if isinstance(structure, np.ndarray):
            indices = np.where(structure == char)[0]
            return int(indices[count - 1]) if count <= len(indices) else -1
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
            if f is not None and f != self.last_len and i >= f:
                return char
        return None

    def bwd(self, i):
        """Edge index j such that edge j points INTO node i. [paper eq. 213-214]"""
        char = self.get_node_last_char(i)
        if char is None:
            return -1
        last = self.last_unpacked
        r = self.rank(last, 1, i) - self.rank(last, 1, self.F[char])
        return self.select(self.W, np.uint8(CHAR_TO_UINT8[char]), r)

    def fwd(self, j):
        """Node index i that edge j points TO. [paper eq. 214]"""
        char_uint8  = self.W[j]
        char        = UINT8_TO_CHAR[int(char_uint8)]
        r           = self.rank(self.W, char_uint8, j)
        last = self.last_unpacked
        target_rank = self.rank(last, 1, self.F[char]) + r
        return self.select(last, 1, target_rank)

    def out_degree(self, v):
        """Number of outgoing edges from node v. [paper Sec. 3.2]"""
        prev_one = -1
        last = self.last_unpacked
        for i in range(v - 1, -1, -1):
            if last[i] == 1:
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

    # Decode each uint8 row back to a string k-mer for display
    nodes_decoded = [''.join(UINT8_TO_CHAR[int(c)] for c in row) for row in dbg.nodes]
    nodes_set = {tuple(row) for row in dbg.nodes}   # for O(1) membership checks
    for node_str in nodes_decoded:
        G.add_node(node_str)

    edge_labels = {}
    for i, u in enumerate(nodes_decoded):
        start, end = dbg.edge_offsets[i], dbg.edge_offsets[i + 1]
        for char in dbg.edge_labels[start:end]:
            char_str = UINT8_TO_CHAR[int(char)]
            v = u[1:] + char_str
            v_key = tuple(CHAR_TO_UINT8[c] for c in v)
            if v_key in nodes_set:
                G.add_edge(u, v, label=char_str)
                edge_labels[(u, v)] = char_str

    pos = nx.kamada_kawai_layout(G, scale=2)
    plt.figure(figsize=(14, 10))
    plt.title(title)

    node_labels = {n: n.replace('$', r'\$') for n in G.nodes()}
    escaped_edge_labels = {k: v.replace('$', r'\$') for k, v in edge_labels.items()}

    nx.draw_networkx_nodes(G, pos, node_size=2000, node_color='lightblue')
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=9)
    nx.draw_networkx_edges(G, pos, arrows=True, arrowsize=20, min_source_margin=25,
                           min_target_margin=25, connectionstyle='arc3,rad=0.2')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=escaped_edge_labels, font_size=9,
                                 bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.7))

    plt.axis('off')
    plt.tight_layout()
    plt.savefig('Standard_DBG_graph.png') 
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
        rowLabels=[k.replace('$', r'\$') for k in f_keys],
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
    last = sdbg.last_unpacked
    rows_data = [[last[i], sdbg.Node[i].replace('$', r'\$'), UINT8_TO_CHAR[int(sdbg.W[i])].replace('$', r'\$')] for i in range(n)]
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
    plt.savefig('Succinct_DBG_table.png')
    plt.show()

# ── Memory comparison ─────────────────────────────────────────────────────────

import sys

def compare_memory(standard, succinct):
    # Rough estimate — Python object overhead inflates both numbers equally
    std_size = (standard.nodes.nbytes
                + standard.edge_labels.nbytes
                + standard.edge_offsets.nbytes)

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
    print(f"List of Nodes (uint8):\n{dbg.nodes}")
    print(f"List of Nodes (decoded): {[''.join(UINT8_TO_CHAR[int(c)] for c in row) for row in dbg.nodes]}")
    print(f"Edge labels (CSR flat):   {dbg.edge_labels}")
    print(f"Edge offsets (CSR):       {dbg.edge_offsets}")


    sdbg = SuccinctDBG(dbg)
    sdbg.convert_to_boss()
    print(f'--- Succinct DBG ---')
    print(f"F array: {sdbg.F}")
    print(f"last array (packed): {sdbg.last}")
    print(f"last array (unpacked): {sdbg.last_unpacked}")
    print(f"Node array: {sdbg.Node}")
    print(f"W array (uint8): {sdbg.W}")
    print(f"W array (chars): {[UINT8_TO_CHAR[int(x)] for x in sdbg.W]}")

    compare_memory(dbg, sdbg)

    # --- Plotting Graphs ---
    display_standard_dbg(dbg)
    display_succinct_dbg(sdbg)

