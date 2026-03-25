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
            print(self.kp1mers[i])
            if (i+1) <= (len(self.kp1mers) - 1):
                same_next = (self.kp1mers[i][:-1] == self.kp1mers[i+1][:-1])  
            else:
                same_next = 1
            self.last.append(0 if same_next else 1)

        # 4. Build F: for each character c, F[c] = first row index whose source
        #    node ends with c.  Because edges are co-lex sorted (primary key =
        #    last char of source node), each character's rows are contiguous.
        #    The last char of the source node of edge e is e[-2].
        self.F['$'] = 0
        for pos in range(1, len(self.kp1mers)):
            diff_before = (self.kp1mers[i][:-2] != self.kp1mers[i-1][:-2])  
            if diff_before and self.kp1mers[i][-2] in ("A", "C", "G", "T"):
                self.F[self.kp1mers[i][:-2]] = pos

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



# ── Memory comparison ─────────────────────────────────────────────────────────

import sys

def compare_memory(standard, succinct):
    # Rough estimate — Python object overhead inflates both numbers equally
    std_size = sys.getsizeof(standard.nodes) + sys.getsizeof(standard.edges)
    for u, mask in standard.edges.items():
        std_size += sys.getsizeof(u) + sys.getsizeof(mask)

    suc_size = (sys.getsizeof(succinct.W)
                + sys.getsizeof(succinct.last)
                + sys.getsizeof(succinct.F))

    print(f'--- Memory Comparison ---')
    print(f'Standard DBG approx : {std_size} bytes')
    print(f'Succinct DBG approx : {suc_size} bytes')
    print(f'Reduction           : {100 - (suc_size / std_size * 100):.2f}%')



# --- Execution ---
if __name__ == '__main__':
    dbg = StandardDBG(k=3)
    reads = ["TACAC", "TACTC", "GACTC"]
    dbg.build_from_reads(reads)
    print(dbg.concat_string(reads))
    print(dbg.nodes)
    print(dbg.edges)

    sdbg = SuccinctDBG(dbg)
    sdbg.convert_to_boss()
    print(sdbg.F)
    print(sdbg.last)
    print(sdbg.Node)
    print(sdbg.W)
  