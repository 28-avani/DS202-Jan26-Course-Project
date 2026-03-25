class StandardDBG:
    def __init__(self, k):
        self.k = k
        self.nodes = set()
        self.edges = dict() 

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
        
        return concat_string

    def build_from_reads(self, reads):
        concat_string = self.concat_string(reads)
        # Build a window of (k+1)-mers that slides across the string
        for i in range(len(concat_string) - (self.k + 1) + 1):
                kp1_mer = concat_string[i : i + self.k + 1]
                self.add_kplus1_mer(kp1_mer)


# --- Execution ---
if __name__ == '__main__':
    dbg = StandardDBG(k=3)
    reads = ["TACAC", "TACTC", "GACTC"]
    dbg.build_from_reads(reads)
    print(dbg.concat_string(reads))
    print(dbg.nodes)
    print(dbg.edges)
  