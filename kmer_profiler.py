class KMerProfiler:
    def __init__(self, sequence, k):
        self.sequence = sequence.upper()
        self.k = k
        self.kmer_counts = {}  # Our Hash Table structure

    def profile_frequencies(self):
        """Slides a window of size K across the sequence in O(N) time."""
        n = len(self.sequence)
        
        # Edge case: sequence must be longer than or equal to window size K
        if n < self.k or self.k <= 0:
            print("[ERROR] Invalid sequence length or K-value.")
            return self.kmer_counts

        # Sliding Window logic
        for i in range(n - self.k + 1):
            kmer = self.sequence[i : i + self.k]
            
            # Update frequency in our Hash Table
            if kmer in self.kmer_counts:
                self.kmer_counts[kmer] += 1
            else:
                self.kmer_counts[kmer] = 1
                
        return self.kmer_counts

    def get_most_frequent(self, top_n=3):
        """Finds the top N repeating genomic signatures."""
        # Sort hash map entries by count in descending order
        sorted_kmers = sorted(self.kmer_counts.items(), key=lambda item: item[1], reverse=True)
        return sorted_kmers[:top_n]


# --- Simulation Run ---
if __name__ == "__main__":
    # Example raw DNA sequence with repeating motifs
    raw_dna = "ACGTGACGTGACGTGATTCGATCGAACGTG"
    k_size = 4
    
    print("=== Genomic K-Mer Frequency Profiler Active ===")
    print(f"Analyzing Sequence Length: {len(raw_dna)} base pairs")
    print(f"Target Substring Window Size (K): {k_size}\n")
    
    profiler = KMerProfiler(raw_dna, k_size)
    all_counts = profiler.profile_frequencies()
    top_signatures = profiler.get_most_frequent(3)
    
    print("=== Processing Complete ===")
    print(f"Total Unique {k_size}-mers Found: {len(all_counts)}")
    print("\nTop 3 Most Frequent Genomic Signatures:")
    for kmer, count in top_signatures:
        print(f"  -> Signature [{kmer}]: Repeated {count} times")