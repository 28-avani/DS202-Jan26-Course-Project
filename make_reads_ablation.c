//Makes reads for MEGAHIT ablation study. Generates a FASTQ file with realistic read lengths, error profiles, and uneven abundance ratios between two input genomes.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_SEQS 4
#define MAX_SEQ_LEN 100000000 // Max 100MB genome per file
#define MAX_READ_LEN 1024

// Input Genomes
const char* FASTA_FILES[MAX_SEQS] = {"S_aureus.fasta"}; // Add more FASTA files as needed
int NUM_FILES = 1;

// Output FASTQ File (Standard MEGAHIT Input)
const char* OUTPUT_FILE = "mock_metagenome_exp1.fq";

// Read Ratios (The "Uneven Depth" Trigger)
// E.g., S. aureus is highly abundant (500,000 reads), L. acidophilus is rare (5,000 reads)
// Array order: [P_aeruginosa, M_tuberculosis, S_aureus, L_acidophilus]
// This gives every genome roughly ~2.5x coverage.
int NUM_READS[1] = {350000};

// Illumina Sequencing Characteristics
double MEAN_READ_LEN = 150.0;
double STD_DEV_LEN = 0.0;

// Noise Parameters (Toggle ERROR_MODE to 0 for perfect reads, 1 for realistic)
int ERROR_MODE = 0; 
double CHIMERA_RATE = 0.00;    // 0% chance of fusing two random sequences
double ERROR_RATE = 0.00;      // 0% substitution sequencing error rate
// ============================================================================

typedef struct {
    char* seq;
    int length;
} Sequence;

Sequence genomes[MAX_SEQS];

// Box-Muller transform to generate normally distributed numbers
double rand_normal(double mean, double stddev) {
    double u1 = ((double)rand() / RAND_MAX);
    double u2 = ((double)rand() / RAND_MAX);
    if (u1 <= 0.0) u1 = 1e-7; 
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return z0 * stddev + mean;
}

// Function to read a standard FASTA file
void read_fasta(const char* filename, int index) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error: Could not open %s. Please ensure the file exists.\n", filename);
        exit(1);
    }

    genomes[index].seq = malloc(MAX_SEQ_LEN);
    genomes[index].length = 0;
    char line[1024];
    
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '>') continue; 
        for (int i = 0; line[i] != '\0'; i++) {
            if (line[i] == 'A' || line[i] == 'C' || line[i] == 'G' || line[i] == 'T' ||
                line[i] == 'a' || line[i] == 'c' || line[i] == 'g' || line[i] == 't') {
                // Convert to uppercase just in case
                genomes[index].seq[genomes[index].length++] = (line[i] & ~32); 
            }
        }
    }
    genomes[index].seq[genomes[index].length] = '\0';
    fclose(file);
    printf("Loaded %s (Length: %d bp)\n", filename, genomes[index].length);
}

// Introduce random substitution errors
void mutate_read(char* read, int len) {
    const char bases[] = "ACGT";
    for (int i = 0; i < len; i++) {
        if (((double)rand() / RAND_MAX) < ERROR_RATE) {
            char original = read[i];
            char mutated;
            do {
                mutated = bases[rand() % 4];
            } while (mutated == original);
            read[i] = mutated;
        }
    }
}

// Grab a random slice from a genome
void get_fragment(int seq_idx, int target_len, char* buffer) {
    if (target_len >= genomes[seq_idx].length) target_len = genomes[seq_idx].length - 1;
    int start_pos = rand() % (genomes[seq_idx].length - target_len);
    strncpy(buffer, genomes[seq_idx].seq + start_pos, target_len);
    buffer[target_len] = '\0';
}

int main() {
    srand(time(NULL));

    printf("--- Generating FASTQ Mock Metagenome for MEGAHIT ---\n");
    
    for (int i = 0; i < NUM_FILES; i++) read_fasta(FASTA_FILES[i], i);

    FILE* out = fopen(OUTPUT_FILE, "w");
    if (!out) return 1;

    int read_counter = 0;
    
    for (int i = 0; i < NUM_FILES; i++) {
        for (int r = 0; r < NUM_READS[i]; r++) {
            
            int read_len = (int)rand_normal(MEAN_READ_LEN, STD_DEV_LEN);
            if (read_len < 50) read_len = 50; 
            if (read_len > MAX_READ_LEN - 1) read_len = MAX_READ_LEN - 1;

            char read_buf[MAX_READ_LEN];

            // Chimera Logic
            if (ERROR_MODE && ((double)rand() / RAND_MAX) < CHIMERA_RATE) {
                int len1 = read_len / 2;
                int len2 = read_len - len1;
                
                // Randomly fuse from either genome
                int seq1_idx = rand() % NUM_FILES;
                int seq2_idx = rand() % NUM_FILES;

                char frag1[MAX_READ_LEN], frag2[MAX_READ_LEN];
                get_fragment(seq1_idx, len1, frag1);
                get_fragment(seq2_idx, len2, frag2);

                strcpy(read_buf, frag1);
                strcat(read_buf, frag2);
            } else {
                get_fragment(i, read_len, read_buf);
            }

            if (ERROR_MODE) mutate_read(read_buf, read_len);

            // Generate Phred+33 Quality String (Standard Illumina representation)
            // 'I' represents a high quality score. 
            char qual_buf[MAX_READ_LEN];
            for(int q=0; q<read_len; q++) qual_buf[q] = 'I';
            qual_buf[read_len] = '\0';

            // Write FASTQ Format
            // Line 1: @ReadID
            // Line 2: Sequence
            // Line 3: +
            // Line 4: Quality Scores
            fprintf(out, "@Read_%d_Source_%d\n%s\n+\n%s\n", read_counter++, i, read_buf, qual_buf);
        }
    }

    fclose(out);
    for(int i=0; i<NUM_FILES; i++) free(genomes[i].seq);

    printf("Successfully generated %d reads to %s\n", read_counter, OUTPUT_FILE);
    return 0;
}