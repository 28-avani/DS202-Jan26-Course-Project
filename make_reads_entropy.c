#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Probability of Base 'A'. We vary this to change Entropy.
// When A=0.25, Entropy is max (2.0). When A=0.9, Entropy is very low.
double a_probs[10] = {0.95, 0.85, 0.75, 0.60, 0.45, 0.35, 0.30, 0.27, 0.26, 0.25};
const char* labels[10] = {"pt1", "pt2", "pt3", "pt4", "pt5", "pt6", "pt7", "pt8", "pt9", "pt10"};

void generate_genome(const char* name, int len, double pA) {
    FILE *f = fopen(name, "w");
    fprintf(f, ">%s\n", name);
    double pOthers = (1.0 - pA) / 3.0;
    for (int i = 0; i < len; i++) {
        double r = (double)rand() / RAND_MAX;
        if (r < pA) fputc('A', f);
        else if (r < pA + pOthers) fputc('C', f);
        else if (r < pA + 2*pOthers) fputc('G', f);
        else fputc('T', f);
        if (i % 80 == 79) fputc('\n', f);
    }
    fclose(f);
}

int main() {
    srand(time(NULL));
    int genome_len = 2000000; // L. acidophilus scale
    for(int i=0; i<10; i++) {
        char fname[32];
        sprintf(fname, "ref_%s.fasta", labels[i]);
        generate_genome(fname, genome_len, a_probs[i]);
        // Simulate reads using your existing wgsim command logic here...
        printf("Generated %s with pA=%.2f\n", fname, a_probs[i]);
    }
    return 0;
}