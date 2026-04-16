This is a readme file describing the overall structure of the files here. 

/images: this folder contains some of the plots which were used within the final presentation.
/results: this contains the megahit output log and all raw data from all the experiments.
/results_entropy: results for the various levels of information entropy in the generated sequences

similar for /results_memory_scaling and /results_taxonomic_distance

/sequences: contains the fasta references downloaded from the NCBI database that were used to generate the mock metagenomes for analysis.

other .fasta or.fq files L_monocytogenes.fasta, S_aureus_8325.fasta, etc. were used in the taxonomic distance experiment or others.

The .png files are some of the plots used in the presentations.

/build_random_genome.py is the python script used to build random genomes (random sequences of ATGCs)

the /eval_megahit1.py, eval_entropy1.py and eval_megahit2.py, eval_taxonomic_distance.py are the python scripts used to generate the plots and to evaluate and analyse the results of the megahit code fort various experiments.

/googleColab.ipynb: The megahit code was primarily run on google colab because it was not compatible with macbook. Some of the commands used to run the megahit codes are there in this notebook.

/make_reads_ablation.c /make_reads_entropy.c, make_reads_taxonomic_distance.c etc. taka reference genomes from the respoective locations and make the number of reads with specified mutations and erros. The relative abundance of reads that best satisfies the test condition was obtained by trial and error.

/memory_scaling.py: script for the memory scaling experiment
/memory_scaling_plots.py: creates plots from the results obtained from memory_scaling.py

/mix_high.fq, /mix_med.fq, /mix_low.fq are the reads obtained for the taxonomic distance experiment on which megahit was run.

/mock_metagenome: this contains the abundance disparity metagenome from the 4 different species, on which the first few experiments were run, using the megahit code.

readme.md contains the overall idea of the project.

the ref_pt1.fasta to ref_pt10.fasta are reads from generated genomes with varying levels of information entropy, on which the megahit code was run.

