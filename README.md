# HPRU_AK_Rna_Seq
RNA-Seq data from 2 patient and one control group.

# Fastq_checks.R
Loads the data from the sample annotation file and splits the samples into lanes. Loads the list of fastq files for each lane and performs qc analysis using ShortRead library and produces a report. It also outputs a fasta file with the top 3 (n = 3 default) recurrent sequences, which can be searched for on ncbi-blast to investigate further, the source of those sequences.

## f_getqa = function(files, title, nseq=3)
### DESC:
takes the list of files, title for the output directoy and number of common sequences to output in a fasta file. It performs the qa in serial instead of parallel.