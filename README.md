# HPRU_AK_Rna_Seq
RNA-Seq data from 2 patient and one control group.

# Fastq_checks.R
Loads the data from the sample annotation file and splits the samples into lanes. Loads the list of fastq files for each lane and performs qc analysis using ShortRead library and produces a report. It also outputs a fasta file with the top 3 (n = 3 default) recurrent sequences, which can be searched for on ncbi-blast to investigate further, the source of those sequences.

## f_getqa = function(files, title, nseq=3)
### DESC:
takes the list of files, title for the output directoy and number of common sequences to output in a fasta file. It performs the qa in serial instead of parallel.  

# analysis01.R
Details of the analysis are here: https://www.evernote.com/shard/s288/nl/38698211/1af31cab-e28c-480b-a91a-b8cee95eaa3b  
DESeq2 is used to call DE genes in the 3 comparisons, out of which the 2 comparisons are used for further analysis. The DE genes are grouped into 3 possible classes based on 2^2-1 combinations. The graphs for each set of genes are built separately, by choosing the apprpriate variable manually. The data is stabalized before calculating correlations.  

# expressed.R
Normalizes the data using DESeq2::estimateSizeFactorsForMatrix. Creates a vector of mean expressions for each group. Plots the log distribution, approximating it as normally distributed. Uses cut to assign groups to each gene based on quantiles of the expression vector. Saves the data.

# find_clusters.R
uses the msd data, removes missing values, scales and clusters the data to find subgroups.


