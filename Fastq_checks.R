# Name: Fastq_checks.R
# Auth: u.niazi@imperial.ac.uk
# Date: 4/12/15
# Desc: create a report for the fastq files



source('Header.R')
library(ShortRead)
#source('../CFastqQuality/CFastqQuality.R')

dfAnnotations = read.csv('Data_external/Sample_information.csv', header=T)

i = nrow(dfAnnotations)
r = 1
nc = ncol(dfAnnotations)
while(r < i){
  dfAnnotations[r+1, 2:nc] = dfAnnotations[r, 2:nc]
  r = r+2
}

dir.create('Results')

summary(dfAnnotations)

f_getqa = function(files, title, nseq=3){
  q = qa(files, BPPARAM=SerialParam())
  report(q, dest=paste('Results/', title, sep=''))
  seq = q[['frequentSequences']]
  # get nseq top sequences from each lane
  seq = split(seq, seq$lane)
  seq.2 = lapply(seq, function(x) as.character(x[1:nseq,'sequence']))
  seq.2 = DNAStringSetList(seq.2)
  Biostrings::writeXStringSet(unlist(seq.2), filepath = paste('Results/', title, '/frequent_seq.fasta', sep=''))
  rm(q)
}


## check samples by factor
fGroup = dfAnnotations$LANE
table(fGroup)

csFiles = paste('Data_external/RNASeq/20151110/FASTQ/', dfAnnotations[fGroup == 1,'Files'], sep='')
f_getqa(csFiles, paste('lane', 1, sep=''))
gc()


csFiles = paste('Data_external/RNASeq/20151110/FASTQ/', dfAnnotations[fGroup == 2,'Files'], sep='')
f_getqa(csFiles, paste('lane', 2, sep=''))
gc()

csFiles = paste('Data_external/RNASeq/20151110/FASTQ/', dfAnnotations[fGroup == 3,'Files'], sep='')
f_getqa(csFiles, paste('lane', 3, sep=''))
gc()

csFiles = paste('Data_external/RNASeq/20151110/FASTQ/', dfAnnotations[fGroup == 4,'Files'], sep='')
f_getqa(csFiles, paste('lane', 4, sep=''))
gc()


# ul = unique(lanes)
# 
# mclapply(seq_along(ul), function(i) {
#   csFiles = paste('Data_external/RNASeq/20151110/FASTQ/', dfAnnotations[lanes == ul[i],'Files'], sep='')
#   qa = f_getqa(csFiles, paste('lane', ul[i], sep=''))
#   f_writeFrequent(qa, paste('lane', ul[i], sep=''))
# })
# 
# for(i in seq_along(ul)){
#   csFiles = paste('Data_external/RNASeq/20151110/FASTQ/', dfAnnotations[lanes == ul[i],'Files'], sep='')
#   qa = f_getqa(csFiles, paste('lane', ul[i], sep=''))
#   f_writeFrequent(qa, paste('lane', ul[i], sep=''))
# }
# 
# f_writeFrequent = function(q, title, n = 3){
#   seq = q[['frequentSequences']]
#   # get n top sequences from each lane
#   seq = split(seq, seq$lane)
#   seq.2 = lapply(seq, function(x) as.character(x[1:n,'sequence']))
#   seq.2 = DNAStringSetList(seq.2)
#   Biostrings::writeXStringSet(unlist(seq.2), filepath = paste('Results/', title, '/frequent_seq.fasta', sep=''))
# }





