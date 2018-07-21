library(Biostrings)



#--------import sequences into individual dataframes to be subset later for alignment

#clinical isolate
clin = readDNAStringSet("clin.fasta")
seq_name = names(clin)
sequence = paste(clin)
clin_df <- data.frame(seq_name, sequence)

#hottub isolate
clin = readDNAStringSet("hotb.fasta")
seq_name = names(clin)
sequence = paste(clin)
hotb_df <- data.frame(seq_name, sequence)

#43290 isolate
clin = readDNAStringSet("43290.fasta")
seq_name = names(clin)
sequence = paste(clin)
ATCC_df <- data.frame(seq_name, sequence)

#Lens isolate
clin = readDNAStringSet("lens.fasta")
seq_name = names(clin)
sequence = paste(clin)
lens_df <- data.frame(seq_name, sequence)

#philly isolate
clin = readDNAStringSet("phil.fasta")
seq_name = names(clin)
sequence = paste(clin)
phil_df <- data.frame(seq_name, sequence)

#TB isolate
clin = readDNAStringSet("thun.fasta")
seq_name = names(clin)
sequence = paste(clin)
thun_df <- data.frame(seq_name, sequence)
#---------------------




#---Reading in filtered groups.txt file to subset the dataframes of each species for the core genes
orthologs<-read.table('groups_filter.txt')
#Subset ATCC DF for only those in the orthologs table
ATCC_df_subset<-ATCC_df[ATCC_df$seq_name %in% orthologs$V2,]
ATCC_df_subset$sequence<-as.character(ATCC_df_subset$sequence)
#Subset clin DF for only those in the orthologs table
clin_df_subset<-clin_df[clin_df$seq_name %in% orthologs$V3,]
clin_df_subset$sequence<-as.character(clin_df_subset$sequence)
#Subset hottub DF for only those in the orthologs table
hotb_df_subset<-hotb_df[hotb_df$seq_name %in% orthologs$V4,]
hotb_df_subset$sequence<-as.character(hotb_df_subset$sequence)
#Subset Lens DF for only those in the orthologs table
lens_df_subset<-lens_df[lens_df$seq_name %in% orthologs$V5,]
lens_df_subset$sequence<-as.character(lens_df_subset$sequence)
#Subset philly DF for only those in the orthologs table
phil_df_subset<-phil_df[phil_df$seq_name %in% orthologs$V6,]
phil_df_subset$sequence<-as.character(phil_df_subset$sequence)
#Subset TB DF for only those in the orthologs table
thun_df_subset<-thun_df[thun_df$seq_name %in% orthologs$V7,]
thun_df_subset$sequence<-as.character(thun_df_subset$sequence)





library(msa)
install.packages('rphast')
library(rphast)





#the problem is that you cant compare each row, but must compare the sequneces listed in the orthologs table 

#modification is to instead of subset i in teh my_sqeuences call, we subset the sequence that is mentioned in teh first row of orthologs 
#as.character(ATCC_df_subset$seq_name) == as.character(orthologs$V2[i])


#---------------------------basis for teh for loop that will perform all teh msas 
my_sequences<-c(ATCC_df_subset$sequence[1], clin_df_subset$sequence[1], hotb_df_subset$sequence[1], lens_df_subset$sequence[1], phil_df_subset$sequence[1], thun_df_subset$sequence[1])
test<-msa::msa(my_sequences, type = 'dna', 'Muscle')
test_convert<-msaConvert(test, type="seqinr::alignment")
write.fasta(test_convert$seq[1], names=strains, file.out='test.fasta')



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#function to iterate over all sequences in teh dataframe and concatenate alignments together below: 
#i<-1
concat_align<-c('','','','','','')
for (i in 1:length(ATCC_df_subset$sequence)){
  my_sequences<-c(ATCC_df_subset$sequence[as.character(ATCC_df_subset$seq_name) == as.character(orthologs$V2[i])], clin_df_subset$sequence[as.character(clin_df_subset$seq_name) == as.character(orthologs$V3[i])], hotb_df_subset$sequence[as.character(hotb_df_subset$seq_name) == as.character(orthologs$V4[i])], lens_df_subset$sequence[as.character(lens_df_subset$seq_name) == as.character(orthologs$V5[i])], phil_df_subset$sequence[as.character(phil_df_subset$seq_name) == as.character(orthologs$V6[i])], thun_df_subset$sequence[as.character(thun_df_subset$seq_name) == as.character(orthologs$V7[i])])
  alignment<-msa::msa(my_sequences, type='dna', order = 'input', 'Muscle')
  alignment_convert<-msaConvert(alignment, type="seqinr::alignment") #as alignments are made they can be pasted together!!!!!!! like this : merge<-paste0(test_convert$seq, test_convert$seq)
  concat_align<-paste0(concat_align, alignment_convert$seq)
}





#save teh concatenated alignments to file 
save(concat_align, file='concat_alignment_NEW_MUSCLE.RData')
#write out each individual concatenated alignment to a seperate file for 
write.fasta(concat_align[1], names = 'ATCC43290', file.out = 'ATCC43290_concat_muscle.fasta')
write.fasta(concat_align[2], names = 'YYC_clinical', file.out = 'YYC_clinical_concat_muscle.fasta')
write.fasta(concat_align[3], names = 'YYC_hottub', file.out = 'YYC_hottub_concat_muscle.fasta')
write.fasta(concat_align[4], names = 'Lens', file.out = 'Lens_concat_muscle.fasta')
write.fasta(concat_align[5], names = 'Philadelphia', file.out = 'philly_concat_muscle.fasta')
write.fasta(concat_align[6], names = 'Thunderbay', file.out = 'TB_concat_muscle.fasta')


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

i<-1
my_sequences<-c(ATCC_df_subset$sequence[i], clin_df_subset$sequence[i], hotb_df_subset$sequence[i], lens_df_subset$sequence[i], phil_df_subset$sequence[i], thun_df_subset$sequence[i])
alignment<-msa::msa(my_sequences, type='dna', order='input')
alignment_convert<-msaConvert(alignment, type="seqinr::alignment") #as alignments are made they can be pasted together!!!!!!! like this : merge<-paste0(test_convert$seq, test_convert$seq)
concat_align<-c('','','','','','')
concat_align<-paste0(concat_align, alignment_convert$seq)