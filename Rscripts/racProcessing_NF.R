library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

in.counts = args[1] # default: 'originals/gene_reads_whole_blood_v8.gct'
in.lengths = args[2] # default: 'originals/gene_lengths.txt'
out.rac = args[3] # default: 'outputs/RNA_avcov.csv'

# Defining function to transform SAMPID to SUBJID format (trim to 9 or 10 chars, sub . for -)
SUBJIDify = function(input){
  input = substr(input, 0, 10)
  input = gsub(".", "-", input, fixed = T)
  input <- ifelse(substr(input, nchar(input), nchar(input)) == "-", 
                  substr(input, 1, nchar(input) - 1), 
                  input)
  return(input)
}

# Reading in counts data and gene length data
wbcounts = read.table(in.counts,
                      skip=2, header = T)
gene_lengths = read.table(in.lengths, header = T)

# Separating counts into mitochondrial and nuclear (wbcounts), by gene description
mtnames_all = wbcounts$Description[grep('^MT-', wbcounts$Description)]
mtcounts = wbcounts[ which(wbcounts$Description %in% mtnames_all), ]
wbcounts = wbcounts[ which(!(wbcounts$Description %in% mtnames_all)), ]

# Formatting, removal of id and Description (gene name) from dfs
wbcounts = select(wbcounts, -c('id','Description'))
mtcounts = select(mtcounts, -c('id','Description'))

# Reformatting gene-lengths in prep for left-join, gene column must be called Name and only
# Name and LENGTH need to be retained
gene_lengths$Name = gene_lengths$GENE
gene_lengths = select(gene_lengths, -c('LENGTH_KB', 'GENE'))

# Use left join to attach Length to data, then immediately extract length to new variable
# without saving joined df. This is necessary to make sure that the IDs match, as they are
# not in the same order.
# Length is not saved in the dfs as it will make operations over the df harder later
wblength = (left_join(wbcounts, gene_lengths, by='Name'))$LENGTH
mtlength = (left_join(mtcounts, gene_lengths, by='Name'))$LENGTH

# Reformatting counts, by moving gene name from a var to rownames, and transposing so
# genes are now columns and samples are rows
rownames(wbcounts) = wbcounts$Name
rownames(mtcounts) = mtcounts$Name
wbcounts = as.data.frame(t(select(wbcounts, -c('Name'))))
mtcounts = as.data.frame(t(select(mtcounts, -c('Name'))))

# calculating average coverage for nuclear and mt genomes for each sample
# (total reads/total region length)
wbcounts$nuc_cov = rowMeans(wbcounts/wblength)
mtcounts$mt_cov = rowMeans(mtcounts/mtlength)

# reformatting in prep for left join, only need the coverage and SUBJID for this
wbcounts$SUBJID = SUBJIDify(rownames(wbcounts))
wbcounts_cov = select(wbcounts, c('SUBJID', 'nuc_cov'))
mtcounts$SUBJID = SUBJIDify(rownames(mtcounts))
mtcounts_cov = select(mtcounts, c('SUBJID', 'mt_cov'))

# left joining and calculating RNA average coverage proportion phenotype
# essentially an RNA analogue of mtDNA copy number
cov_df = left_join(wbcounts_cov, mtcounts_cov, by='SUBJID')
cov_df$prop = (cov_df$mt_cov/cov_df$nuc_cov)*2

# selecting only relevant columns
avcov_df = select(cov_df, c('SUBJID','prop'))

# writing RNA average coverage csv
write.csv(avcov_df, out.rac)
