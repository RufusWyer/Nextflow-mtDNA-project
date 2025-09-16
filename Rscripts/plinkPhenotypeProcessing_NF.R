library(dplyr)

args = commandArgs(trailingOnly = TRUE)

in.rac = args[1] # default: 'outputs/RNA_avcov.csv
in.ids = args[2] # default: 'outputs/plink/gpca_analysis/fullset_LDgpca_50-5-.5.eigenvec'
out.log = args[3] # default: 'outputs/plink/inputs/rac_pheno_log.txt'
out.nolog = args[4] # default: 'outputs/plink/inputs/rac_pheno_nolog.txt'

# Defining function to transform SAMPID to SUBJID format (trim to 9 or 10 chars, sub . for -)
SUBJIDify = function(input){
  input = substr(input, 0, 10)
  input = gsub(".", "-", input, fixed = T)
  input <- ifelse(substr(input, nchar(input), nchar(input)) == "-", 
                  substr(input, 1, nchar(input) - 1), 
                  input)
  return(input)
}

# Reading in RNA average coverage table, and plink IDs (from gpca file, first two cols)
rac_df = read.csv(in.rac, row.names = 1)
plink_ids = read.table(in.ids, header=T)[1:2]

# log transforming the phenotype (for us in heritability estimate)
rac_df$log_mt.RNA.average.coverage = log(rac_df$prop)

# Joining to IDs by SUBJID
plink_ids$SUBJID = plink_ids$FID
rac_pheno = left_join(plink_ids, rac_df, by='SUBJID')

# Separating the log/standard phenotypes and dropping unnecessary columns
raclog_pheno = select(rac_pheno, -c('SUBJID','prop'))
racnolog_pheno = select(rac_pheno, -c('SUBJID','log_mt.RNA.average.coverage'))

# writing final files
write.table(raclog_pheno, out.log, row.names=F, sep="\t", quote = F)
write.table(racnolog_pheno, out.nolog, row.names=F, sep="\t", quote = F)