library(dplyr)

args = commandArgs(trailingOnly = TRUE)

in.covs = args[1] # default: 'outputs/peer_analysis/wb_covs_celltypes.csv'
in.pfs = args[2] # default: 'outputs/peer_analysis/peer_results/FSPF.factors.csv'
in.gpca = args[3] # default: 'outputs/plink/gpca_analysis/fullset_LDgpca_50-5-.5.eigenvec'
out.pmc = args[4] # default: 'outputs/plink/inputs/plinkMasterCovariates.txt'

# Defining function to transform SAMPID to SUBJID format (trim to 9 or 10 chars, sub . for -)
SUBJIDify = function(input){
  input = substr(input, 0, 10)
  input = gsub(".", "-", input, fixed = T)
  input <- ifelse(substr(input, nchar(input), nchar(input)) == "-", 
                  substr(input, 1, nchar(input) - 1), 
                  input)
  return(input)
}

# Reading in covs/ct, peerfactors, gpca
# For PFs, removing first column (mean covariate), and second column for PF (mt covariate)
wbcovs = read.csv(in.covs, row.names = 1)
gpca = read.table(in.gpca, header=T)
PF = read.csv(in.pfs, row.names = 1)[2:11]

# formatting PFs, renaming columns and setting SUBJID
colnames(PF) = paste('PF', seq(1,10), sep='')
PF$SUBJID = SUBJIDify(rownames(PF))

# Replace 1:40 cols with new CT ids, keep 41:end the same
# Formatting, easier to read and to avoid formatting errors later down the line
ids = paste('CT', seq(1,40), sep='')
colnames(wbcovs) = c(ids, colnames(wbcovs[41:length(wbcovs)]))

# joining all 3 tables together by SUBJID
master_covs = left_join(PF, wbcovs, by='SUBJID')
gpca$SUBJID = gpca$IID
master_covs = left_join(gpca, master_covs, by='SUBJID')

# Removing unnecessary ID columns
master_covs = select(master_covs, -c('SUBJID', 'SAMPID'))

# Writing final file
write.table(master_covs,out.pmc, row.names=F, sep="\t", quote = F)