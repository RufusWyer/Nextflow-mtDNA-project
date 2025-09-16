library(dplyr)

args = commandArgs(trailingOnly = TRUE)


in.pmc = args[1] # default: 'fullset_pmc_plinkMasterCovariates_fstrain.txt'
in.covlist = args[2] # default: 'covlistVIF10.txt'
out.covs = args[3] # default: 'outputs/fusion/fstrain_covariates.txt'
out.quant = args[4] # default: 'outputs/fusion/fstrain_covariates_quant.txt'
out.disc = args[5] # default: 'outputs/fusion/fstrain_covariates_disc.txt'

# reading in pmc covariates, and VIF trimmed covlist
pmc = read.table(in.pmc, header=T)
covlist = read.table(in.covlist, header=T, sep = ',')

# select for only the IDs, and the names of the covars in the covlist file
covs = select(pmc, c('FID','IID', colnames(covlist)))

# Writing the full covars file
write.table(covs, out.covs, row.names=F, sep="\t", quote = F)

# Removing SEX from quantitative covars, and placing them in discrete covars
# Heritability estimate requires the discrete and quantitative covars to be separated
covs_quant = select(covs, -c('SEX'))
covs_disc = select(covs, c('FID','IID','SEX'))

# writing final split files
write.table(covs_quant, out.quant, row.names=F, sep="\t", quote = F)
write.table(covs_disc, out.disc, row.names=F, sep="\t", quote = F)
