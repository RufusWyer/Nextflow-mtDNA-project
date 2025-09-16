library(dplyr)

args = commandArgs(trailingOnly = TRUE)

in.pmc = args[1] # default: 'plinkMasterCovariates.txt'
out.train = args[2] # default: 'train_subset.txt'
out.test = args[3] # default: 'test_subset.txt'

# setting seed to ensure consistent test-train split
set.seed(123)

# pmc (covars) read in for IDs (FID, IID)
plink_ids = read.table(in.pmc, header=T)

# ISC. TIME extracted as some are NA. Therefore omitting all IDs with NA ISC.TIME
plink_ids = data.frame(FID = plink_ids$FID, IID = plink_ids$IID, ISC = plink_ids$SMTSISCH)
plink_ids = na.omit(plink_ids)

# Removing ISC.TIME after NAs are filtered
plink_ids = data.frame(FID = plink_ids$FID, IID = plink_ids$IID)

# splitting IDs into train 80:20 train
plink_ids$split = 1:nrow(plink_ids)
train = plink_ids %>% sample_frac(0.80)
test  = anti_join(plink_ids, train, by = 'split')

# Removing the split column, leaving only ID columns
train = train[1:2]
test = test[1:2]

# Writing final ID files
write.table(out.train, , row.names=F, sep="\t", quote = F, col.names = F)
write.table(out.test, , row.names=F, sep="\t", quote = F, col.names = F)