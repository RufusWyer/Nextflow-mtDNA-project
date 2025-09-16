library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

in.weights = args[1] # Output from the FUSIO.WEIGHTS script
in.pheno = args[2] # this is the fullset pheno, not just train
in.covs = args[3] #'outputs/fusion/fs_covariates.txt'
in.ids = args[4] #'outputs/plink/inputs/fs_test_subset.txt'
out.res = args[5]

# reading in all inputs
pheno = read.table(in.pheno, header=T)
lasso = read.table(in.weights, header=T)
covar = read.table(in.covs, header=T)
test_ids = read.table(in.ids, header=F)

# Removing NAs from covar file
covar = na.omit(covar)

# filtering pheno file to match the covars (removing NAs)
pheno = pheno[ which(pheno$IID %in% covar$IID), ]

# regressing out the covariates from the phenotype and rescaling, allowing it to be used in comparison with the lasso sumscore results
reg = summary(lm( pheno[,3] ~ as.matrix(covar[,3:ncol(covar)]) ))
new_pheno = data.frame(IID = pheno$IID, pheno = scale(reg$resid))

# Splitting new phenotype into train and test sets
train_pheno = new_pheno[which(!(new_pheno$IID %in% test_ids$V1)), ]
test_pheno = new_pheno[which(new_pheno$IID %in% test_ids$V1), ]

# Joining the lasso results and calculated new test pheno into a single table
lasso = left_join(test_pheno, lasso, by='IID')

# Compute mean/SD of genetic scores in the TRAINING set
test_mean <- mean(test_pheno$pheno)
test_sd <- sd(test_pheno$pheno)

# Standardize test set scores
lasso$SCORE_std = (lasso$SCORESUM - test_mean) / test_sd

# Extracting and saving stats for final results table (published)
lasso_model = lm(lasso$pheno ~ lasso$SCORE_std)
lasso_sum = summary(lasso$SCORE_std)
lasso_rmse <- sqrt(mean((lasso$pheno - lasso$SCORE_std)^2))
lasso_cor = cor.test(lasso$pheno, lasso$SCORE_std)
lasso_rsq = (lasso_cor$estimate)^2
lasso_pval = lasso_cor$p.value
res_rows = c('R^2','RMSE','pval')
lasso_res = c(lasso_rsq, lasso_rmse, lasso_pval)
res_df = data.frame(LASSO = lasso_res, row.names = res_rows)

# Writing final results table
write.table(res_df, quote=F, row.names=T, col.names=T,file=out.res)
