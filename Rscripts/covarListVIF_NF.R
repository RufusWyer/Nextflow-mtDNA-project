library(rvif)
library(dplyr)


in.covs = args[1] # default: 'outputs/peer_analysis/wb_covs_celltypes.csv'
in.pfs = args[2] # default: 'outputs/peer_analysis/peer_results/FSPF.factors.csv'
in.gpca = args[3] # default: 'outputs/plink/gpca_analysis/fullset_LDgpca_50-5-.5.eigenvec'
out.vif = args[4] # default: 'outputs/plink/inputs/CTmtcPFASRI_v10_list.txt'

# Defining function to transform SAMPID to SUBJID format (trim to 9 or 10 chars, sub . for -)
SUBJIDify = function(input){
  input = substr(input, 0, 10)
  input = gsub(".", "-", input, fixed = T)
  input <- ifelse(substr(input, nchar(input), nchar(input)) == "-", 
                  substr(input, 1, nchar(input) - 1), 
                  input)
  return(input)
}

# Reading in all covar files
mtcPF = read.csv(in.pfs, row.names = 1)[2:11]
wbcovs = read.csv(in.covs, row.names = 1)
gpca = read.table(in.gpca, header=T)

# Reformatting PFs and CTs
colnames(mtcPF) = paste('PF', seq(1,10), sep='')
mtcPF$SUBJID = SUBJIDify(rownames(mtcPF))
ids = paste('CT', seq(1,40), sep='')
colnames(wbcovs) = c(ids, colnames(wbcovs[41:length(wbcovs)]))

# Joining tables and removing IDs that are not present in both RNA and WGS
wbcovs = left_join(wbcovs, mtcPF, by='SUBJID')
wbcovs = wbcovs[ which((wbcovs$SUBJID) %in% gpca$FID), ]
gpca = gpca[ which((gpca$FID) %in% wbcovs$SUBJID), ]
gpca$SUBJID = gpca$FID
wbcovs = left_join(wbcovs, gpca, by='SUBJID')

# Removing unnecessary columns
wbcovs = select(wbcovs, -c('SAMPID','SUBJID', 'FID', 'IID'))


vif_df = na.omit(wbcovs)
while (TRUE){
    # calculate VIF values, no intercept (all indepedent variables, not modelled)
    vif = RVIF(vif_df, intercept = F)
    # vif loses var names as default
    rownames(vif) = colnames(vif_df)
    if (max(vif$RVIF) < 10){
      break
    }
    # extract highest VIF value and save values
    max_name <- rownames(vif)[which.max(vif$RVIF)]
    # If it's 'PF1', exclude it and get the next highest
    if (max_name == 'PF1') {
      vif_excl <- vif[!(rownames(vif) %in% 'PF1'), ]
      max_name <- rownames(vif_excl)[which.max(vif_excl$RVIF)]
    }
    # remove highest VIF val variable
    vif_df = select(vif_df, -(all_of(max_name)))
}

# Writing file with list of all cvoars still remaining at VIF10
cov_string = paste(colnames(vif_df), collapse=',')
write(cov_string, out.vif)



