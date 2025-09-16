# to be run on the HPC only
# following this, run plink.sh

library(dplyr)
library(readr)

args = commandArgs(trailingOnly = TRUE)

print(args)

# default: 'outputs/peer_analysis/wb_covs_celltypes.csv'
in.covs = args[1]

# default: 'originals/GTEx_v8.maf1hwe.ALLinds.MT.fam'
in.fam = args[2]

# default: 'outputs/plink/inputs/fullset_ids_subset.txt'
out.ids = args[3]


# loading in covariate data to extract subject ids
covs = read.csv(in.covs, row.names = 1)

#keep_ids = as.data.frame(read.csv('outputs/keep_ids.csv', header = T))


print('covs read')
ids = as.factor(covs$SUBJID)
print(length(ids))
print('ids set')

fam = read.table(in.fam, header=F)
print(dim(fam))

# creating subset text file for plink filtering
# V1 and V2 are the family and within-family ids, both identical and equivalent to subjid
fam %>%
  select(V1, V2) %>%
  filter(V1 %in% ids) %>%
  write_delim(out.ids, col_names = F)
print('subset txt written')