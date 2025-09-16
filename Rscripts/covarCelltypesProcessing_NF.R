
args <- commandArgs(trailingOnly = TRUE)

in.tpm = args[1] # RNAseq tpm file
in.samples = args[2] # metadata sample attributes file
in.phs = args[3] # metadat subject phenotypes
out.wbcovs_ct = args[4] # output: 'wb_covs_celltypes.csv'

# Installing packages
install.packages('devtools')
devtools::install_github('AlmogAngel/xCell2')
library(xCell2)
print('xCell2 LOADED')

library(dplyr)
print('dplyr LOADED')

# First, must download the pre-trained reference
# Immune compendium reference dataset chosen, immune/blood cells, homo sapiens
# See more info at section 3 of this doc: https://aran-lab.com/xcell2-vignette#introduction-to-cell-type-enrichment-analysis-with-xcell-2.0
# Set the URL of the pre-trained reference
ref_url <- "https://dviraran.github.io/xCell2refs/references/ImmuneCompendium.xCell2Ref.rds"
# Set the local filename to save the reference
local_filename <- "ImmuneCompendium.xCell2Ref.rds"
# Download thes reference file
download.file(ref_url, local_filename, mode = "wb")
# Load the downloaded reference
ImmuneCompendium.xCell2Ref <- readRDS(local_filename)

# Reading GTEx whole blood tpm data, skip=2 necessary for GCT files
wbtpm = read.table(in.tpm, skip=2, header = T)

# There are multiple rows with duplicate descriptions (gene names), with different IDs
# xCell2 requires gene names as rownames in input matrix, and rownames cannot be duplicated
# This sums the TPM data for each duplicate gene name row
wbtpm_ct = aggregate(wbtpm[4:length(wbtpm)], by=wbtpm['Description'], sum)
# Setting rownames to gene name (description), then removing description column
row.names(wbtpm_ct) = wbtpm_ct$Description
wbtpm_ct = wbtpm_ct[-1]

# Running xCell2 analysis using whole blood TPM matrix, and chosen reference
celltypes <- xCell2::xCell2Analysis(mix = wbtpm_ct,
                                         xcell2object = ImmuneCompendium.xCell2Ref)

# Now celltypes are calculated, need to extract the relevant metadata info for other covariates
# reformat tpm file to remove unnecessary columns
wbtpm = t(wbtpm[4:length(wbtpm)])

# reformat sample names to match the SUBJID (first 10 chars only, subbed . for -)
SUBJID = substr(rownames(wbtpm), 0, 10)
SUBJID = gsub(".", "-", SUBJID, fixed = T)

# some samples are only 9 characters (with - as 10th char), so need adjusting as well
SUBJID <- ifelse(substr(SUBJID, nchar(SUBJID), nchar(SUBJID)) == "-", 
                 substr(SUBJID, 1, nchar(SUBJID) - 1), 
                 SUBJID)
print('SUBJID CREATED')

# reformat sample names to match the SAMPID (sub . for -)
SAMPID = gsub(".", "-", rownames(wbtpm), fixed = T)
print('SAMPID CREATED')

# create df with these new ID variables
wbcovs = data.frame(SUBJID, SAMPID, row.names = rownames(wbtpm))

# read in GTEx metadata files
samples = read.table(in.samples, sep="\t", header = T,
                     quote = "", fill = TRUE, comment.char = "")
print('samples READ')

phs = read.table(in.phs, sep="\t", header = T,
                     quote = "", fill = TRUE, comment.char = "", skip=10)
print('phs READ')

# left-join tables, keeping only the required columns for covariate analysis
# Extracting RIN (SMRIN), Sample ID (SAMPID), and Sample Ischemic Time (SMTSISCH)
# one-to-one relationships specified
samples = data.frame(SAMPID = samples$SAMPID,SMRIN = samples$SMRIN,
                     SMTSISCH = samples$SMTSISCH)
wbcovs = left_join(wbcovs, samples, by='SAMPID', relationship = 'one-to-one')
print('samples JOINED')

# Extracting Age, Subject ID (SUBJID), and Sex
phs = data.frame(SUBJID = phs$SUBJID, AGE = phs$AGE, SEX = phs$SEX)
wbcovs = left_join(wbcovs, phs, by='SUBJID', relationship = 'one-to-one')
print('phs JOINED')

# Incorporating celltype data from earlier into covariates
# Transposing, and creating SAMPID from rownames for join
celltypes = as.data.frame(t(celltypes))
celltypes$SAMPID = rownames(celltypes)
celltypes$SAMPID = gsub(".", "-", celltypes$SAMPID, fixed = T)
print('SAMPID CT CREATED')

wbcovs_ct = left_join(celltypes, wbcovs, by='SAMPID', relationship = 'one-to-one')
print('CT JOINED')

# Checking for NAs
print(colSums(is.na(wbcovs_ct)))

# writing final file, with CTs and metadata covars
write.csv(wbcovs_ct, out.wbcovs_ct)