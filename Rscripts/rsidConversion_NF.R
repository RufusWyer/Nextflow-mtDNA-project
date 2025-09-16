library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)

in.bim = args[1] # '/project/outputs/plink/clumping/fstrain.1e-1.filtered.bim'
in.wgt = args[2] # '/project/outputs/fusion/fusion.weights.fstrain.1e-1.log.wgt.RDat'
out.snplist = args[3] # '/project/outputs/plink/clumping/fstrain.1e-1.snplist.rsid'
out.snpmap= args[4] # '/project/outputs/plink/clumping/fstrain.1e-1.snpmap.rsid'
out.wgt = args[5] # '/project/outputs/fusion/fusion.weights.fstrain.1e-1.log.rsid.wgt.RDat'

# reading in and reformatting the plink SNPs from bim file
plink_bim = read.table(in.bim)
colnames(plink_bim)[2] <- 'snp_id'
plink_snps = plink_bim

# Plink ID system is made up of chr, pos, ref, and alt in succesion. This parses them to separate these into individual values
# Parse PLINK IDs into chr, pos, ref, alt
plink_snps$chr <- sub("^chr([0-9XYM]+)_.*", "\\1", plink_snps$snp_id)
plink_snps$pos <- as.numeric(sub("^chr\\w+_(\\d+)_.*", "\\1", plink_snps$snp_id))
plink_snps$ref.x <- sub("^chr\\w+_\\d+_([ACGT]+)_.*", "\\1", plink_snps$snp_id)
plink_snps$alt.x <- sub("^chr\\w+_\\d+_[ACGT]+_([ACGT]+)_.*", "\\1", plink_snps$snp_id)

plink_snps$chr[plink_snps$chr == "26"] = 'MT'

# To convert into rsid formats, the positions and allels of the SNPs to be cross-referenced with a known database
# Create a GRanges object
snp_ranges <- GRanges(
  seqnames = plink_snps$chr,
  ranges = IRanges(plink_snps$pos, plink_snps$pos)
)
# Load dbSNP database
snp_db <- SNPlocs.Hsapiens.dbSNP155.GRCh38
# Find overlapping dbSNP records
matched <- snpsByOverlaps(snp_db, snp_ranges)
# Convert matched SNPs to data.frame
matched_df <- as.data.frame(matched)
colnames(matched_df)[colnames(matched_df) == "start"] <- "pos"  # rename for merging
colnames(matched_df)[colnames(matched_df) == "seqnames"] <- "chr"
# Merge PLINK SNPs with dbSNP
merged <- merge(plink_snps, matched_df, by = c("chr", "pos"))

# Double checking matches with ref alleles against IUPAC ambiguity codes
iupac_ambiguity_map <- list(
  A = "A",
  C = "C",
  G = "G",
  T = "T",
  R = c("A", "G"),
  Y = c("C", "T"),
  S = c("G", "C"),
  W = c("A", "T"),
  K = c("G", "T"),
  M = c("A", "C"),
  B = c("C", "G", "T"), # B = C/G/T (not A)
  D = c("A", "G", "T"), # D = A/G/T (not C)
  H = c("A", "C", "T"), # H = A/C/T (not G)
  V = c("A", "C", "G"), # V = A/C/G (not T)
  N = c("A", "C", "G", "T") # N = A/C/G/T (any base)
)

merged$ref_match <- mapply(function(ambig_code, ref) {
  possible_alleles <- iupac_ambiguity_map[[ambig_code]]
  return(ref %in% possible_alleles)
}, merged$alleles_as_ambig, merged$ref.x)


# Keep only good matches, and reformat
final_snps <- merged[merged$ref_match == TRUE, ]
final_snps = data.frame(snp_id = final_snps$snp_id, rsid = final_snps$RefSNP_id)

# reformat plink_bim and SNPlist to replace old IDs with new
plink_bim = merge(plink_bim, final_snps, by = 'snp_id')
snplist = as.data.frame(plink_bim$snp_id)
plink_bim$snp_id = plink_bim$rsid
plink_bim = subset(plink_bim, select = -c(rsid))
plink_bim = plink_bim[, c("V1", setdiff(names(plink_bim), "V1"))]

# Loading weights file, and subbing in rsid SNPs for old SNPs
load(file = in.wgt)
wgt.matrix = as.data.frame(wgt.matrix)
wgt.matrix$snp_id = row.names(wgt.matrix)
wgt.matrix = merge(wgt.matrix, final_snps, by = "snp_id")
rownames(wgt.matrix) = wgt.matrix$rsid
wgt.matrix = subset(wgt.matrix, select = -c(snp_id, rsid))
wgt.matrix = as.matrix(wgt.matrix)

snps$snp_id = snps$V2
snps = merge(snps, final_snps, by = "snp_id")
snps$V2 = snps$rsid
snps = subset(snps, select = -c(snp_id, rsid))

final_snps <- final_snps[!duplicated(final_snps$snp_id), ]

# Writing final files
write.table(snplist, file = out.snplist, quote = F, sep = "\t", row.names = F, col.names = F)
write.table(final_snps, file = out.snpmap, quote = F, sep = "\t", row.names = F, col.names = F)
save(wgt.matrix , snps , cv.performance , hsq, hsq.pv, N.tot, 
     file = out.wgt)
