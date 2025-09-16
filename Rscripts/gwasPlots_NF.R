library(qqman)

args <- commandArgs(trailingOnly = TRUE)

in.glm = args[1] # glm.linear file
out.man = args[2] # 'mtDNA.manPlot.png'
out.qq = args[3] # 'mtDNA.qqPlot.png'

# reading in the glm.linear file of the gwas results
df = read.table(in.glm, header=T, comment.char = "")

# Renaming columns to match with qqman settings
colnames(df)[1:3] = c('CHR', 'BP', 'SNP')

# Removing X,Y,XY, and MT chromosome data. Only testing non-sex nuclear chromosomes
df <- df[!df$CHR %in% c("X", "Y", "XY", "MT"), ]

# Cleaning CHR column for non-numeric and NAs
df$CHR = as.numeric(df$CHR)
df = na.omit(df)

# Writing QQ plot
png(out.qq, width = 800, height = 800)
qq(df$P, main = 'Q-Q plot of mitochondrial expression GWAS (RNA average coverage phenotype)')
dev.off()
  
# Writing Manhattan plot, with suggestive line off
png(out.man, width = 800, height = 600)
manhattan(df, main = 'Manhattan plot of mitochondrial expression GWAS (RNA average coverage phenotype)', suggestiveline = F)
dev.off()


