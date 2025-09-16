args <- commandArgs(trailingOnly = TRUE)

in.trimmed = args[1] # ".*\\.trimmed\\.tsv$"
in.covid = args[2] # "GCST90454541.tsv"
in.alz = args[3] # "GCST007320_GRCh37.tsv.gz"

# 12 of 14 files can all be processed in the same way, as they are from the same study
for (filename in list.files(in.trimmed,full.names = T)){
  df = read.table(filename,
                     header = TRUE,
                     sep = "",             
                     fill = TRUE,          
                     comment.char = "",    
                     quote = "")
  # Generate Z value from pvalue and odds ratio
  df$odds_ratio = as.numeric(df$odds_ratio)
  df$p_value = as.numeric(df$p_value)
  df = na.omit(df)
  df$Z = qnorm(df$p_value/ 2, lower.tail = FALSE)
  df$Z = ifelse(log(df$odds_ratio) < 0, -df$Z, df$Z)
  
  # Rename ref, alt, and ID to A1, A2, SNP, and Z to match with FUSION's reqs
  colnames(df)[3:4] = c('A1','A2')
  colnames(df)[9] = 'SNP'
  df = subset(df, select = c('A1','A2','SNP','Z'))
  
  # Formatting filename and saving
  save_name = basename(filename)
  save_name = substr(save_name, 1, nchar(save_name)-12)
  save_name = paste0(save_name, '.processed.txt', sep='')
  write.table(df, file = save_name,
            quote = F, sep = "\t", row.names = F, col.names = T)
}


# Long COVID19 sumstats are formatted differently, so require different processing
df = read.table(in.covid,
                  header = TRUE,
                  sep = "",             # allow flexible whitespace
                  fill = TRUE,          # fill missing fields
                  comment.char = "",    # do not treat # as comment
                  quote = "")

# Calculating Z value from beta and SE instead of odds ratio
df$Z = df$beta/df$standard_error
df = na.omit(df)

# Rename ref, alt, and ID to A1, A2, SNP, and Z to match with FUSION's reqs
colnames(df)[3:4] = c('A1','A2')
colnames(df)[10] = 'SNP'
df = subset(df, select = c('A1','A2','SNP','Z'))

# Saving final file
write.table(df, file = 'longCOVID.non.processed.txt',
            quote = F, sep = "\t", row.names = F, col.names = T)

# Alzheimers sumstats are formatted differently, so require different processing
df = read.table(in.alz,
                header = TRUE,
                sep = "",             # allow flexible whitespace
                fill = TRUE,          # fill missing fields
                comment.char = "",    # do not treat # as comment
                quote = "")

# Calculating Z value from beta and SE instead of odds ratio
df$Z = df$beta/df$standard_error
df = na.omit(df)
# Rename ref, alt, and ID to A1, A2, SNP, and Z to match with FUSION's reqs
colnames(df)[3:5] = c('A1','A2','SNP')
df = subset(df, select = c('A1','A2','SNP','Z'))

# Saving final file
write.table(df, file = 'alzheimers.non.processed.txt',
            quote = F, sep = "\t", row.names = F, col.names = T)




