# This is only to be run using a docker image - See peer_docker_readme.Rmd for instructions
# This is currently set up to be run on the HPC

library(peer)

args = commandArgs(trailingOnly = TRUE)

print(args)

# default: 'originals/gene_tpm_whole_blood_v8.gct'
in.tpm = args[1]

# default: 'outputs/RNA_avcov.csv'
in.rac = args[2]

# default: 'outputs/gpca.eigenvec'
in.ids = args[3]

# default: '/project/outputs/mtcPF.factors.csv'
out.pf = args[4]

SUBJIDify = function(input){
  input = substr(input, 0, 10)
  input = gsub(".", "-", input, fixed = T)
  input <- ifelse(substr(input, nchar(input), nchar(input)) == "-", 
                  substr(input, 1, nchar(input) - 1), 
                  input)
  return(input)
}

keep_ids = read.table(in.ids)[1]
colnames(keep_ids) = c('SUBJID')

wbtpm = read.table(in.tpm, skip=2, header = T)
row.names(wbtpm) = wbtpm$Name
wbtpm = wbtpm[4:length(wbtpm)]
wbtpm = as.data.frame(t(wbtpm))
rownames(wbtpm) = SUBJIDify(rownames(wbtpm))
wbtpm = wbtpm[ which(SUBJIDify(rownames(wbtpm)) %in% keep_ids$SUBJID), ]

head(rownames(wbtpm))

rpdf = read.csv(in.rac, header= T, row.names = 2)[2]
head(rpdf)

common_rows = intersect(rownames(rpdf), rownames(wbtpm))
rpdf = rpdf[common_rows, , drop = FALSE]
print('identical:')
identical(rownames(rpdf), rownames(wbtpm))
identical(sort(rownames(rpdf)), sort(rownames(wbtpm)))


model = PEER()
PEER_setPhenoMean(model,as.matrix(wbtpm))
PEER_setNk(model,10)
PEER_setNmax_iterations(model, 5000)
PEER_setCovariates(model, as.matrix(rpdf$prop))
PEER_update(model)

factors = as.data.frame(PEER_getX(model))
rownames(factors) = rownames(wbtpm)
write.csv(factors, out.pf)
