# Nextflow-mtDNA-project
This project is an investigation into the genetic regulation of the mitochondrial transcriptome in whole blood, and its role in common disease risk. This is the culmination of several months of work under project supervisor Alan Hodgkinson, as an extended research project during my MSc in Applied Bioinformatics at King's College London. 

The entire workflow is contained within a single Nextflow pipeline to maximise organisation, while being an exercise to improve my skills in this language. Unfortunately this project requires sensitive data that cannot be shared on Github, however the information for all requirements is detailed below.

A brief summary of the workflow:
Firstly, we explored the relationship between nuclear genetic variation and overall mitochondrial expression using GWAS modelling with PLINK, provising associated nuclear variants  could be candidates for regulating mitochondrial expression. Secondly, we built a machine learning model to impute overall mitochondrial expression in whole blood from variation in the nuclear genome, using FUSION. This is achieved by leveraging scores from our association model to filter out noise from variants below a certain association threshold, allowing us to predict whole blood mitochondrial expression in new individuals solely through their genetic variation. Through this, we provided proof of concept for applying expression imputation models to mtDNA and produce quantifiable evidence of the regulatory function of the nuclear genome on mitochondrial expression. Thirdly, we tested for associations between common diseases and our imputed nuclear-controlled mitochondrial expression, isolating the potential causal role of mtDNA in this. We achieved this by conducting a series of transcriptome wide association studies (TWAS) between our expression imputation model and summaries of previous genome wide association studies (GWAS) on common diseases, again using FUSION.

The workflow outputs a full set of results for the initial GWAS (QQ and Manhattan plots, and full GWAS scores by SNP), an evaluation of the FUSION model prediction scores on new individuals (including R-squared and Root Mean Squared Error), and association results from TWAS between the FUSION model and 14 common diseases (7 blood-based, 7 other).

REQUIRED SOFTWARE
Added to path:
Plink v1.9 (https://www.cog-genomics.org/plink/)
Plink2 v2.00 (https://www.cog-genomics.org/plink/2.0/)
Singularity

Downloaded into /software:
FUSION (2022/02/01 release - http://gusevlab.org/projects/fusion/)
GCTA v1.95.0 (https://yanglab.westlake.edu.cn/software/gcta/#Download)

R packages:
dplyr v1.1.4
xCell2 v2.0 (installed from github AlmogAngel/xCell2)
r-peer (via singularity, run singularity/setup.sh)
readr v2.1.5
rvif v3.1
qqman v0.1.9
SNPlocs.Hsapiens.dbSNP155.GRCh38 v	0.99.24
GenomicRanges v1.60.0

REQUIRED INPUTS (detailed in nextflow.config)
The following are downloadable from GTEx: (https://gtexportal.org/home/)
tpm = "path/to/gene_tpm_whole_blood_v8.gct"
samples = "path/to/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
phs = "path/to/phs000424.v9.pht002742.v9.p2.c1.GTEx_Subject_Phenotypes.GRU.txt"
counts = "path/to/gene_reads_whole_blood_v8.gct"

The following WGS data underwent additional pre-processing prior to the start of this project, as detailed here. WGS data were acquired as VCF (nuclear genome only) for 838 individuals via dbGaP (phs000424.v8.p2). The nuclear VCF data then underwent additional processing using Plink (v1.90), filtering for variant and genotype Phred score (>40), minor allele frequency (>5%),
genotype missingness rate (<1%), and adherence to Hardy-Weinberg equilibrium (P>0.001), outputting to a bfile trio.
wgsfam = "path/to/GTEx_v8.maf1hwe.ALLinds.MT.fam"
wgsbim = "path/to/GTEx_v8.maf1hwe.ALLinds.MT.bim"
wgsbed = "path/to/GTEx_v8.maf1hwe.ALLinds.MT.bed"
