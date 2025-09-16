#!/usr/bin/env nextflow

process sumstatsProcessing {
    output:
    path("*non.processesd.txt"), emit: nonwb
    path("*wb.processesd.txt"), emit: wb

    script:
	"""
// Downloading GWAS sumstats for required diseases, unzipping, and trimming unnecessary columns to prevent errors
// Non-blood diseases
// type 2 diabetes
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479885/GCST90479885.tsv.gz
gunzip GCST90479885.tsv.gz
cut -f1-9 GCST90479885.tsv > type2Diabetes.non.trimmed.tsv

// Coronary Artery Disease
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479543/GCST90479543.tsv.gz
gunzip GCST90479543.tsv.gz
cut -f1-9 GCST90479543.tsv > coronaryArtery.non.trimmed.tsv

// Chronic kidney disease
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90480001-GCST90481000/GCST90480376/GCST90480376.tsv.gz
gunzip GCST90480376.tsv.gz
cut -f1-9 GCST90480376.tsv > chronicKidney.non.trimmed.tsv

// Parkinson's disease
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90480001-GCST90481000/GCST90480008/GCST90480008.tsv.gz
gunzip GCST90480008.tsv.gz
cut -f1-9 GCST90480008.tsv > parkinsons.non.trimmed.tsv

// Multiple Sclerosis
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90480001-GCST90481000/GCST90480014/GCST90480014.tsv.gz
gunzip GCST90480014.tsv.gz
cut -f1-9 GCST90480014.tsv > multipleSclerosis.non.trimmed.tsv

// Alzheimer's disease
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007320/GCST007320_GRCh37.tsv.gz
gunzip GCST007320_GRCh37.tsv.gz

// Long COVID
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90454001-GCST90455000/GCST90454541/GCST90454541.tsv

// Blood diseases
// Thrombocytopenia
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479981/GCST90479981.tsv.gz
gunzip GCST90479981.tsv.gz
cut -f1-9 GCST90479981.tsv > thrombocytopenia.wb.trimmed.tsv

// Anemia
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479972/GCST90479972.tsv.gz
gunzip GCST90479972.tsv.gz
cut -f1-9  GCST90479972.tsv > anemia.wb.trimmed.tsv

// Eosinophilia
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479988/GCST90479988.tsv.gz
gunzip GCST90479988.tsv.gz
cut -f1-9 GCST90479988.tsv > eosinophilia.wb.trimmed.tsv

// Myeloid Leukemia
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479832/GCST90479832.tsv.gz
gunzip GCST90479832.tsv.gz
cut -f1-9 GCST90479832.tsv > mLeukemia.wb.trimmed.tsv

// Lymphoid Leukemia
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479831/GCST90479831.tsv.gz
gunzip GCST90479831.tsv.gz
cut -f1-9 GCST90479831.tsv > lLeukemia.wb.trimmed.tsv

// Non-Hodgkins Lymphoma
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479828/GCST90479828.tsv.gz
gunzip GCST90479828.tsv.gz
cut -f1-9 GCST90479828.tsv > nLymphoma.wb.trimmed.tsv

// Hodgkins Lymphoma
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479827/GCST90479827.tsv.gz
gunzip GCST90479827.tsv.gz
cut -f1-9 GCST90479827.tsv > hLymphoma.wb.trimmed.tsv

Rscript Rscripts/sumstatsProcessing_NF.R ".*\\.trimmed\\.tsv$" "GCST90454541.tsv" "GCST007320_GRCh37.tsv.gz"

	"""
}
