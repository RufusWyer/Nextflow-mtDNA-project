#!/usr/bin/env nextflow

process plinkGPCA {
    input:
    tuple path(maf5geno1_bed), path(maf5geno1_bim), path(maf5geno1_fam)

    output:
    path "plink_gpca.eigenvec", emit: gpca
    tuple path("plink_LD.pruned.bed"), path("plink_LD.pruned.bim"), path("plink_LD.pruned.fam"), emit: pruned

    script:
	"""
	BASENAME=$(basename ${maf5geno1_bed} .bed)

	plink --bfile $BASENAME --indep-pairwise 50 5 0.5 --out plink_LD
	
	plink --bfile $BASENAME --extract plink_LD.prune.in --make-bed --out plink_LD.pruned
	
	plink --bfile plink_LD.pruned --pca 10 header --out plink_gpca
	"""
}
