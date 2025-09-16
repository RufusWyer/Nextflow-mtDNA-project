#!/usr/bin/env nextflow

process plinkAssoc {
    input:
    tuple path(maf5geno1_bed), path(maf5geno1_bim), path(maf5geno1_fam), path(plinkMasterCovariates), path(plink_phenotype), path(plink_covlist)

    output:
    path "assoc*.glm.linear", emit: assoc
    tuple path("assoc*.glm.linear.adjusted"), path( "assoc*.glm.linear") emit: publish

    script:
	"""
	BASENAME=$(basename ${maf5geno1_bed} .bed)
	COVLIST=$(head -n 1 ${plink_covlist})

	plink2 --bfile $BASENAME \
        --linear hide-covar --covar-variance-standardize --ci 0.95 \
        --pheno ${plink_phenotype} \
        --covar-name  --adjust \
        --covar $COVLIST \
	--out assoc

	"""
}

