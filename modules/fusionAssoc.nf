#!/usr/bin/env nextflow

process fusionAssoc {
    publishDir "results/fusionAssoc", mode: 'copy'

    input:
    path sumstat
    path weights
    tuple path(ref_bed), path(ref_bim), path(ref_fam)

    output:
    path("*.results.dat")

    script:
	"""
	echo "WGT     ID      CHR     P0      P1" > wgtlist
        echo "${weights}   mt_rac  1       1       1" >> wgtlist


        SS_BASENAME=$(basename "${sumstat}" .processed.txt)
	REF_BASENAME=$(basename "${ref_bed}" .bed)

        Rscript Rscripts/FUSION.assoc_test.nochr.R \
                --sumstats ${sumstat} \
                --force_model lasso \
                --weights wgtlist \
                --weights_dir ./ \
                --ref_ld_chr $REF_BASENAME \
                --out "$SS_BASENAME".results.dat
	
	"""
}
