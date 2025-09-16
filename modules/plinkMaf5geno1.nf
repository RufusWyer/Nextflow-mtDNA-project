#!/usr/bin/env nextflow

process plinkMaf5geno1 {
    input:
    tuple path(wgs_bed), path(wgs_bim), path(wgs_fam), path(keep_ids)

    output:
    tuple path("maf5geno1.bed"), path("maf5geno1.bim"), path("maf5geno1.fam")

    script:
    """
    BASENAME=$(basename ${wgs_bed} .bed)
    plink --bfile $BASENAME --keep ${keep_ids} --geno 0.01 --maf 0.05 --nonfounders --allow-no-sex --make-bed --out #maf5geno1
    """
}
