#!/usr/bin/env nextflow

process rsidConversion {
    input:
    tuple path(train_bed), path(train_bim), path(train_fam), path(weights)

    output:
    path "weights.wgt.rsid.Rdat", emit: weights
    tuple path("train_filtered.rsid.bed"), path("train_filtered.rsid.bim"), path("train_filtered.rsid.fam"), emit: ref


    script:
    """
    Rscript Rscripts/rsidConversion_NF.R ${train_bim} ${weights} "train_snplist_rsid.txt" "train_snpmap_rsid.txt" "weights.wgt.rsid.Rdat"

    BASENAME=$(basename ${train_bed} .bed)

    plink --bfile $BASENAME \
        --extract train_snplist_rsid.txt \
        --make-bed --out train_filtered.ext

    plink --bfile train_filtered.ext \
        --update-name train_snpmap_rsid.txt \
        --make-bed --out train_filtered.rsid

    """
}
