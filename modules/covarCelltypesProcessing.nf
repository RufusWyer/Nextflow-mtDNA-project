#!/usr/bin/env nextflow

process covarCelltypesProcessing {
    input:
    tuple path(tpm), path(samples), path(phs)

    output:
    path "wb_covs_celltypes.csv"

    script:
    """
    Rscript Rscripts/covarCelltypesProcessing_NF.R $tpm $samples $phs wb_covs_celltypes.csv
    """
}
