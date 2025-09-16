#!/usr/bin/env nextflow

process racProcessing {
    input:
    tuple path(counts), path(lengths)

    output:
    path "RNA_avcov.csv"

    script:
    """
    Rscript Rscripts/racProcessing_NF.R $counts $lengths RNA_avcov.csv
    """
}
