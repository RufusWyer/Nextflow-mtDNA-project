#!/usr/bin/env nextflow

process plinkMasterCovariates {
    input:
    tuple path(wbcovs_celltypes), path(pfs), path(gpca)

    output:
    path "plinkMasterCovariates.csv"

    script:
    """
    Rscript Rscripts/plinkMasterCovariates_NF.R $wbcovs_celltypes $pfs $gpca plinkMasterCovariates.csv
    """
}
