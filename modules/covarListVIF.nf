#!/usr/bin/env nextflow

process covarListVIF {
    input:
    tuple path(wbcovs_celltypes), path(pfs), path(gpca)

    output:
    path "covlist_VIF10.txt"

    script:
    """
    Rscript Rscripts/covarListVIF_NF.R $wbcovs_celltypes $pfs $gpca covlist_VIF10.txt
    """
}
