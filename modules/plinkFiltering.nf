#!/usr/bin/env nextflow

process plinkFiltering {
    input:
    tuple path(wbcovs_celltypes), path(fam)

    output:
    path "keep_ids.txt"

    script:
    """
    Rscript Rscripts/plink_filtering_NF.R $wbcovs_celltypes $fam keep_ids.txt
    """
}
