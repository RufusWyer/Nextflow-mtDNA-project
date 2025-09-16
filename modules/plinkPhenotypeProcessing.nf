#!/usr/bin/env nextflow

process plinkPhenotypeProcessing {
    input:
    tuple path(RNA_avcov), path(gpca)

    output:
    path "rac_pheno_log.txt" emit: log
    path "rac_pheno_nolog.txt", emit: nolog

    script:
    """
    Rscript Rscripts/plinkPhenotypeProcessing_NF.R $RNA_avcov $gpca rac_pheno_log.txt rac_pheno_nolog.txt
    """
}
