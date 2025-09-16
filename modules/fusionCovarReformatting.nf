#!/usr/bin/env nextflow

process fusionCovarReformatting {
    input:
    tuple path(pmc), path(covlist)

    output:
    path "fusion_covs.txt", emit: covs
    path "fusion_covs_quant.txt", emit: quant
    path "fusion_covs_disc.txt", emit: disc

    script:
    """
    Rscript Rscripts/fusionCovarReformatting_NF.R $pmc $covlist fusion_covs.txt fusion_covs_quant.txt fusion_covs_disc.txt
    """
}
