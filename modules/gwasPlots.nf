#!/usr/bin/env nextflow

process gwasPlots {
    publishDir "results/gwasPlots", mode: 'copy'

    input:
    path glm

    output:
    path(mtDNA.manPlot.png)
    path(mtDNA.qqPlot.png)

    script:
	"""
        Rscript Rscripts/gwasPlots_NF.R "${glm}" "mtDNA.manPlot.png" "mtDNA.qqPlot.png"
	
	"""
}
