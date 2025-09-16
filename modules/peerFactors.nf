#!/usr/bin/env nextflow

process peerFactors {
	container 'singularity/r-peer_1.0.0.sif'

	input:
	tuple path(tpm), path(RNA_avcov), path(gpca)
	path script

	output:
	path "mtcPF.factors.csv"

	script:
	"""
        Rscript $script $tpm $RNA_avcov $gpca mtcPF.factors.csv


	"""
}
