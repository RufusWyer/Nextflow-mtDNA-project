#!/usr/bin/env nextflow

process calculateHsq {
    input:
    tuple path(pruned_bed), path(pruned_bim), path(pruned_fam), path(rac_pheno_log), path(fusion_covs_quant), path(fusion_covs_disc)

    output:
    path "hsq_results.HEreg"

    script:
	"""
	BASENAME=$(basename ${pruned_bed} .bed)

	software/gcta-1.94.4-linux-kernel-3-x86_64/gcta64 \
	  --bfile $BASENAME \
	  --make-grm \
	  --out pruned_grm

	software/gcta-1.94.4-linux-kernel-3-x86_64/gcta64 \
	  --grm        pruned_grm \
	  --pheno      rac_pheno_log  \
	  --covar       fusion_covs_disc \
	  --qcovar      fusion_covs_quant \
	  --HEreg \
	  --out outputs/fusion/hsq_results

	"""
}
