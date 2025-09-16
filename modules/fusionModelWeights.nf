#!/usr/bin/env nextflow

process fusionModelWeights {
    publishDir "results/fusionModelWeights", mode: 'copy'

    input:
    tuple path(clumped_train_bed), path(clumped_train_bim), path(clumped_train_fam), path(wgs_bed), path(wgs_bim), path(wgs_fam), path(test_ids), path(hsq), path(pheno_train), path(pheno_all), path(covs_train), path(covs_all)

    output:
    path "weights.wgt.Rdat", emit: weights
    path "fusionPredictionsResults.txt", publish: true

    script:
	"""
	// extract hsq result
	HSQ=$(awk 'NR==4 {print $2}' ${hsq})

	// set extract basenames from staged bed file trios, both training and testing sets
	TRAIN_BASENAME=$(basename ${clumped_train_bed} .bed)
	WGS_BASENAME=$(basename ${wgs_bed} .bed)
	
	// create a temp directory for fusion to work in
	mkdir fusion_tmp
	
	Rscript software/fusion/fusion_twas-master/FUSION.compute_weights.R \
	--bfile $"TRAIN_BASENAME".filtered \
	--tmp fusion_tmp/  \
	--out weights \
	--pheno ${pheno_train} \
	--crossval 0 \
	--hsq_set $HSQ \
	--covar ${covs_train} \
	--models lasso \
	
	// remove the temp directory fusion was working in
	rm -r fusion_tmp
	
	// convert the weights from an .Rdat file to a more usable and readable txt file
	Rscript Rscripts/fusionWeightsProcessing_NF.R weights.wgt.Rdat "weights.processed.txt"
	
	// filter the original WGS data to match the testing set IDs and the same SNP list as the model
	plink2 --bfile $WGS_BASENAME --keep ${test_ids} --extract ${snplist} --make-bed --out test.filtered
	
	plink \
	  --bfile test.filtered \
	  --score weights.processed.txt sum \
	  --out weights.results



	Rscript Rscripts/fusionModelPredictions_NF.R weights.processed.txt ${covs_all} ${pheno_all} ${test_ids} "fusionPredictionsResults.txt"


	"""
}
