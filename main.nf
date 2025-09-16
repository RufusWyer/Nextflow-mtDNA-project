#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// params defined in config

// include modules
include { covarCelltypesProcessing } from './modules/covarCelltypesProcessing.nf'
include { racProcessing } from './modules/racProcessing.nf'
include { peerFactors } from './modules/peerFactors.nf'
include { plinkFiltering } from './modules/plinkFiltering.nf'
include { plinkMaf5geno1 } from './modules/plinkMaf5geno1.nf'
include { plinkGPCA } from './modules/plinkGPCA.nf'
include { plinkMasterCovariates } from './modules/plinkMasterCovariates.nf'
include { plinkPhenotypeProcessing } from './modules/plinkPhenotypeProcessing.nf'
include { covarListVIF } from './modules/covarListVIF.nf'
include { plinkAssoc } from './modules/plinkAssoc.nf'
include { gwasResults } from './modules/gwasResults.nf'
include { gwasPlots } from './modules/gwasPlots.nf'
include { testTrainSplit } from './modules/testTrainSplit.nf'
include { plinkClumping } from './modules/plinkClumping.nf'
include { fusionCovarReformatting } from './modules/fusionCovarReformatting.nf'
include { calculateHsq } from './modules/calculateHsq.nf'
include { fusionModelWeights } from './modules/fusionModelWeights.nf'
include { rsidConversion } from './modules/rsidConversion.nf'
include { sumstatsProcessing } from './modules/sumstatsProcessing.nf'
include { fusionAssoc } from './modules/fusionAssoc.nf'


workflow {

	// Processing GTEx metadata covariates and estimating whole blood cell type proportions (using RNAseq tpm data)
	covarCelltypesProcessing_in = Channel.of(tuple(file(params.tpm), file(params.samples), file(params.phs)))
	covarCelltypes_ch = covarCelltypesProcessing(covarCelltypesProcessing_in)

	// Calculating mitochondrial expression phenotype, mtRNA average coverage/nuclear RNA average coverage (RAC)
	racProcessing_in = Channel.of(tuple(file(params.counts), file(params.lengths)))
	RNA_avcov_ch = racProcessing(racProcessing_in)
	
	// Generate plink ID list of individuals in both RNAseq and WGS datasets
	wgs_ch = Channel.of(tuple(file(params.wgsbed), file(params.wgsbim), file(params.wgsfam)))
	plinkFiltering_in = covarCelltypes_ch.combine(wgs_ch)
        keepIds_ch = plinkFiltering(plinkFiltering_in)
	
	// Filter WGS file to ID subset, and filter to minor allele freq 5% and geno missingness 1%, keeping only common variants
	plinkMaf5geno1_in = fam_ch.combine(keepIds_ch)
	plinkMaf5geno1_ch = plinkMaf5geno1(plinkMaf5geno1_in)	
	plinkMaf5geno1_ch = plinkMaf5geno1_ch.map { bed, bim, fam -> tuple(bed, bim, fam) }
	
	// Calculate genetic principal components from processed WGS data, undergoing linkage-disequilibrium pruning first
	// Pruned WGS data saved for use later in heritability calculation
	plinkGPCA_ch = plinkGPCA(plinkMaf5geno1_ch)
	gpca_ch = plinkGPCA_ch.gpca
	pruned_ch = plinkGPCA_ch.pruned.map { bed, bim, fam -> tuple(bed, bim, fam) }

	// Calculating peerfactors of RNAseq tpm data, with mt RAC phenotype as a covariate to avoid capturing expression information, GPCA used for ID filtering
        peerFactors_script = Channel.of(file('Rscripts/peerFactors_mtCovariate_NF.R'))
        peerFactors_in = Channel.of(file(params.tpm)).combine(RNA_avcov_ch).combine(gpca_ch)
        peerFactors_ch = peerFactors(peerFactors_in, peerFactors_script)
	
	// Collating and formatting all covariates for use in plink GWAS models
	plinkMasterCovariates_in = covarCelltypes_ch.combine(peerFactors_ch).combine(gpca_ch)
	plinkMasterCovariates_ch = plinkMasterCovariates(plinkMasterCovariates_in)
	
	// Formatting RNA avcov phenotype for use in plink GWAS models
	// Both log-transformed and non-log-transformed variations, for use in different processes
	plinkPhenotypeProcessing_in = RNA_avcov_ch.combine(gpca_ch)
	plinkPhenotypeProcessing_ch = plinkPhenotypeProcessing(plinkPhenotypeProcessing_in)
	rac_log_ch   = plinkPhenotypeProcessing_ch.log
	rac_nolog_ch = plinkPhenotypeProcessing_ch.nolog

	// VIF trimming of covariate list to avoid overfitting, VIF threshold set to 10
	// Output is only list of covariate names, used as argument in plink GWAS model
	covarListVIF_in = covarCelltypes_ch.combine(peerFactors_ch).combine(gpca_ch)
	covarListVIF10_ch = plink_covs_VIF(plink_covs_VIF_in)

	// Running first plink GWAS test (GLM)
	// Results for this are used as overall GWAS results for phenotype. A second GWAS is run after splitting into testing/training sets to build the predictive ML model with Fusion
	plinkAssoc_in = plinkMafgeno1_ch.combine(plinkMasterCovariates_ch).combine(rac_nolog_ch).combine(covarListVIF10_ch)
	plink_assoc_ch = plink_assoc(plink_assoc_in)
	assocResults_ch = plinkAssoc_ch.assoc
	// Publishing results of initial GWAS to results/gwasResults
	gwasResults(plinkAssoc_ch.publish)

	// Creating Manhattan and QQ plots for GWAS model
	// Plots are published to results/gwasPlots
	gwasPlots(assocResults_ch)

	// Splitting IDs into testing and training set to build the predictive ML model with Fusion, and allow its validation.
	testTrainSplit_ch = testTrainSplit(plinkMasterCovariates_ch)
	testIds_ch = testTrainSplit_ch.test
	trainIds_ch = testTrainSplit_ch.train
	
	// Now have to repeat the above pre-processing and GWAS steps with the training set, to avoid data leakage into the training set model from the testing set
	// First, filter WGS file to id subset, and filter to minor allele freq 5% and geno missingness 1% as before
        trainMaf5geno1_in = fam_ch.combine(trainIds_ch)
        trainMaf5geno1_ch = plinkMaf5geno1(trainMaf5geno1_in)
        trainMaf5geno1_ch = trainMaf5geno1_ch.map { bed, bim, fam -> tuple(bed, bim, fam) }

	// Calculate genetic principal components from processed WGS data, undergoing linkage-disequilibrium pruning first
        trainGPCA_ch = plinkGPCA(trainMaf5geno1_ch)
	trainGPCA_ch = trainGPCA_ch.gpca

	// Calculating peerfactors of RNAseq tpm data, with mt RAC phenotype as a covariate to avoid capturing expression information, GPCA used for ID filtering
        trainPeerFactors_in = Channel.of(file(params.tpm)).combine(RNA_avcov_ch).combine(trainGPCA_ch)
        trainPeerFactors_ch = peerFactors(trainPeerFactors_in, peerFactors_script)

        // Collating and formatting all covariates for use in plink GWAS models
        trainMasterCovariates_in = covarCelltypes_ch.combine(trainPeerFactors_ch).combine(trainGPCA_ch)
        trainMasterCovariates_ch = plinkMasterCovariates(trainMasterCovariates_in)

        // Formatting RNA avcov phenotype for use in plink GWAS models
        // Both log-transformed and non-log-transformed variations, for use in different processes
        trainPhenotypeProcessing_in = RNA_avcov_ch.combine(trainGPCA_ch)
        trainPhenotypeProcessing_ch = plinkPhenotypeProcessing(trainPhenotypeProcessing_in)
        train_rac_log_ch   = rac_phenos_ch.log
        train_rac_nolog_ch = rac_phenos_ch.nolog

        // Running plink GWAS (GLM) with training set, using original covariate list (no chance of data leakage but likely to perform better)
        trainPlinkAssoc_in = trainMafgeno1_ch.combine(trainMasterCovariates_ch).combine(train_rac_nolog_ch).combine(covarListVIF10_ch)
        trainPlinkAssoc_ch = plinkAssoc(train_plink_assoc_in)
        trainPlinkAssoc_ch = trainPlinkAssoc_ch.assoc
	
	// Linkage-disequilibrium clumping of the training set using the pvalue of the training set GWAS, eliminating both LD-based issues and poorly associated SNPs (pval threshold = 0.1)
	plinkClumping_in = trainMafgeno1_ch.combine(trainPlinkAssoc_ch)
	plinkClumping_ch = plinkClumping(plinkClumping_in)
	clumpedSnplist_ch = plinkClumping.snplist
	clumpedFiltered_ch = plinkClumping.filtered
	clumpedFiltered_ch = clumpedFiltered_ch.map { bed, bim, fam -> tuple(bed, bim, fam) }

	// Reformatting training set covariates for use in Fusion and GCTA fro heritability estimate. GCTA also requires covars to be split into discrete and quantitative files.
	// GCTA Hsq calculation can be done with entire dataset, but Fusion requires just training set data, so this must be repeated once for each set.
	fusionCovarReformatting_in = plinkMasterCovariates_ch.combine(covarListVIF10_ch)
	fusionCovarReformatting_ch = fusionCovarReformatting(fusionCovarReformatting_in)
	fusionCovar_quant = fusionCovarReformatting_ch.quant
	fusionCovar_disc = fusionCovarReformatting_ch.disc
	fusionCovar = fusionCovarReformatting_ch.covs
	
	trainFusionCovarReformatting_in = trainMasterCovariates_ch.combine(covarListVIF10_ch)
	trainFusionCovarReformatting_ch = trainFusionCovarReformatting(trainFusionCovarReformatting_in)
        trainFusionCovar = trainFusionCovarReformatting_ch.covs

	// Calculating heritability (Hsq) for use as input in fusion, using reformatted full-set covariates, split into quantitative and discrete covars.
	// Haseman-Elston Regression (HEreg) is used here over other methods due to relatively low sample size
	calculateHsq_in = pruned_ch.combine(rac_log_ch).combine(fusionCovar_quant).combine(fusionCovar_disc)
	calculateHsq_ch = calculateHsq(calculateHsq_in)

	// The predictive ML model is built using LASSO regression in Fusion, taking covariates into account by linearly regressing them out of the phenotype
	// The model is then tested by using genotype sum-scoring the testing set WGS data, and running a pearson's correlation test on the Z-score standardised results versus the true testing set phenotype (also with covariates regressed out).
	// Performance results for the test are published to results/fusionModelWeights
	fusionModelWeights_in = clumpedFiltered_ch.combine(wgs_ch).combine(testIds_ch).combine(calculateHsq_ch).combine(train_rac_log_ch).combine(rac_log_ch).combine(trainFusionCovar).combine(fusionCovar)
	fusionModelWeights_ch = fusionModelWeights(fusionModelWeights_in)
	fusionModelWeights_ch = fusionModelWeights_ch.weights

	// Now the model is created, it can be used in a series of TWAS tests against common diseases, using GWAS summary statistics
	// The ML model SNPs are converted into RSID formats, to enable compatabiltiy with sumstats
	// The SNPs from the model are used later as a reference genome for the Fusion TWAS tests
	rsidConversion_in = clumpedFiltered_ch.combine(fusionModelWeights_ch)
	rsidConversion_ch = rsidConversion(rsidConversion_in)
	rsidWeights_ch = rsidConversion_ch.weights
	rsidRef_ch = rsidConversion_ch.ref
	rsidRef_ch = rsidRef_ch.map { bed, bim, fam -> tuple(bed, bim, fam) }

	
	// Download and process the sumstats files to make them compatible, uniform, and error-free
	// These are split into whole-blood and non-whole-blood GWAS sumstats 
	sumstatsProcessing_ch = sumstatsProcessing()
	wbSumstats = sumstatsProcessing_ch.wb
	nonwbSumstats = sumstatsProcessing_ch.nonwb

	// Running Fusion TWAS test, once with whole-blood GWAS sumstats, and then with non-whole-blood GWAS sumstats
	// Results for each sumstat file are published to results/fusionAssoc
	fusionAssoc(wbSumstats, rsidWeights_ch, rsidRef_ch)
	fusionAssoc(nonwbSumstats, rsidWeights_ch, rsidRef_ch)
	
}

