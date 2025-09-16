#!/usr/bin/env nextflow

process testTrainSplit {
    input:
    path pmc

    output:
    path "train_subset.txt", emit: train
    path "test_subset.txt", emit: test

    script:
    """
    Rscript Rscripts/testTrainSplit_NF.R $pmc train_subset.txt test_subset.txt
    """
}
