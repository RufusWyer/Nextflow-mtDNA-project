#!/usr/bin/env nextflow

process gwasResults {
    publishDir "results/gwasResults", mode: 'copy'

    input:
    path results

    output:
    path results
}
