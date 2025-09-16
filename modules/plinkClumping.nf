#!/usr/bin/env nextflow

process plink_assoc {
    input:
    tuple path(maf5geno1_bed), path(maf5geno1_bim), path(maf5geno1_fam), path(assoc)

    output:
    path "ld1e-1.snplist", emit: snplist
    tuple path("ld1e-1.filtered.bed"), path("ld1e-1.filtered.bim"), path("ld1e-1.filtered.fam"), emit: filtered



    script:
	"""

	BASENAME=$(basename ${maf5geno1_bed} .bed)

	plink2 --bfile $BASENAME --clump ${assoc} \
	--clump-p1 0.1 \
	--clump-p2 0.5 \
	--clump-kb 250 \
	--out ld1e-1

	awk 'NR > 1 {print $3}' ld1e-1.clumps > ld1e-1.snplist

	plink2 --bfile $BASENAME --extract ld1e-1.snplist --make-bed --out ld1e-1.filtered
	
	"""
}

