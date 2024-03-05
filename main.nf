#!/usr/bin/env nextflow

params.output = "./output/"
params.meta_file = '' // Path to meta.txt
params.count_file = '' // Path to count.txt

params.gene_column = '' // Gene column with gene names
params.organism = '' // Organsim of genes in gene column

params.logFC = true
params.logFC_up = 1
params.logFC_down = -1
params.p_adj = true
params.alpha = 0.05

// scripts
go_enrichment_script = Channel.fromPath("${projectDir}/GOEnrichment.R")

process go_enrichment {
    container 'kadam0/goenrichment:0.0.1'
    publishDir params.output, mode: "copy"

    input:
    path script_file from go_enrichment_script
    path meta_file from Channel.fromPath(params.meta_file)
    path count_file from Channel.fromPath(params.count_file)

    output:                                
    path "*"

    script:
    """
    Rscript $script_file --meta_file ${meta_file} --gene_column ${params.gene_column} --organism ${params.organism} --count_file ${count_file} --out_dir ${params.output} --logFC ${params.logFC} --logFC_up ${params.logFC_up} --logFC_down ${params.logFC_down} --p_adj ${params.p_adj} --alpha ${params.alpha}
    """
}

workflow {
  go_enrichment()
}
