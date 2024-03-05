
# README for GO Enrichment Analysis Pipeline

## Overview
This repository contains a Nextflow pipeline specifically designed for Gene Ontology (GO) Enrichment Analysis. It is tailored to analyze gene expression data, facilitating the identification of significantly enriched GO terms based on differential expression profiles.

## Usage
To execute the pipeline, several parameters must be provided:

1. `--meta_file`: Path to the metadata file containing information like experimental conditions.
2. `--count_file`: Path to the count file that includes gene expression data.
3. `--gene_column`: The column name in the count file that contains gene identifiers.
4. `--organism`: The organism for the genes (e.g., 'human', 'mouse').
5. `--logFC`: Boolean flag to apply a log fold change (logFC) threshold (default: `true`).
6. `--logFC_up`: Upper log2 fold change threshold for upregulated genes (default: `1`).
7. `--logFC_down`: Lower log2 fold change threshold for downregulated genes (default: `-1`).
8. `--p_adj`: Boolean flag to use adjusted p-values (default: `true`).
9. `--alpha`: Significance threshold for p-values (default: `0.05`).
10. `--output`: Directory for storing output files (default: `"./output/"`).

Example command:
```
./nextflow run main.nf --meta_file path/to/meta.txt --count_file path/to/count.txt --gene_column GeneSymbol --organism human --output output_directory
```

## Output Description
The GO Enrichment Analysis pipeline generates various outputs, including:

- **GO Enrichment Results**: Detailed results of the GO enrichment analysis, indicating significantly enriched GO terms across different biological processes, cellular components, and molecular functions.
- **Statistical Summaries**: Summaries of statistical analyses, including p-values, log fold changes, and adjusted p-values.
- **Visualization Files**: Graphical representations of enriched GO terms, aiding in the interpretation and communication of the results.

## Getting Started
To begin using this pipeline, clone the repository and ensure that Nextflow is installed on your system. Prepare your metadata and count files in the specified format. Customize the input parameters as needed and run the pipeline using the example command provided. The results will be outputted in the designated directory.

---
**Note**: Adjust paths, filenames, and parameters to align with the specific needs of your data and analysis.