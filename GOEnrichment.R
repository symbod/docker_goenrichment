#!/usr/bin/env Rscript

## Script name: GOEnrichment.R
##
## Purpose of script: GO Enrichment 
##
## Author: Klaudia Adamowicz
## Email: klaudia.adamowicz@uni-hamburg.de
##
## Date Created: 2024-03-04
##
## Copyright (c) Dr. Tanja Laske, 2024

# Load required libraries --------------------------
suppressPackageStartupMessages({
  # List of required CRAN and Bioconductor packages
  required_cran_packages <- c("optparse", "argparse", "data.table", "tibble", 
                              "igraph", "tidyverse", "rjson", "ggplot2", 
                              "plotly", "htmlwidgets", "dplyr")
  required_bioc_packages <- c("limma", "clusterProfiler", 
                              "org.Rn.eg.db", "org.Hs.eg.db")
  
  # Function to check and install CRAN packages
  check_install_cran_packages <- function(packages) {
    for (package in packages) {
      if (!require(package, character.only = TRUE)) {
        install.packages(package)
      }
      library(package, character.only = TRUE)
    }
  }
  
  # Function to check and install Bioconductor packages
  check_install_bioc_packages <- function(packages) {
    for (package in packages) {
      if (!require(package, character.only = TRUE, quietly = TRUE)) {
        BiocManager::install(package)
      }
      library(package, character.only = TRUE)
    }
  }
  
  # Check and install required packages
  check_install_cran_packages(required_cran_packages)
  check_install_bioc_packages(required_bioc_packages)
})


# Define Methods --------------------------

#' Method performing limma
#'
#' @param set_condition Vector of experimental design specifying the condition(s)
#' @param set_counts Counts (rows = proteins, columns = samples)
#' @param set_comparisons Vector of comparisons
#'
#' @return Results of the limma analysis
#'
#' @export
perform_de <- function(set_condition, set_counts, set_comparisons){
  # create design matrix
  groupsM <- as.factor(set_condition)
  designM <- model.matrix(~0+groupsM) 
  colnames(designM) <-levels(groupsM) 
  fit <- lmFit(set_counts, designM)
  
  # create contrasts
  contr <- makeContrasts(contrasts = set_comparisons, levels = colnames(coef(fit)))
  
  fit2 <- contrasts.fit(fit, contr)
  ebfit <- eBayes(fit2, trend = TRUE)
  return(ebfit)
}

#' Method to extract results of fit object
#'
#' @param set_fit Fit object of the perform_de method
#' @param set_comparisons Vector of comparisons
#' @param out_dir Output directory
#' @param lfc_up Log fold change threshold (upregulated)
#' @param lfc_down Log fold change threshold (downregulated)
#' @param alpha p-value or p.adjust threshold (depending on padj)
#' @param strict TRUE if < and > should be used for threshold comparison, FALSE if <= and >=
#' @param padj TRUE if alpha should be applied to padj, FALSE if alpha applied to p-value
#' @param logFC_thr TRUE if a threshold for logFC should be set, FALSE if not
#'
#' @return Extracted results from the fit object
#'
#' @export
compare_de_expr <- function(set_fit, set_comparisons, out_dir, lfc_up = args$logFC_up, lfc_down = args$logFC_down, alpha = args$alpha, 
                            strict = FALSE, padj = args$p_adj, logFC_thr=args$logFC, write_output = TRUE){
  de_results <- list()
  for (i in 1:length(set_comparisons)){
    # get results for each comparison
    top.table <- topTable(set_fit, sort.by = "P", number=Inf, coef=c(i))
    gene_reg <- setDT(top.table, keep.rownames = "gene") # save row names as column
    
    # different threshold settings
    if (logFC_thr){
      if (strict){
        if (padj){
          gene_reg$Change <- ifelse(gene_reg$logFC > lfc_up & gene_reg$adj.P.Val < alpha, 1, ifelse(gene_reg$logFC < lfc_down & gene_reg$adj.P.Val < alpha, 1, 0))
        } else {
          gene_reg$Change <- ifelse(gene_reg$logFC > lfc_up & gene_reg$P.Value < alpha, 1, ifelse(gene_reg$logFC < lfc_down & gene_reg$P.Value < alpha, 1, 0))
        }
      } else {
        if (padj){
          gene_reg$Change <- ifelse(gene_reg$logFC >= lfc_up & gene_reg$adj.P.Val <= alpha, 1, ifelse(gene_reg$logFC <= lfc_down & gene_reg$adj.P.Val <= alpha, 1, 0))
        } else {
          gene_reg$Change <- ifelse(gene_reg$logFC >= lfc_up & gene_reg$P.Value <= alpha, 1, ifelse(gene_reg$logFC <= lfc_down & gene_reg$P.Value <= alpha, 1, 0))
        }
      }
    } else {
      if (strict){
        if (padj){
          gene_reg$Change <- ifelse(gene_reg$adj.P.Val < alpha, 1, 0)
        } else {
          gene_reg$Change <- ifelse(gene_reg$P.Value < alpha, 1, 0)
        }
      } else {
        if (padj){
          gene_reg$Change <- ifelse(gene_reg$adj.P.Val <= alpha, 1, 0)
        } else {
          gene_reg$Change <- ifelse(gene_reg$P.Value <= alpha, 1, 0)
        }
      }
    }
    if (write_output){
      write.table(gene_reg, file = paste0(out_dir,"de_data/","DEdata.", set_comparisons[[i]],".txt"),sep=" ",row.names = FALSE) 
    }
    de_results[[set_comparisons[[i]]]] <- gene_reg
  }
  return(de_results)
}

#' Prepare matrix for differential expression analysis
#'
#' @param counts Counts (rows = proteins, columns = samples)
#' @param md Metadata
#' @param condition_name Name of the condition in metadata
#' @param comparisons Vector of comparisons
#' @param out_dir Output directory
#' @param network_proteins Proteins in the network
#' @param lfc_up Log fold change threshold for upregulated proteins
#' @param lfc_down Log fold change threshold for downregulated proteins
#' @param alpha p-value or p.adjust threshold (depending on padj)
#' @param strict TRUE if < and > should be used for threshold comparison, FALSE if <= and >=
#' @param padj TRUE if alpha should be applied to padj, FALSE if alpha applied to p-value
#' @param logFC_thr TRUE if a threshold for logFC should be set, FALSE if not
#' @param plot TRUE to generate plots, FALSE otherwise
#' @param write_output TRUE to write output files, FALSE otherwise
#'
#' @return Prepared matrix for differential expression analysis
#'
#' @export
prepare_matrix <- function(counts, md, condition_name, comparisons, out_dir,
                           lfc_up = args$logFC_up, lfc_down = args$logFC_down, alpha = args$alpha, 
                           strict = FALSE, padj = args$p_adj, logFC_thr=args$logFC, plot = TRUE, write_output = TRUE){
  if (write_output){
    # create output directories
    dir.create(out_dir, showWarnings = FALSE) #stops warnings if folder already exists
    dir.create(file.path(out_dir,"de_data"), showWarnings = FALSE) #stops warnings if folder already exists
    
  }
  
  # get condition vector
  condition <- md[, condition_name]
  
  # perform DE analysis pairwise
  pairwise_de <- perform_de(set_condition = condition, set_counts = counts, set_comparisons = comparisons)
  # extract results
  de_results <- compare_de_expr(set_fit = pairwise_de, set_comparisons = comparisons, out_dir = out_dir, lfc_up = lfc_up, 
                                lfc_down = lfc_down, alpha = alpha, strict = strict, padj = padj, logFC_thr = logFC_thr, write_output = write_output)
  mat = data.frame(matrix(ncol=1,nrow=0, dimnames=list(NULL, c("gene"))))
  for (comp in comparisons) {
    cur_df <- de_results[[comp]][,c("gene","Change"), drop = FALSE]
    cur_df <- cur_df[!is.na(cur_df$Change), ]
    colnames(cur_df) <- c("gene", comp)
    mat = merge(mat, cur_df, all=TRUE)
  }
  
  # Fill NA with 0 
  mat[is.na(mat)] <- 0
  
  # filter for only expressed rows 
  mat <- mat[rowSums(mat==0) != (ncol(mat)-1), ]
  # explode id column
  mat <- mat %>% mutate(gene=strsplit(gene, ";")) %>% unnest(gene)
  return(mat)
}


protein_to_gene_mapping <- function(count_data, ind_mat, gene_col) {
  # Split and explode the Protein.IDs column
  long_format <- count_data[, .(Protein.ID = unlist(strsplit(Protein.IDs, ";"))), by = .(Gene = get(gene_col))]
  
  # Remove duplicate rows
  long_format <- unique(long_format)
  
  # Rename column in ind_mat
  setnames(ind_mat, old = "gene", new = "Protein.ID")
  
  # Merge with ind_mat
  merged_data <- merge(ind_mat, long_format, by = "Protein.ID", all.x = TRUE)
  
  return(merged_data)
}


save_meta_data <- function(comp, cond, meta_data, out_dir, filename_prefix) {
  meta_all <- list()
 
  sample_one <- strsplit(comp, split="-")[[1]][1]
  sample_two <- strsplit(comp, split="-")[[1]][2]
  
  # save meta data into json 
  sample_group <- if (cond == "TimeCond") list("Timepoint", "Condition") else if (cond == "TimeCondLoc") list("Timepoint", "Condition", "Location") else list("Condition")
  meta_all[[comp]] = list(
    samples_group_A = meta_data[meta_data[[cond]] %in% c(sample_one)]$Column_name,
    samples_group_B = meta_data[meta_data[[cond]] %in% c(sample_two)]$Column_name,
    group_A = str_split(sample_one,"_")[[1]],
    group_B = str_split(sample_two,"_")[[1]],
    sample_group = sample_group
  )
 
  # Add links to the entire meta_all list
  meta_all_links <- list(
    druglist = file.path(out_dir, paste0(filename_prefix, "_drugs.csv")),
    genelist = file.path(out_dir, paste0(filename_prefix, "_genelist.csv")),
    graph_df = file.path(out_dir, paste0(filename_prefix, "_graph.csv")),
    graph_file = file.path(out_dir, paste0(filename_prefix, "_network.graphml")))
  
  # Combine metadata with links
  meta_all_combined <- list(meta_data = meta_all, links = meta_all_links)
  # Write to JSON
  write(toJSON(meta_all_combined), file.path(out_dir, paste0(filename_prefix,".txt.json")))
}

## ------------- Prepare Network Enrichment ------------------

### Parse arguments -----

# set up arguments
parser <- OptionParser()
parser <- add_option(parser, c("-m","--meta_file"), help="Meta data description file")
parser <- add_option(parser, c("-c","--count_file"), help="Preprocessed count file")
parser <- add_option(parser, c("-o","--out_dir"), help="Directory for output files", default="")

# Adding option for column selection

parser <- add_option(parser, c("-g","--gene_column"), help="Gene column with gene names.")
parser <- add_option(parser, c("--organism"), help="Organsim of genes in gene column.")

# Adding new options for thresholds with defaults
parser <- add_option(parser, c("--logFC"), help = "Boolean specifying whether to apply a logFC threshold (TRUE) or not (FALSE)", type = "logical", default = TRUE)
parser <- add_option(parser, c("--logFC_up"), help = "Upper log2 fold change threshold (dividing into upregulated)", type = "numeric", default = 1)
parser <- add_option(parser, c("--logFC_down"), help = "Lower log2 fold change threshold (dividing into downregulated)", type = "numeric", default = -1)
parser <- add_option(parser, c("--p_adj"), help = "Boolean specifying whether to apply a threshold on adjusted p-values (TRUE) or on raw p-values (FALSE)", type = "logical", default = TRUE)
parser <- add_option(parser, c("--alpha"), help = "Threshold for adjusted p-values or p-values", type = "numeric", default = 0.05)

# get command line options, if help option encountered print help and exit
args <- parse_args(parser)
# check if mandatory options are set
check_options <- function(tags){
  for (tag in tags){
    if (is.null(args[tag])){
      print_help(parser)
      stop("Missing mandatory option.", call.=FALSE)
    }
  }
}
check_options(c('meta_file','count_file','gene_column'))

# save arguments
meta_file_path <- args$meta_file
count_file_path <- args$count_file
out_dir <- args$out_dir
gene_column <- args$gene_column
organism <- args$organism

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) #stops warnings if folder already exists

#### Load data --------
meta_data <- fread(meta_file_path)
count_data <- fread(count_file_path)

## Prepare data ----

#### Correct data ----
ID_column <- "Animal"

# remove ref
meta_data <- meta_data[meta_data[[ID_column]] != "ref",]

# convert timepoint column
meta_data[, Timepoint := sapply(Timepoint, function(tp) ifelse(grepl("pre", tp), -as.numeric(gsub("pre", "", tp)), ifelse(grepl("post", tp), as.numeric(gsub("post", "", tp)), as.numeric(tp))))]
meta_data[, Sample_name := if ("Sample_name" %in% names(meta_data)) Sample_name else if ("Label" %in% names(meta_data)) Label else NULL]
meta_data[, Column_name := if ("Column_name" %in% names(meta_data)) Column_name else if ("Column" %in% names(meta_data)) Column else NULL]
#meta_data[, Timepoint := as.numeric(Timepoint)]

### Rename Columns ---

# rename from file name to sample name
names(count_data) <- plyr::mapvalues(names(count_data), from = meta_data$Column_name, to = meta_data$Sample_name, warn_missing=FALSE) 
meta_data <- subset(meta_data, meta_data$Sample_name %in% names(count_data))

# remove columns that are not in meta_data
columns_to_keep <- c("Protein.IDs", gene_column, meta_data$Sample_name)
existing_columns <- columns_to_keep[columns_to_keep %in% names(count_data)]
count_data <- count_data[, existing_columns, with = FALSE]

#### Encode Columns ----

mappings <- list(
  Condition = list(setNames(paste0("C", seq_along(unique(meta_data$Condition))), unique(meta_data$Condition))),
  Location = if("Location" %in% names(meta_data) && length(unique(meta_data$Location)) > 1) 
    list(setNames(paste0("L", seq_along(unique(meta_data$Location))), unique(meta_data$Location))),
  Timepoint = list(setNames(paste0("T", seq_along(sort(unique(meta_data$Timepoint)))), sort(unique(meta_data$Timepoint))))
)

write(toJSON(mappings), file.path(out_dir, "Metadata_encoding.txt.json")) 

reverse_mappings <- list(
  Condition = if (!is.null(mappings$Condition)) setNames(names(mappings$Condition[[1]]), unlist(mappings$Condition[[1]])) else NULL,
  Location = if (!is.null(mappings$Location)) setNames(names(mappings$Location[[1]]), unlist(mappings$Location[[1]])) else NULL,
  Timepoint = if (!is.null(mappings$Timepoint)) setNames(names(mappings$Timepoint[[1]]), unlist(mappings$Timepoint[[1]])) else NULL
)

write(toJSON(reverse_mappings), file.path(out_dir, "Metadata_reverse_encoding.txt.json")) 

# rename condition
meta_data[, Condition := mappings$Condition[[1]][Condition]]
# rename timepoint
meta_data[, Timepoint := as.character(Timepoint)]
meta_data[, Timepoint := mappings$Timepoint[[1]][Timepoint]]
# rename location
if(!is.null(mappings$Location)) {
  meta_data[, Location := mappings$Location[[1]][Location]]
}

### Create indicator matrix ----
#### Prepare input data ----

# create count matrix 
counts <- as.data.frame(count_data[,c(meta_data$Sample_name), with=F], )
rownames(counts) <- count_data$Protein.IDs
counts <- counts[rowSums(is.na(counts)) != ncol(counts), ] # remove where full row is NA


#### Assign active and inactive cases ----
##### Conditions ----
###### Save comparisons ----
comparisons <- c()
elements <- sort(unique(meta_data$Condition))
for (index_a in 1:(length(elements)-1)){
  for (index_b in (index_a+1):length(elements)){
    comparisons <- c(comparisons, paste0(elements[index_a], "-", elements[index_b]))
  }
}
ind_mat_cond <- prepare_matrix(counts=counts, md=as.data.frame(meta_data), condition_name="Condition", 
                               comparisons=comparisons, out_dir = out_dir, write_output = FALSE)
ind_mat_cond <- protein_to_gene_mapping(count_data = count_data, ind_mat = ind_mat_cond, gene_col = gene_column)
write.table(ind_mat_cond, file.path(out_dir, "indicator_matrix_cond.tsv"), sep="\t", row.names=FALSE)


##### Conditions and Timepoints ----
# save sample groups
target_column <- ifelse("Location" %in% names(meta_data) && length(unique(meta_data$Location)) > 1, 
                        "TimeCondLoc", 
                        "TimeCond")
# Create the target column based on the presence of 'Location'
meta_data[, (target_column) := if(target_column == "TimeCondLoc") {
  paste(meta_data$Timepoint, meta_data$Condition, meta_data$Location, sep = "_")} else {
    paste(meta_data$Timepoint, meta_data$Condition, sep = "_")}]

###### Save comparisons ----
comparisons <- c()
elements <- sort(unique(meta_data[[target_column]]))
for (index_a in 1:(length(elements)-1)){
  for (index_b in (index_a+1):length(elements)){
    # only if time point or condition are the same!
    split_a = str_split(elements[index_a],"_")[[1]]
    split_b = str_split(elements[index_b],"_")[[1]]
    # Count the number of differences
    differences <- sum(split_a != split_b)
    # Only add to comparisons if exactly one of timepoint, condition, or location is different
    if (differences == 1){
      comparisons <- c(comparisons, paste0(elements[index_a], "-", elements[index_b]))
    }
  }
}

ind_mat_timecond <- prepare_matrix(counts=counts, md=as.data.frame(meta_data), condition_name=target_column, 
                                   comparisons=comparisons, out_dir = out_dir, write_output = FALSE)
ind_mat_timecond <- protein_to_gene_mapping(count_data = count_data, ind_mat = ind_mat_timecond, gene_col = gene_column)
write.table(ind_mat_timecond, file.path(out_dir,"indicator_matrix_timecond.tsv"), sep="\t", row.names=FALSE)


# Run GO Enrichment -----------------------

## Functions ----

perform_GO_enrichment <- function(genes, organism, case, out_dir) {
  # Choose the correct OrgDb based on the organism
  OrgDb <- switch(organism,
                  "human" = org.Hs.eg.db,
                  "rat" = org.Rn.eg.db,
                  stop("Unsupported organism. Please use 'human' or 'rat'."))
      
  # Perform GO enrichment analysis for all ontologies
  ego <- enrichGO(gene = genes,
                  OrgDb = OrgDb,
                  keyType = "SYMBOL",
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  
  # Check results and create visualizations
  if (nrow(ego) > 0) {
    # Save results to a CSV file
    result_file <- file.path(out_dir, paste0("GO_Enrichment_", organism, "_", case, ".csv"))
    write.csv(ego, result_file, row.names = FALSE)
    
    dotplot_ego <- dotplot(ego) + facet_grid(~ONTOLOGY)
    # Convert ggplot object to an interactive plotly object
    p_interactive <- ggplotly(dotplot_ego)
    
    # Save the interactive plot to an HTML file
    plot_file <- file.path(out_dir, paste0("GO_Enrichment_DotPlot_", organism, "_", case, ".html"))
    saveWidget(p_interactive, file = plot_file)
    
    return(ego)
  } else {
    warning(paste0("No significant GO terms found for ", case))
    return(NULL)
  }
}

## condition vs condition ----

go_results <- list()

if (nrow(ind_mat_cond) > 0) {
  unique_genes <- unique(unlist(strsplit(na.omit(gsub("^$|^NA$", "", ind_mat_cond$Gene)), ";")))
  go_result <- perform_GO_enrichment(genes = unique_genes, organism = organism,
                                     case = names(ind_mat_cond)[2], out_dir = out_dir)
  
  go_results[["cond"]][[names(ind_mat_cond)[2]]] <- list("go_result" = go_result)
  save_meta_data(comp = names(ind_mat_cond)[2], cond = "Condition", meta_data = meta_data,
                 out_dir = out_dir, filename_prefix = "cond")
  
} else {
  warning("No data found. Skipping condition vs condition analysis.")
}
## time point condition vs time point condition ----

go_results[["tp"]] <- list()

for ( case in setdiff(names(ind_mat_timecond), c("Protein.ID", "Gene"))){
  if (sum(ind_mat_timecond[[case]]) > 0) {
    df <- ind_mat_timecond[c("Protein.ID", "Gene", case)]  %>% filter(.data[[case]] != 0)
    unique_genes <- unique(unlist(strsplit(na.omit(gsub("^$|^NA$", "", df$Gene)), ";")))
    if ( length(unique_genes) > 0) {
      go_result <- perform_GO_enrichment(genes = unique_genes, organism = organism, 
                                         case = case, out_dir = out_dir)
      
      go_results[["tp"]][[case]] <- list("go_result" = go_result)
    } 
    if (!is.null(go_results[["tp"]][[case]])){
      save_meta_data(comp = case, cond = "Condition", meta_data = meta_data,
                     out_dir = out_dir, filename_prefix = paste0("timepoint_",case))
    }
  } else {
    warning(paste0("No data found. Skipping ", case, " analysis."))
  }
}

# Overview statistics ----

# Define a threshold for significance, e.g., adjusted p-value < 0.05
significance_threshold <- 0.05

# Initialize data frames and list outside of the loop for efficiency
go_term_counts <- data.frame(comparison = character(), sig_go_terms = integer())
unique_go_terms <- data.frame(GO_ID = character(), ontology = character(), descritpion = character(), comparison = character())
gene_go_list <- list()

# Loop through each main part and sub-part of go_result
for (main_part in names(go_results)) {
  for (sub_part in names(go_results[[main_part]])) {
    go_result <- as.data.frame(go_results[[main_part]][[sub_part]]$go_result)
    
    # Check for valid enrichResult
    if (!is.null(go_result)) {
      # Get significant GO terms
      significant_go <- go_result[go_result$p.adjust < significance_threshold, ]
      sig_count <- nrow(significant_go)
      
      # Append to go_term_counts
      go_term_counts <- rbind(go_term_counts, data.frame(sub_part, sig_go_terms = sig_count))
      
      # Append to unique_go_terms
      if (sig_count > 0) {
        temp_df <- data.frame(GO_ID = as.character(significant_go$ID), 
                              ontology = as.character(significant_go$ONTOLOGY),
                              description = as.character(significant_go$Description), 
                              comparison = rep(sub_part, sig_count))
        unique_go_terms <- rbind(unique_go_terms, temp_df)
        # Update gene_go_list
        significant_genes <- strsplit(as.character(significant_go$geneID), "/")  # Splitting gene IDs
        for (i in seq_along(significant_genes)) {
          genes <- significant_genes[[i]]
          for (gene in genes) {
            if (!is.null(gene_go_list[[gene]])) {
              gene_go_list[[gene]] <- unique(c(gene_go_list[[gene]], significant_go$ID[i]))
            } else {
              gene_go_list[[gene]] <- significant_go$ID[i]
            }
          }
        }
      }
    } else {
      # Append zero count to go_term_counts
      go_term_counts <- rbind(go_term_counts, data.frame(main_part, sub_part, sig_go_terms = 0))
    }
  }
}

# Summarize unique GO terms
unique_go_summary <- unique_go_terms %>%
  group_by(GO_ID, Ontology, Description) %>%
  summarize(Subpart_Count = n_distinct(Subpart))

# Initialize a list to store data frames for each gene
gene_annotation_counts_list <- vector("list", length(gene_go_list))

# Iterate over the genes and store each result in the list
for (i in seq_along(gene_go_list)) {
  gene <- names(gene_go_list)[i]
  gene_annotations <- gene_go_list[[gene]]
  gene_annotation_counts_list[[i]] <- data.frame(
    Gene = gene,
    Comparison_Count = length(gene_annotations),
    Unique_GO_Terms = length(unique(gene_annotations))
  )
}

# Combine all data frames into one
gene_annotation_counts <- do.call(rbind, gene_annotation_counts_list)


## Top comparisons ----

go_term_counts <- go_term_counts[order(-go_term_counts$sig_go_terms),]

write.table(go_term_counts, file.path(out_dir,"top-comparisons.tsv"), sep="\t", row.names=FALSE)

## Top IDs ----

unique_go_summary <- unique_go_summary[order(-unique_go_summary$Subpart_Count),]

write.table(unique_go_summary, file.path(out_dir,"top-gos.tsv"), sep="\t", row.names=FALSE)

gene_annotation_counts <- gene_annotation_counts[order(-gene_annotation_counts$Comparison_Count),]

write.table(gene_annotation_counts, file.path(out_dir,"top-ids.tsv"), sep="\t", row.names=FALSE)
