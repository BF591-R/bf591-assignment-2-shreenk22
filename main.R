#!/usr/bin/Rscript
## Author: Taylor Falk
## tfalk@bu.edu
## BU BF591
## Assignment Week 2

#### Bioconductor ####
# it is standard among R packages to define libraries and packages at the 
# beginning of a script. Also note that a package should NOT be installed every 
# time a script runs.
# The bioconductor repository has installation instructions for biomaRt: 
# https://bioconductor.org/install/

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install(version = "3.14")
}
if (!require("biomaRt", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
# load tidyverse and your new bioconductor package
library(tidyverse)
library(biomaRt)

#### Loading and processing data ####
#' Load Expression Data
#'
#' @param filepath A text string of the full filepath to the file to load.
#'
#' @return A tibble containing the data loaded from the CSV in `filepath`. 
#' 
#' @details Note that not all CSVs are created equal, and there are often cases where 
#' the data will not load in correctly on the first try. You may want to write this functon to 
#' adjust the CSV file being loaded into this assignment so that it can be formed into a 
#' tibble correctly.
#'
#' @examples 
#' `data <- load_expression('/project/bf528/project_1/data/example_intensity_data.csv')`
load_expression <- function(filevariable) {
  expr_mat <-read_delim(filevariable,delim=" ",show_col_types = FALSE)
  return(expr_mat)
}

#' Filter 15% of the gene expression values.
#'
#' @param tibble A tibble of expression values, rows by probe and columns by sample.
#'
#' @return A tibble of affymetrix probe names from the input expression data tibble. 
#' These names match the rows with 15% or more of the expression values about log2(15).
#' 
#' @details This is similar to the filters being implemented in BF528's project 1. 
#' We do not necessarily want to capture all parts of the assay in our analysis, so 
#' filters like this serve to reduce the noise and amount of data to examine.
#'
#' @examples `samples <- filter_15(data_tib)`
#' `> str(samples)`
#' `tibble [40,158 × 1] (S3: tbl_df/tbl/data.frame)`
#' `$ probeids: chr [1:40158] "1007_s_at" "1053_at" "117_at" "121_at" ...`
filter_15 <- function(tibble){
  percent_function <- function(row) {
    filtered_row <- row > log2(15)
    percent_row <- sum(filtered_row)/length(row)
    return(percent_row)
  }

  row_new <- apply(tibble[2:ncol(tibble)], percent_function, MARGIN = 1) 
  filtered_rows <- which(row_new > 0.15) #above the expression value of 15
  return(tibble[filtered_rows, 1])
}

#### Gene name conversion ####

#' Convert affymetrix array names into hgnc_symbol IDs using biomaRt. Inputs and 
#' outputs will likely not be the same size.
#'
#' @param affy_tib A single column tibble of strings containing array names.
#'
#' @return A 2 column tibble that contains affy IDs in the first column,
#' and their corresponding HGNC gene ID in the second column. Note that not all affy IDs 
#' will necessarily correspond with a gene ID, and one gene may have multiple affy IDs.
#' 
#' @details Connecting to ensembl via biomaRt can be...hit or miss...so you may 
#' want to check if data was correctly returned (or if it was just empty). The 
#' `getBM()` function may not accept a tibble, so you might need to convert your 
#' input into a flat vector.
#'
#' @examples 
#' `> affy_to_hgnc(tibble(c('202860_at', '1553551_s_at')))`
#' `affy_hg_u133_plus_2 hgnc_symbol`
#' `1        1553551_s_at      MT-ND1`
#' `2        1553551_s_at       MT-TI`
#' `3        1553551_s_at       MT-TM`
#' `4        1553551_s_at      MT-ND2`
#' `5           202860_at     DENND4B`
affy_to_hgnc <- function(vector_affy) {
  vector_affy <- pull(vector_affy)
  use <- useMart(host = 'https://www.ensembl.org', biomart = 'ENSEMBL_MART_ENSEMBL')
  ensembl <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  newvalues <- getBM(attributes = c('affy_hg_u133_plus_2', 'hgnc_symbol'), 
                    filters = c('affy_hg_u133_plus_2'),
                    values = vector_affy, 
                    mart = ensembl)
  return(as_tibble(newvalues))
}

#### ggplot ####

#' Reduce a tibble of expression data to only the rows in good_genes or bad_genes.
#'
#' @param expr_tibble A tibble holding the expression data, each row corresponding to
#' one affymetrix probe ID and each column to a sample.
#' @param names_ids A two column tibble that associates affy IDs with HGNC gene IDs. 
#' Generated `with affy_to_hgnc()`.
#' @param good_genes A list of gene names stored as a vector of strings.
#' @param bad_genes A list of gene names stored as a vector of strings.
#'
#' @return A tibble with two additional columns added:
#' 1. HGNC gene IDs 
#' 2. Does the gene is this row fall into "good" or "bad" genes?
#' This tibble should be reduced to only rows present in good or bad genes. All
#' other rows can be discarded.
#' 
#' @details In order to plot only our genes of interest, we need to rearrange our 
#' data to include only the elements we want to see. We also want to add to columns, 
#' one that associates the probeids with the HGNC gene name, and one that says if 
#' that gene is in the good or bad sets of genes.
#'
#' @examples 
#' `plot_tibble <- reduce_data(expr_tibble = expr, names_ids = sample_names,`
#' `                           goodGenes, badGenes)`
#' `> head(plot_tibble)`
#' `A tibble: 6 × 38`
#' `  probeids    hgnc    gene_set    GSM972389 ...`
#' `  <chr>       <chr>   <chr>       <dbl>     ...`
#' `1 202860_at   DENND4B good        7.16      ...`
#' `2 204340_at   TMEM187 good        6.40      ...`
reduce_data <- function(expr_tibble, names_ids, good_genes, bad_genes){
  expr_tibble <- add_column(expr_tibble, 
                            .after = "probeids",
                            names_ids[match(expr_tibble$probeids,
                                            names_ids$affy_hg_u133_plus_2), 2])
  good_tib <- expr_tibble[which(expr_tibble$hgnc_symbol %in% good_genes),]
  good_tib <- add_column(good_tib, .after = "hgnc_symbol", gene_set = "good")
  bad_tib <- expr_tibble[which(expr_tibble$hgnc_symbol %in% bad_genes),]
  bad_tib <- add_column(bad_tib, .after = "hgnc_symbol", gene_set = "bad")
  return(rbind(good_tib, bad_tib))
}


#' Plot a boxplot of good and bad genes.
#'
#' @param tibble A reduced tibble of expression data, with information about
#' good and bad genes and gene names.
#'
#' @return A ggplot object which contains a boxplot of the genes and samples we 
#' are interested in.
#' 
#' @details This function performs one additional step before using `ggplot()`: 
#' converting the _wide_ format of the input tibble to a _long_ format.
#'
#' @examples `p <- plot_ggplot(plot_tibble)`
plot_ggplot <- function(tibble) {
  #converting the _wide_ format of the input tibble to a _long_ format.
  tibble<- tibble %>% pivot_longer(starts_with('GSM'), names_to = 'SAMPLES',names_sep = NULL,names_pattern = NULL,names_transform = NULL) 
  #extracting the hgyc symbol col in tibble
  x <- tibble$hgnc_symbol
  x <- factor(x, levels = unique(x,ordered = TRUE))
  m <- ggplot(tibble, aes(x=hgnc_symbol, y=value)) +geom_boxplot(aes(fill = gene_set)) +theme_gray(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),legend.position = "right") +
    scale_fill_manual(values = c("#4DB3E6", "#AA4371"),aesthetics = "colour") +
    ggtitle("Boxplot of randomly chosen genes good or bad") +
    xlab("Genes") +
    ylab("Expression levels of genes")
  return(m)
  
}

