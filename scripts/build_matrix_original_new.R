#!/usr/bin/env Rscript

#######################################
## scASfind build a splicing index
## ysong 14 Apr 2022
#######################################

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(biomaRt))

option_list <- list(
  make_option(c("-n", "--data_name"), type = "character", default = NULL, help = "Name of dataset"),
  make_option(c("-p", "--pseudobulk_psi"),
    type = "character", default = NULL, help = "Combined pseudobulk psi matrix as input, usually with .psi.tsv extensions, tab-deliminated"
  ),
  make_option(c("-i", "--node_info"),
    type = "character", default = NULL, help = "Node information input from Whippet, tab-deliminated"
  ),
  make_option(c("-o", "--output"),
    type = "character",
    default = NULL, help = "Directory where multiple output RDS file will be written"
  ),
  make_option(c("-r", "--num_reads_min"),
    type = "numeric", default = 10,
    help = "Minimum number of total reads covering node, which will be included in the output. This is for ensure meaningful psi quantification, default = 10"
  ),
  make_option(c("-d", "--psi_diff_cutoff"), type = "numeric", default = 0.2, help = "Minimum PSI difference from dataset average which will lead to the node being kept in the index, default = 0.2"),
  make_option(c("-s", "--species"),
    type = "character", default = NULL, help = "ENSEMBL name of species to get gene annotations"
  )
)

# parse input
opt <- parse_args(OptionParser(option_list = option_list))

NAME <- opt$data_name
INPUT <- opt$pseudobulk_psi
NODE <- opt$node_info
OUTPUT <- opt$output
NUM_READS_MIN <- opt$num_reads_min
PSI_DIFF_CUTOFF <- opt$psi_diff_cutoff
SPECIES <- opt$species

######################## Functions build various inputs for scASfind

# build original PSI matrix from MicroExonator output file *.psi.tsv

message("Start building inputs for scASfind index")

data <- readr::read_tsv(INPUT, col_names = TRUE, progress = show_progress())

message("Read PSI input, generating matrces...")

matrix.original <- data %>%
  tidyr::unite("Gene_node", Gene, Node, sep = "_") %>%
  dplyr::group_by(Sample) %>%
  dplyr::distinct(Gene_node, .keep_all = TRUE) %>%
  dplyr::select(Sample, Gene_node, Total_Reads, Psi) %>%
  dplyr::mutate(Filter = case_when(Total_Reads < NUM_READS_MIN ~ "DROP", Total_Reads >= NUM_READS_MIN ~ "KEEP", TRUE ~ NA_character_)) %>% ## keep only those PSI qualifications by >= NUM_READS_MIN as they have an acceptable confidence interval
  dplyr::filter(Filter == "KEEP") %>%
  dplyr::select(Sample, Gene_node, Psi) %>%
  dplyr::ungroup(Sample) %>%
  dplyr::group_by(Gene_node) %>%
  tidyr::pivot_wider(names_from = Sample, values_from = Psi)

message("Matrices constructed, performing scailing")

df <- data.frame(matrix.original, row.names = matrix.original$Gene_node)
dm <- as.matrix(df[, -1])
mean <- rowMeans(dm, na.rm = TRUE)

# get differential from dataset mean PSI per node

matrix.scaled_diff <- dm - mean

# times a scale factor 100 for accurate encoding of absolute PSI values
matrix.scaled_diff <- matrix.scaled_diff * 100

# drop all-na rows

matrix.scaled_diff <- matrix.scaled_diff[which(rowSums(is.na(matrix.scaled_diff)) < ncol(matrix.scaled_diff)), ]
matrix.scaled_diff_selected <- data.frame(row.names = rownames(matrix.scaled_diff))

# initialize diff_cut matrix to store index of nodes with PSI = dataset mean, to distinguish from nodes with PSI not quantified
diff_cut <- matrix(0, nrow = nrow(matrix.scaled_diff), ncol = ncol(matrix.scaled_diff))

# keep only values with above psi_diff_cutoff

for (cell in seq(1, ncol(matrix.scaled_diff))) {
  tv <- matrix.scaled_diff[, cell]
  tvd <- diff_cut[, cell]
  temp <- which(abs(tv) < psi_diff_cutoff * 100)
  tv[temp] <- NA
  tvd[temp] <- 1
  matrix.scaled_diff_selected[[colnames(matrix.scaled_diff)[cell]]] <- tv
  diff_cut[, cell] <- tvd
}

# 1 in diff_cut means dropped calls due to equal to mean

rownames(diff_cut) <- rownames(matrix.scaled_diff_selected)
colnames(diff_cut) <- colnames(matrix.scaled_diff_selected)

# remove all na rows
matrix.scaled_diff_selected <- matrix.scaled_diff_selected %>% filter(if_any(everything(), ~ !is.na(.)))

diff_cut <- diff_cut[rownames(matrix.scaled_diff_selected), ]
diff_cut <- Matrix(diff_cut, sparse = TRUE)

matrix.above <- data.frame(row.names = rownames(matrix.scaled_diff_selected))

# set PSI of NA and below ones to zero for the above index

for (cell in seq(1, ncol(matrix.scaled_diff_selected))) {
  tv <- matrix.scaled_diff_selected[, cell]
  temp <- which(is.na(tv) | tv < 0)
  tv[temp] <- 0
  matrix.above[[colnames(matrix.scaled_diff_selected)[cell]]] <- tv
}

matrix.below <- data.frame(row.names = rownames(matrix.scaled_diff_selected))

# set NA and above ones to zero for the below indes

for (cell in seq(1, ncol(matrix.scaled_diff_selected))) {
  tv <- matrix.scaled_diff_selected[, cell]
  temp <- which(is.na(tv) | tv > 0)
  tv[temp] <- 0
  matrix.below[[colnames(matrix.scaled_diff_selected)[cell]]] <- tv
}
matrix.below <- matrix.below * (-1)

# get mean and SD per node

df <- data.frame(matrix.original, row.names = matrix.original$Gene_node)
dm <- as.matrix(df[, -1])
mean <- transform(dm, mean = apply(dm, 1, mean, na.rm = TRUE))
sd <- transform(dm, SD = apply(dm, 1, sd, na.rm = TRUE))
mean$SD <- sd$SD
mean <- mean[order(mean$SD), ]
stats <- mean[, c("mean", "SD")]

stats <- stats[which(rownames(stats) %in% rownames(matrix.scaled_diff_selected)), ]

# get node annotations from ENSEMBL

node_list <- rownames(matrix.scaled_diff_selected)

message("Read node information...")

ni <- readr::read_tsv(NODE)
ni$Node_id <- paste(ni$Gene, ni$Node, sep = "_")

node_list_all <- ni[which(ni$Node_id %in% node_list), ]
node_list_all$Gene_num <- gsub("\\..*$", "", node_list_all$Gene)

# install.packages('XML', repos = 'http://www.omegahat.net/R') BiocManager::install('biomaRt')
library("biomaRt")

ensembl <- useMart("ensembl")
ensembl <- useDataset(paste0(SPECIES, "_gene_ensembl"), mart = ensembl)

# takes a few minutes to match gene to name
gene_name <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters = "ensembl_gene_id", values = node_list_all$Gene_num, mart = ensembl)

# de-duplicate and match gene name using gene id
gene_name <- gene_name[!duplicated(gene_name$ensembl_gene_id), ]

gene_node_all <- merge(node_list_all, gene_name, by.x = "Gene_num", by.y = "ensembl_gene_id", all.x = TRUE)

names(gene_node_all)[names(gene_node_all) == "Gene"] <- "Gene_id"
names(gene_node_all)[names(gene_node_all) == "external_gene_name"] <- "Gene_name"
names(gene_node_all)[names(gene_node_all) == "Gene_node"] <- "Node_id"
gene_node_all$Node_name <- paste(gene_node_all$Gene_name, gene_node_all$Node, sep = "_")

message("Save results...")

saveRDS(matrix.scaled_diff_selected, paste(OUTPUT, "_matrix_scaled_diff_selected.rds", sep = ""))
saveRDS(matrix.above, paste(OUTPUT, "_matrix_above.rds", sep = ""))
saveRDS(matrix.below, paste(OUTPUT, "_matrix_below.rds", sep = ""))
saveRDS(diff_cut, paste(OUTPUT, "_diff_cut.rds", sep = ""))
saveRDS(stats, paste(OUTPUT, "_stats.rds", sep = ""))
saveRDS(gene_node_all, paste(OUTPUT, "_gene_node_all.rds", sep = ""))
