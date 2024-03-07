#!/usr/bin/env Rscript

#######################################
## scASfind scan for mutually exclusive exons
## ysong 14 Apr 2022
#######################################

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scASfind))

option_list <- list(
  make_option(c("-n", "--data_name"),
    type = "character", default = NULL,
    help = "Name of dataset"
  ),
  make_option(c("-i", "--index"),
    type = "character", default = NULL,
    help = "Path to a complete scASfind index"
  ),
  make_option(c("-t", "--node_types", default = "CE,AD,AA,NA,RI"),
    type = "character", default = NULL,
    help = "Node types to consider when detecting blocks, split by comma"
  ),
  make_option(c("-m", "--mean_cutoff"),
    type = "numeric", default = 0.1,
    help = "Potential MXE Nodes are first filtered by having a difference of dataset-wise mean within the cutoff, default 0.1"
  ),
  make_option(c("-s", "--SD_cutoff"),
    type = "numeric", default = 0.1,
    help = "Potential MXE Nodes are first filtered by having a difference of dataset-wise SD within the cutoff, default 0.1"
  ),
  make_option(c("-p", "--pval_cutoff"),
    type = "numeric", default = 0.05,
    help = "Cut-off of P value of hypergeometric test for cell type significance of MXE pairs"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = NULL,
    help = "Directory where discovered node blocks in all genes will be written as an .rds file"
  )
)

# parse input
opt <- parse_args(OptionParser(option_list = option_list))

name <- opt$data_name
index_path <- opt$index
output <- opt$output
types <- opt$node_types
mean_cutoff <- opt$mean_cutoff
SD_cutoff <- opt$SD_cutoff
pval_cutoff <- opt$pval_cutoff

index <- scASfind::loadObject(index_path)

stats <- index@metadata$stats

stats$Node_id <- rownames(stats)

node.list <- index@metadata$node_list

types <- str_split(types, ",") %>% flatten_chr()

a <- merge(stats, node.list,
  by = "Node_id",
  all.x = FALSE,
  all.y = FALSE
) %>% unique()


d <- data.frame()

message("start processing all genes to find MXEs")

node_num <- 1

for (gene in levels(factor(a$Gene_name))[nzchar(levels(factor(a$Gene_name)))]) {
  message(gene)
  nodes <- geneNodes(index, gene, "Gene_name")
  nodes <- nodes %>% filter(Node_name %in% a$Node_name)
  nodes$Type <- as.character(nodes$Type)

  ## a gene does not have nodes
  if (nrow(nodes) == 0) {
    next
  }

  ## the gene has more than one node of specified types to check
  if (nrow(nodes) > 1) {
    nodes_check <- nodes[which(nodes$Type %in% types), "Node_id"]
    message(nodes_check)

    ## test all possible pairs of nodes of specified types
    if (length(nodes_check) >= 2) {
      pairs <- as.data.frame(combn(nodes_check, 2))

      for (i in seq(1, ncol(pairs))) {
        skip_to_next <- FALSE
        node_1 <- pairs[1, i]
        node_2 <- pairs[2, i]

        a_1 <- a[which(a$Node_id %in% node_1), ]

        a_2 <- a[which(a$Node_id %in% node_2), ]


        if (nrow(a_1 > 0) & nrow(a_2) > 0) {
          # two mutually exclusive patterns

          test_comb <- c(node_1, paste("-", node_2, sep = ""))
          message(test_comb)

          test_comb_2 <- c(node_2, paste("-", node_1, sep = ""))

          # first filter by mean and SD

          if ((1 - mean_cutoff) <= (a_1[1, "mean"] + a_2[1, "mean"]) &
            (a_1[1, "mean"] + a_2[1, "mean"]) <= (1 + mean_cutoff) &
            abs(a_1[1, "SD"] - a_2[1, "SD"]) <= SD_cutoff) {
            # then test cell type specificity: at least in one cell type, one pattern is specific

            condition <- tryCatch(
              {
                suppressMessages(sum(scASfind::hyperQueryCellTypes(index, test_comb, datasets = "above")$pval < pval_cutoff) >= 1 | sum(scASfind::hyperQueryCellTypes(index, test_comb_2, datasets = "above")$pval < pval_cutoff) >= 1)
              },
              error = function(e) {
                skip_to_next <<- TRUE
              }
            )

            if (skip_to_next) {
              next
            } else if (condition == TRUE) {
              # this is a potential mutually exclusive exon

              message("find a mutually exclusive exon pair that is cell type specific")

              if (sum(scASfind::hyperQueryCellTypes(index, test_comb, datasets = "above")$pval < pval_cutoff) >= 1 & sum(scASfind::hyperQueryCellTypes(index, test_comb_2, datasets = "above")$pval < pval_cutoff) >= 1) {
                sig_cell_types <- suppressMessages(rbind(scASfind::hyperQueryCellTypes(index, test_comb, datasets = "above") %>% filter(pval < pval_cutoff), scASfind::hyperQueryCellTypes(index, test_comb_2, datasets = "above") %>% filter(pval < pval_cutoff)))
              } else if (sum(scASfind::hyperQueryCellTypes(index, test_comb, datasets = "above")$pval < pval_cutoff) >= 1) {
                sig_cell_types <- suppressMessages(scASfind::hyperQueryCellTypes(index, test_comb, datasets = "above") %>% filter(pval < pval_cutoff))
              } else if (sum(scASfind::hyperQueryCellTypes(index, test_comb_2, datasets = "above")$pval < pval_cutoff) >= 1) {
                sig_cell_types <- suppressMessages(scASfind::hyperQueryCellTypes(index, test_comb_2, datasets = "above") %>% filter(pval < pval_cutoff))
              }

              add <- rbind(a_1, a_2)
              add <- merge(add, sig_cell_types)
              add$node_num <- node_num

              # count the occurance of MXE pairs in this gene
              node_num <- node_num + 1
              d <- rbind(d, add)
            }
          }
        }
      }
    }
  }
}

message("finish screening for all MXEs")
saveRDS(d, paste(output, "/", name, "_mutually_exclusive_exons.rds", sep = ""))
