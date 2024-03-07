#!/usr/bin/env Rscript

#######################################
## scASfind detect coordinated spliced-in events
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


index <- loadObject(index_path)

all_nodes = scfindNodes(index)

all_genes <- levels(factor(nodeDetails(index, node.list = all_nodes)$Gene_name))

all_blocks <- data.frame()

node_types <- str_split(types, ",") %>% flatten_chr()

block_num <- 1

message("Start looking for coordinated splice-in events")

for (gene in all_genes) {

  # get all encoded nodes per gene
  nodes <- geneNodes(index, gene, "Gene_name")
    
  nodes <- nodes %>% filter(Node_id %in% all_nodes)

  # how many nodes are in between consecutive records
  tbl <- index@metadata$stats[which(rownames(index@metadata$stats) %in% nodes$Node_id), ] %>%
    tibble::rownames_to_column("Node_id") %>%
    mutate(node_num = as.numeric(gsub("^.*_", "", Node_id))) %>%
    arrange(node_num) %>%
    mutate(node_num_diff = ave(node_num, FUN = function(x) c(0, diff(x))))

  details <- nodeDetails(index, tbl$Node_id)

  tbl <- merge(tbl, details, by = "Node_id") %>%
    arrange(node_num) %>%
    filter(Type %in% node_types)

  new_block <- data.frame()

  block_num <- 1

  # if a gene has more than 2 nodes, enables it to test for node blocks
  if (nrow(tbl) > 2) {

    # first node in gene
    i <- 1
    avg <- tbl[i, "mean"]
    SD <- tbl[i, "SD"]

    # start a new node block record
    new_block <- data.frame()

    while (i < nrow(tbl)) {

      # extend the node block by the next nodes

      # first test for block pattern using mean and SD

      if (abs(tbl[i + 1, "mean"] - avg) < mean_cutoff &
        abs(tbl[i + 1, "SD"] - SD < SD_cutoff)) {

        # add i+1 node to the block
        add_block <- tbl[i + 1, ]

        add_block$block_num <- block_num

        # update node block average PSI and SD
        avg <- (avg + add_block$mean) * 0.5
        SD <- sqrt(SD^2 + add_block$SD^2)

        if (nrow(new_block) == 0) {
          # the node is the first to add to the block
          first_in_block <- tbl[i, ]
          first_in_block$block_num <- block_num
          new_block <- rbind(first_in_block, add_block)
        } else {
          # the node is to extend an existing block
          new_block <- rbind(new_block, add_block)
        }
        # proceed to next node in gene, keep adding nodes to this block
        i <- i + 1
      } else {

        # the next node does not pass the pattern filter, nothing more to add for this block

        if (nrow(new_block) == 0) {

          # nothing passed the filter: go to next node in this gene
          i <- i + 1
          avg <- tbl[i, "mean"]
          SD <- tbl[i, "SD"]
        } else if (nrow(new_block) > 0) {

          # potential new block to be added
          # check cell type specificity

          block_now <- block_num

          test_comb <- new_block %>%
            dplyr::filter(block_num == block_now) %>%
            select(Node_id)

          skip_to_next <- FALSE

          condition <- tryCatch(
            {
              sum(scASfind::hyperQueryCellTypes(index, test_comb$Node_id)$pval < pval_cutoff) >= 1
            },
            error = function(e) {
              skip_to_next <<- TRUE
            }
          )

          # if the block is significant in some cell types, add it to the results

          if (condition == TRUE) {
            message(paste("add block No. ", block_num, " of gene: ", gene, sep = ""))

            sig_cell_types <- scASfind::hyperQueryCellTypes(index, test_comb$Node_id) %>% filter(pval < 0.05)

            add_new_block <- merge(new_block, sig_cell_types, all = TRUE)

            all_blocks <- rbind(all_blocks, add_new_block)

            # start a new block in this gene

            block_num <- block_num + 1

            new_block <- data.frame()

            i <- i + 1

            avg <- tbl[i, "mean"]

            SD <- tbl[i, "SD"]
          } else {
            # the block is not specific in any cell types
            # proceed to the next node

            i <- i + 1

            avg <- tbl[i, "mean"]

            SD <- tbl[i, "SD"]

            # no need to update block num because no new block is added in this gene

            new_block <- data.frame()

            next
          }
        }
      }
    }
  }
}

all_blocks <- unique(all_blocks)

message("Finish detecting coordinated spliced-in events")

saveRDS(all_blocks, paste(output, "/", name, "_all_node_blocks.rds", sep = ""))
