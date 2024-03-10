#!/usr/bin/env Rscript

#######################################
## scASfind detect all marker nodes
## ysong 8 Mar 2024
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
  make_option(c("-m", "--num_markers"),
    type = "numeric", default = 2000,
    help = "How mant top marker nodes to return"
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
num_markers <- opt$num_markers

obj <- loadObject(index_path)

all_ct_above = cellTypeNames(obj, datasets = 'above')

mk_above = data.frame()

for (ct_now in all_ct_above) {

add = cellTypeMarkers(object = obj, cell.types = ct_now, background.cell.types = cellTypeNames(obj, datasets = 'above'), 
                top.k = num_markers , sort.field = 'f1') %>% arrange(desc(f1))
mk_above = rbind(mk_above, add)    
    
    }

all_ct_below = cellTypeNames(obj, datasets = 'below')

mk_below = data.frame()

for (ct_now in all_ct_below) {

add = cellTypeMarkers(object = obj, cell.types = ct_now, background.cell.types = cellTypeNames(obj, datasets = 'below'), 
                top.k = num_markers , sort.field = 'f1') %>% arrange(desc(f1))
mk_below = rbind(mk_below, add)    
    
    }

rbind(mk_above, mk_below) %>% mutate(mk_type = gsub("\\..*", "", cellType)) %>% 
mutate(cell_type = gsub("above\\.|below\\.", "", cellType)) %>% 
write_csv(paste0(output, "/", name, "_top_", num_markers, "_markers.csv"))
