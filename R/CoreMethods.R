#' Builds an \code{SCFind} object from a \code{matrix}
#'
#' This function will index a \code{matrix} including alternative splicing information for each cell as an SCFind index.
#'
#' @param psival a dataframe object containing psi value of each node of a gene and cell with gene_node as row.names and cell id as col.names
#' @param dataset.name name of the dataset, i.e. tissue name
#' @param metadata metadata of the alternative splicing dataset, must include cell.type information for the dataset
#' @param column.label the cell.type in metadata that will be used for the index
#' @param qb number of bits per cell that are going to be used for quantile compression of the expression data
#'
#' @name buildAltSpliceIndex
#'
#' @return an SCFind object
#'
#' @importFrom hash hash
#' @importFrom methods new
#'
#' @importFrom Rcpp cpp_object_initializer
#' @useDynLib scASfind
#'

buildAltSpliceIndex.PSI <- function(psival, metadata, dataset.name, column.label, qb = 2) {
  if (missing(dataset.name)) {
    stop("Please name your dataset with dataset.name")
  }
  if (grepl(dataset.name, ".")) {
    stop("The dataset name should not contain any dots")
  }

  cell.types.all <- as.factor(metadata[, column.label])
  cell.types <- levels(cell.types.all)
  new.cell.types <- hash(keys = cell.types, values = paste0(dataset.name, ".", cell.types))
  node.names <- unique(rownames(psival))

  if (length(cell.types) > 0) {
    non.zero.cell.types <- c()
    index <- hash()
    message(paste("Generating index for", dataset.name))

    ef <- new(EliasFanoDB)
    qb.set <- ef$setQB(qb)
    if (qb.set == 1) {
      stop("Setting the quantization bits failed")
    }

    for (cell.type in cell.types) {
      inds.cell <- which(cell.type == cell.types.all) # which cells belong to this cell type
      if (length(inds.cell) < 1) {
        message(paste("Skipping", cell.type))
        next
      }
      non.zero.cell.types <- c(non.zero.cell.types, cell.type)
      message(paste("\tIndexing", cell.type, "as", new.cell.types[[cell.type]], " with ", length(inds.cell), " cells."))

      # now build index
      ## order cells in the psival matrix as you order the cellls in metadata
      cell.type.psi.scaled <- psival[, inds.cell]

      if (is.matrix(psival)) {
        ef$indexMatrix(new.cell.types[[cell.type]], as.matrix(cell.type.psi.scaled))
      } else {
        ef$indexMatrix(new.cell.types[[cell.type]], as.matrix(cell.type.psi.scaled))
      }
    }
  }
  index <- new("SCFind", index = ef, datasets = dataset.name)
  return(index)
}



#' @rdname buildAltSpliceIndex
#' @aliases buildAltSpliceIndex
setMethod("buildAltSpliceIndex",
  definition = buildAltSpliceIndex.PSI
)


#' Add necessary elements of the metadata slot of an altervnative splicing SCFind index object
#'
#' @param object an SCFind class object built by the method "builtAltSpliceIndex"
#' @param stats the statistics matrix built by the function "buildMatrix.stats"
#' @param node_list the node information matrix built by the function "buildMatrix.node_list"
#' @param diff_cut sparse matrix with 1 denoting the removal of node PSI due to not sufficiently deviate from tissue mean
#' @name addIndexMeta
#'
#' @return an SCFind object with classic metadata components
#' @useDynLib scASfind
#'

addIndexMeta.classic <- function(object, stats, node_list, diff_cut) {
  object@metadata$stats <- stats
  object@metadata$node_list <- node_list
  object@metadata$diff_cut <- diff_cut
  return(object)
}

#' @rdname addIndexMeta
#' @aliases addIndexMeta
setMethod("addIndexMeta",
  definition = addIndexMeta.classic
)



#' Runs a query and performs the hypergeometric test for the retrieved cell types
#'
#' @name hyperQueryCellTypes
#' @param object the \code{SCFind} object
#' @param node.list AS nodes to be searched in the node.list.index
#' (Operators: "-gene" to exclude a gene | "*gene" either gene is expressed
#' "*-gene" either gene is expressed to be excluded)
#' @param datasets the datasets vector that will be tested as background for the hypergeometric test
#'
#' @return a DataFrame that contains all cell types with the respective cell cardinality and the hypergeometric test

cell.types.phyper.test.AS <- function(object, node.list, datasets) {
  continue <- FALSE
  node.list.2 <- gsub("[\\*\\-]", "", node.list)

  if (is.null(object@metadata$node_list)) {
    question1 <- readline("Warning: missing node_list metadata in index, can not verify existance of query nodes in index. This may cause the session to terminate. \nAre you sure the node is in index and would like to continue query? (Y/N)")

    if (regexpr(question1, "y", ignore.case = TRUE) == 1) {
      continue <- TRUE
    } else if (regexpr(question1, "n", ignore.case = TRUE) == 1) {
      return("Exit query")
    }
  } else {
    if (!all(node.list.2 %in% scfindNodes(object))) {
      stop("Query nodes not in index, please change your query")
    } else {
      # message("Verified all query nodes are in index, generating results...")
      continue <- TRUE
    }
  }
  if (continue == TRUE) {
    result <- findCellTypes.geneList(object, node.list, datasets)
    if (!identical(result, list())) {
      return(phyper.test(object, result, datasets))
    } else {
      # message("No Cell Is Found!")
      return(data.frame(cell_type = c(), cell_hits = c(), total_cells = c(), pval = c()))
    }
  }
}

#' @rdname hyperQueryCellTypes
#' @aliases hyperQueryCellTypes
setMethod("hyperQueryCellTypes",
  signature(
    object = "SCFind",
    node.list = "character"
  ),
  definition = cell.types.phyper.test.AS
)



#' This function retrieves node details by node id
#'
#' @name nodeDetails
#' @param object the \code{SCFind} object
#' @param node.list AS node ids to find their details
#' @return a dataframe that contains details of the query nodes

node.details <- function(object, node.list) {
  if (is.null(object@metadata$node_list)) stop("Missing node details in index metadata")

  details <- object@metadata$node_list[which(as.character(object@metadata$node_list[["Node_id"]]) %in% node.list), ]

  details_ordered <- details[match(node.list, details$Node_id), ]

  return(details_ordered)
}

#' @rdname nodeDetails
#' @aliases nodeDetails
setMethod("nodeDetails",
  signature(
    object = "SCFind",
    node.list = "character"
  ),
  definition = node.details
)


#' This function retrieves node details in an index by gene id or gene name
#'
#' @name geneNodes
#' @param object the \code{SCFind} object
#' @param gene.list gene id or gene name list to find nodes
#' @param query.type either "Gene", "external_gene_name", "Gene_node" to use in query
#' @return a dataframe that contains nodes for gene.list

gene.nodes <- function(object, gene.list, query.type) {
  if (is.null(object@metadata$node_list)) stop("Missing node details in index metadata")
  if (!query.type %in% c("Gene_id", "Gene_name", "Node_id", "Node_name")) stop("query.type must be \"Gene_id\", \"Gene_name\", \"Node_id\" or \"Node_name\"")
  node.list <- subset(object@metadata$node_list, as.character(object@metadata$node_list[[query.type]]) %in% gene.list, )
  if (nrow(node.list) == 0) {
    warning("No node is found in this index, please change your query")
    return(data.frame())
  }
  return(node.list)
}

#' @rdname geneNodes
#' @aliases geneNodes
setMethod("geneNodes",
  signature(
    object = "SCFind",
    gene.list = "character",
    query.type = "character"
  ),
  definition = gene.nodes
)

#' This function finds coordinated node sets for a gene
#'
#' @name findNodeSets
#' @param object the \code{SCFind} object
#' @param gene.list gene id or gene name to find coordinated node sets
#' @param query.type either "gene_id" or "external_gene_name" to use in query
#' @param node.types types of splicing nodes to consider in the gene.list
#' @return a dataframe that contains nodes for gene.list

gene.node.sets <- function(object, gene.list, query.type, node.types = c("CE", "AA", "AD", "RI", NA, "NA")) {
  nodes <- gene.nodes(object, gene.list, query.type)
  nodes.new <- nodes[which(as.character(nodes$Type) %in% node.types), ]

  if (nrow(nodes.new) != 0) {
    nodes <- nodes.new
  } else {
    warning(paste("there is no ", node.types, " node in this gene, please change query", sep = ""))
    return(NA)
  }

  markers <- markerNodes(object, as.character(nodes.new$Node_id))
  if (nrow(markers) == 0) stop("No gene pattern is found")
  sets <- data.frame()

  query <- strsplit(as.character(markers[which.max(markers$tfidf), "Query"]), ",")[[1]]
  result <- cell.types.phyper.test.AS(object, query)

  message("running hyperQueryCellType using")
  print(query)
  print(result)


  for (j in seq(1, nrow(result))) {
    if (result$pval[[j]] <= 0.05) {
      message("find a cell type specific node set \n")
      print(query)
      sets <- rbind(sets, result[j, ])
    }
  }

  if (nrow(sets) == 0) {
    message("This query is not cell type specific")
  }

  return(sets)
}


#' @rdname findNodeSets
#' @aliases findNodeSets
setMethod("findNodeSets",
  signature(
    object = "SCFind",
    gene.list = "character",
    query.type = "character",
    node.types = "character"
  ),
  definition = gene.node.sets
)


#' This function finds coordinated node sets for a gene
#'
#' @name getCoordinatedNodes
#' @param object the \code{SCFind} object
#' @param gene.name Gene_name of the interested gene from geneNodes
#' @return a dataframe that contains top queries
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange slice_head filter slice_max
#'
get_coordinated_nodes <- function(object, gene.name) {
  query <- markerNodes(index, geneNodes(index, gene.name, query.type = "Gene_name")$Node_id) %>%
    arrange(desc(tfidf)) %>%
    slice_head(n = 30) %>%
    filter(Cells == Mode(Cells)) %>%
    slice_max(Genes, n = 5)

  return(query)
}

#' @rdname getCoordinatedNodes
#' @aliases getCoordinatedNodes
setMethod("getCoordinatedNodes",
  signature(
    object = "SCFind",
    gene.name = "character"
  ),
  definition = get_coordinated_nodes
)



#' This function gets the raw PSI of a node in a cell type
#'
#' @name getRawPsi
#' @param object the \code{SCFind} object
#' @param node.list several nodes that we wish to get the raw PSI
#' @param cell.types cell type tp query, can only query one cell type at once
#' @return a dataframe that contains raw psi value in the queried cell type of the gene.list
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#'

get.raw.psi <- function(object, node.list, cell.types) {
  cell_types_all <- cell.types

  gene_nodes_all <- node.list

  raw_psi <- data.frame()


  for (cell_type in cell_types_all) {
    message(cell_type)

    raw_psi_ct <- get.cell.type.raw.psi(object, gene_nodes_all, cell_type)

    colnames(raw_psi_ct) <- c(paste(cell_type, seq(1, ncol(raw_psi_ct)), sep = "_"))
    raw_psi_add <- raw_psi_ct %>%
      rownames_to_column("Node_id") %>%
      pivot_longer(-c(Node_id), names_to = "Cell_id", values_to = "raw_psi")
    if (nrow(raw_psi) == 0) {
      raw_psi <- raw_psi_add
    } else {
      raw_psi <- rbind(raw_psi, raw_psi_add)
    }

    message(paste("Added", cell_type, sep = " "))
  }

  details <- nodeDetails(object, gene_nodes_all)
  raw_psi_full <- merge(raw_psi, details, by = "Node_id")
  return(raw_psi_full)
}



#' @rdname getRawPsi
#' @aliases getRawPsi
setMethod("getRawPsi",
  signature(
    object = "SCFind",
    node.list = "character",
    cell.types = "character"
  ),
  definition = get.raw.psi
)


#' This function plots correlation of splicing PSI and return significantly correlated nodes
#'
#' @name plotRawPsiCorr
#' @param raw_psi raw psi matrx, each row is a gene, each col is a cell
#' @param node.list node list whose PSI is to be plotted, default 'all_nodes'
#' @param cell.types cell types to consider, default 'all_cell_types'
#' @importFrom Hmisc rcorr
#' @import ggplot2
#' @return a list object with heatmap, pos corr and neg corr
#'

plot.raw.psi.corr <- function(raw_psi, node.list = "all_nodes", cell.types = "all_cell_types") {

  # in raw_psi, each row is a node, each col is a cell type
  # we need cor among nodes, so transform

  if (node.list == "all_nodes") node.list <- rownames(raw_psi)

  if (cell.types == "all_cell_types") cell.types <- colnames(raw_psi)

  raw_psi <- raw_psi[which(rownames(raw_psi) %in% node.list), which(colnames(raw_psi) %in% cell.types)]


  if (nrow(raw_psi) == 0 | ncol(raw_psi) == 0) {
    warning("No value to plot, please change query")
    return(NA)
  }

  corr <- cor(t(raw_psi), method = "pearson", use = "pairwise.complete.obs")


  corr[which(is.na(corr))] <- 0


  # Reorder the correlation matrix
  cormat <- reorder_cormat(corr)
  upper_tri <- get_upper_tri(cormat)

  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)

  labels <- melted_cormat$Var1

  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white",
      midpoint = 0, limit = c(-1, 1), space = "Lab",
      name = "Pearson\nCorrelation"
    ) +
    theme_minimal() + # minimal theme
    theme(
      axis.title = element_text(size = 12),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 0.5),
      axis.text.y = element_text(size = 10)
    ) +
    coord_fixed() +
    scale_x_discrete(
      breaks = melted_cormat$Var1,
      labels = as.character(melted_cormat$Var1)
    ) +
    labs(title = "PSI correlation of query nodes using query cell types", x = "Gene_node", y = "Gene_node")

  # Print the heatmap
  # Takes a while
  print(ggheatmap)

  # get pos and neg correlated nodes
  melted_cormat$Var1 <- as.character(melted_cormat$Var1)
  melted_cormat$Var2 <- as.character(melted_cormat$Var2)


  pos_corr <- melted_cormat[which(melted_cormat$value > 0 & melted_cormat$value < 1), ]

  neg_corr <- melted_cormat[which(melted_cormat$value < 0), ]

  # test the significance of correlation
  # library(Hmisc)
  res <- suppressWarnings(Hmisc::rcorr(as.matrix(t(raw_psi))))

  round_p <- round(res$P, 3)

  round_p <- get_upper_tri(round_p)

  melted_p <- melt(round_p, na.rm = TRUE)

  sig <- melted_p[which(melted_p$value < 0.05), ]

  # print(levels(sig$Var1))

  pos_sig <- pos_corr[which(pos_corr$Var1 %in% levels(sig$Var1)), ]
  message("significantly positively correlated nodes :")
  print(pos_sig)

  neg_sig <- neg_corr[which(neg_corr$Var1 %in% levels(sig$Var1)), ]

  message("significantly negatively correlated nodes :")
  print(neg_sig)

  final <- list("heatmap" = ggheatmap, "sig_pos_corr" = pos_sig, "sig_neg_corr" = neg_sig)
  return(final)
}



#' @rdname plotRawPsiCorr
#' @aliases plotRawPsiCorr
setMethod("plotRawPsiCorr",
  signature(
    raw_psi = "data.frame",
    node.list = "character",
    cell.types = "character"
  ),
  definition = plot.raw.psi.corr
)





#' This function plots heatmap of PSI values
#'
#' @name plotRawPsiHeatmap
#' @param raw_psi the raw_psi matrix from getRawPsi
#' @param node.list a list of nodes whose raw psi is to be plotted
#' @param cell.types a list of cell types whose raw psi is to be plotted
#' @return a heatmap ggplot object
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange mutate
#' @importFrom rquery natural_join
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom forcats fct_inorder

plot.raw.psi.heatmap <- function(raw_psi, node.list, cell.types) {
  ggheatmap <- raw_psi %>%
    mutate(node_num = as.numeric(Node)) %>%
    arrange(node_num) %>%
    mutate(Node_name = factor(Node_name)) %>%
    mutate(Cell_type = gsub("_\\d*$", "", Cell_id)) %>%
    group_by(Cell_type, Node_id) %>%
    mutate(mean_psi = mean(raw_psi, na.rm = TRUE)) %>%
    ungroup() %>%
    ggplot(aes(x = fct_inorder(Node_name), y = Cell_type, fill = mean_psi)) +
    geom_tile() +
    theme_minimal() +
    labs(
      title = "Raw PSI value for query nodes",
      x = "node_name",
      y = "cell_type"
    ) +
    theme(axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 0.5
    )) +
    scale_fill_gradient(
      low = "#D01B1B", high = "#47abd8",
      limits = c(0, 1), na.value = "white"
    )
  return(ggheatmap)
}

#' @rdname plotRawPsiHeatmap
#' @aliases plotRawPsiHeatmap
setMethod("plotRawPsiHeatmap",
  signature(
    raw_psi = "data.frame",
    node.list = "character",
    cell.types = "character"
  ),
  definition = plot.raw.psi.heatmap
)





#' This function serializes the DB and save the object as an rds file
#'
#' This function can be used to enable the user save the loaded file in a database
#' to avoid re-indexing and re-merging individual assays.
#'
#' After serializing and saving it clears the redundant bytestream from memory
#' because the memory is already loaded in memory
#' @param object an SCFind object
#' @param file the target filename that the object is going to be stored
#'
#' @return the \code{SCFind} object
#' @name saveObject
save.serialized.object <- function(object, file) {
  object@serialized <- object@index$getByteStream()
  a <- saveRDS(object, file)
  # Clear the serialized stream
  object@serialized <- raw()
  gc()
  return(object)
}

#' @rdname saveObject
#' @aliases saveObject
setMethod("saveObject", definition = save.serialized.object)


#' This function loads a saved \code{SCFind} object and deserializes
#' the object and loads it into an in-memory database.
#'
#' After loading the database it clears the loaded bytestream from the memory.
#'
#' @param filename the filepath of a specialized serialized scfind object
#'
#' @return an \code{SCFind} object
#' @name loadObject
#'
#' @useDynLib scASfind
load.serialized.object <- function(filename) {
  object <- readRDS(filename)
  # Deserialize object
  object@index <- new(EliasFanoDB)
  success <- object@index$loadByteStream(object@serialized)
  object@serialized <- raw()
  gc()
  ## Dirty hack so we do not have to rebuild again every scfind index
  if (is.null(object@metadata)) {
    object@metadata <- list()
  }
  return(object)
}

#' @rdname loadObject
#' @aliases loadObject
setMethod("loadObject", definition = load.serialized.object)



#' Merges an external index into the existing object
#'
#' This function is useful to merge \code{SCFind} indices.
#' After this operation object that was merged can be discarded.
#'
#' The only semantic limitation for merging two databases is to
#' have different dataset names in the two different indices.
#' If that is not case user may run into problems masking datasets
#' from the different datasets while there is a possibility of having
#' different cell types under the same name. This will most likely cause
#' undefined behavior during queries.
#'
#' @param object the root scfind object
#' @param new.object external scfind object to be merged
#'
#' @name mergeDataset
#' @return the new extended object
#'
merge.dataset.from.object <- function(object, new.object) {
  common.datasets <- intersect(new.object@datasets, object@datasets)

  message(paste("Merging", new.object@datasets))
  if (length(common.datasets) != 0) {
    warning("Common dataset names exist, undefined merging behavior, please fix this...")
  }

  object@index$mergeDB(new.object@index)
  object@datasets <- c(object@datasets, new.object@datasets)
  return(object)
}

#' Used to merge multiple eliasfanoDB
#'
#'
#' @rdname mergeDataset
#' @aliases mergeDataset mergeObjects
setMethod(
  "mergeDataset",
  signature(
    object = "SCFind",
    new.object = "SCFind"
  ),
  merge.dataset.from.object
)



#' Query Optimization Function for SCFind objects.
#'
#' This function can be used with quite long gene lists
#' that otherwise would have no cell hits in the database
#'
#' @param object SCFind object
#' @param gene.list A list of nGenes existing in the database
#' @param datasets the datasets of the objects to be considered
#' @param log.message whether to print a verbose message
#'
#' @name markerNodes
#' @return hierarchical list of queries and their respective scores
find.marker.nodes <- function(object, gene.list, datasets, log.message = 0) {
  datasets <- select.datasets(object, datasets)
  results <- object@index$findMarkerGenes(as.character(caseCorrect(object, gene.list)), as.character(datasets), 5, log.message)
  return(results)
}


#' @rdname markerNodes
#' @aliases markerNodes
setMethod(
  "markerNodes",
  signature(
    object = "SCFind",
    gene.list = "character"
  ),
  find.marker.nodes
)

#' Find marker genes for a specific cell type
#'
#' @name cellTypeMarkers
#'
#' @param object SCFind object
#' @param cell.types the cell types that we want to extract the marker genes
#' @param background.cell.types the universe of cell.types to consider
#' @param top.k how many genes to retrieve
#' @param sort.field the dataframe will be sorted according to this field
#'
#' @return a data.frame that each row represent a gene score for a specific cell type
cell.type.marker <- function(object, cell.types, background.cell.types, top.k, sort.field) {
  if (missing(background.cell.types)) {
    background.cell.types <- cellTypeNames(object)
  }

  all.cell.types <- object@index$cellTypeMarkers(cell.types, background.cell.types)

  if (!(sort.field %in% colnames(all.cell.types))) {
    message(paste("Column", sort.field, "not found"))
    sort.field <- "f1"
  }
  all.cell.types <- all.cell.types[order(all.cell.types[[sort.field]], decreasing = T)[1:top.k], ]

  if (!is.null(object@metadata$node_list)) {
    all.cell.types <- merge(all.cell.types, object@metadata$node_list, by.x = "nodes", by.y = "Node_id", all.x = TRUE, all.y. = FALSE)
  }

  return(all.cell.types)
}


#' @rdname cellTypeMarkers
#' @aliases cellTypeMarkers
setMethod(
  "cellTypeMarkers",
  signature(
    object = "SCFind",
    cell.types = "character"
  ),
  cell.type.marker
)


#' Return a vector with all existing cell type names in the database
#'
#' @name cellTypeNames
#' @param object SCFind object
#' @param datasets individual datasets to consider
#'
#' @return a character list
get.cell.types.names <- function(object, datasets) {
  if (missing(datasets)) {
    return(object@index$getCellTypes())
  } else {
    return(object@index$getCellTypes()[lapply(strsplit(object@index$getCellTypes(), "\\."), `[[`, 1) %in% datasets])
  }
}
#' @rdname cellTypeNames
#' @aliases cellTypeNames
setMethod(
  "cellTypeNames",
  signature(
    object = "SCFind"
  ),
  get.cell.types.names
)


#' Evaluate a user specific query by calculating the precision recall metrics
#'
#' @name evaluateMarkers
#' @param object the \code{SCFind} object
#' @param gene.list the list of genes to be evaluated
#' @param cell.types a list of cell types for the list to evaluated
#' @param background.cell.types the universe of cell.types to consider
#' @param sort.field the dataframe will be sorted according to this field
#'
#' @return a DataFrame that each row represent a gene score for a specific cell type
#'
evaluate.cell.type.markers <- function(object, gene.list, cell.types, background.cell.types, sort.field) {
  if (missing(background.cell.types)) {
    message("Considering the whole DB..")
    background.cell.types <- cellTypeNames(object)
  }
  all.cell.types <- object@index$evaluateCellTypeMarkers(cell.types, caseCorrect(object, gene.list), background.cell.types)

  if (!(sort.field %in% colnames(all.cell.types))) {
    message(paste("Column", sort.field, "not found"))
    sort.field <- "f1"
  }
  all.cell.types <- all.cell.types[order(all.cell.types[[sort.field]]), ]
  return(all.cell.types)
}

#' @rdname evaluateMarkers
#' @aliases evaluateMarkers
setMethod(
  "evaluateMarkers",
  signature(
    object = "SCFind",
    gene.list = "character"
  ),
  evaluate.cell.type.markers
)





#' Find cell types associated with a given gene list. All cells
#' returned express all of the genes in the given gene list
#'
#' @param object the \code{SCFind} object
#' @param gene.list genes to be searched in the gene.index
#' (Operators: "-gene" to exclude a gene | "*gene" either gene is expressed
#' "*-gene" either gene is expressed to be excluded)
#' @param datasets the datasets that will be considered
#'
#'
#' @importFrom utils setTxtProgressBar stack unstack tail
#'
#' @name findCellTypes
#' @return a named numeric vector containing p-values
findCellTypes.geneList <- function(object, gene.list, datasets) {
  datasets <- if (missing(datasets)) object@datasets else select.datasets(object, datasets)

  if (length(grep("^-|^\\*", gene.list)) == 0) {
    return(object@index$findCellTypes(caseCorrect(object, gene.list), datasets))
  } else {
    pos <- caseCorrect(object, grep("^[^-\\*]", gene.list, value = T))
    excl.or <- grep("^-\\*|^\\*-", gene.list, value = T)
    or <- caseCorrect(object, sub("\\*", "", setdiff(grep("^\\*", gene.list, value = T), excl.or)))
    excl <- caseCorrect(object, sub("-", "", setdiff(grep("^-", gene.list, value = T), excl.or)))
    excl.or <- caseCorrect(object, sub("\\*-||-\\*", "", grep("^-\\*|^\\*-", gene.list, value = T)))

    if (length(c(intersect(pos, or), intersect(pos, excl), intersect(pos, excl.or), intersect(or, excl), intersect(or, excl.or), intersect(excl, excl.or))) != 0) {
      message("Warning: Same gene labeled with different operators!")
      message("There is a priority to handle operators:")
      message(paste(
        "Cells with", paste(pos, collapse = " ^ "), "expression will be included.",
        if (length(or) != 0) "Then cells with", paste(or, collapse = " v "), "expression will be included."
      ))
      message(paste(
        "The result will be excluded by", paste(excl, collapse = " ^ "),
        if (length(excl.or != 0)) paste("and further be excluded by", paste(excl.or, collapse = " v "))
      ))
      cat("\n")
    }


    cell.to.id <- NULL

    # Using pair.id to create unique variable for each cell by pairing cell types to cell ID

    if (length(pos) == 0 && length(or) == 0 && (length(excl) != 0 || length(excl.or) != 0)) {
      # When no positive selection, include all cells
      cell.to.id <- lapply(as.list(object@index$getCellTypeSupport(cellTypeNames(object, datasets))), seq)
      names(cell.to.id) <- cellTypeNames(object, datasets)
      cell.to.id <- pair.id(cell.to.id)
    }


    if (length(or) != 0) {
      # Include any cell expresses gene in OR condition
      gene.or <- c()
      for (i in 1:length(or))
      {
        tmp.id <- pair.id(object@index$findCellTypes(c(pos, or[i]), datasets))
        if (length(pos) != 0 && !is.null(tmp.id)) message(paste("Found", length(tmp.id), if (length(tmp.id) > 1) "cells" else "cell", "co-expressing", paste(c(pos, or[i]), collapse = " ^ ")))
        if (!is.null(tmp.id)) {
          cell.to.id <- unique(c(cell.to.id, tmp.id))
          # Store used query
          gene.or <- c(gene.or, or[i])
        } else {
          cell.to.id <- cell.to.id
        }
      }
      if (length(pos) == 0 && length(gene.or) != 0) message(paste("Found", length(cell.to.id), if (length(cell.to.id) > 1) "cells" else "cell", "expressing", paste(gene.or, collapse = " v ")))
    } else {
      cell.to.id <- if (length(pos) != 0) pair.id(object@index$findCellTypes(pos, datasets)) else cell.to.id
      if (length(pos) != 0) message(paste("Found", length(cell.to.id), if (length(pos) > 1) "cells co-expressing" else "cell expressing", paste(pos, collapse = " ^ ")))
    }

    count.cell <- length(cell.to.id)
    gene.excl <- NULL

    if (length(excl.or) != 0) {
      # Negative select cell in OR condition
      for (i in 1:length(excl.or))
      {
        ex.tmp.id <- pair.id(object@index$findCellTypes(c(excl, excl.or[i]), datasets))

        message(paste(
          "Excluded", sum(cell.to.id %in% ex.tmp.id),
          if (sum(cell.to.id %in% ex.tmp.id) > 1) "cells" else "cell",
          if (length(excl) != 0) paste("co-expressing", paste(c(excl, excl.or[i]), collapse = " ^ ")) else paste("expressing", excl.or[i])
        ))

        if (!is.null(ex.tmp.id)) {
          cell.to.id <- setdiff(cell.to.id, ex.tmp.id)
          gene.excl <- c(gene.excl, excl.or[i])
        } else {
          cell.to.id <- cell.to.id
        }
      }
      count.cell <- count.cell - length(cell.to.id)
      if (count.cell > 0 && length(gene.excl) == 0) message("Excluded", count.cell, if (count.cell > 1) "cells" else "cell", "expressing", paste(excl, collapse = " ^ "))
    } else {
      if (length(excl) != 0) {
        # Negative selection
        cell.to.id <- setdiff(cell.to.id, pair.id(object@index$findCellTypes(excl, datasets)))
        count.cell <- count.cell - length(cell.to.id)
        if (count.cell > 0) message(paste("Excluded", count.cell, if (count.cell > 1) "cells" else "cell", if (length(excl) > 1) "co-expressing" else "expressing", paste(excl, collapse = " ^ "))) else message("No Cell Is Excluded!")
      }
    }

    # Generate a new list
    df <- do.call(rbind, strsplit(as.character(cell.to.id), "#"))
    if (!is.null(df)) {
      result <- as.list(setNames(as.numeric(split(df[, 2], seq(nrow(df)))), df[, 1]))
      if (length(unique(df[, 1])) == nrow(df)) {
        return(result)
      } else {
        if (length(unique(names(result))) == 1) {
          tmp <- list(stack(result)$values)
          names(tmp) <- unique(names(result))
          return(tmp)
        } else {
          return(unstack(stack(result)))
        }
      }
    } else {
      message("No Cell Is Found!")
      return(list())
    }
  }
}


#' @rdname findCellTypes
#' @aliases findCellTypes
setMethod(
  "findCellTypes",
  signature(
    object = "SCFind",
    gene.list = "character"
  ),
  findCellTypes.geneList
)

#' Get all nodes in the database
#'
#' @name scfindNodes
#'
#' @param object the \code{scfind} object
#'
#' @return the list of genes present in the database
scfind.get.nodes.in.db <- function(object) {
  return(object@index$genes())
}


#' @rdname scfindNodes
#' @aliases scfindNodes
setMethod("scfindNodes", signature(object = "SCFind"), scfind.get.nodes.in.db)


#' Find out how many cell-types each gene is found
#'
#' @param object the \code{SCFind} object
#' @param gene.list genes to be searched in the gene.index
#' @param datasets the datasets that will be considered
#' @param min.cells threshold of cell hit of a cell type
#' @param min.fraction portion of total cell as threshold
#'
#' @name findCellTypeSpecificities
#' @return the list of number of cell type for each gene
cell.type.specificity <- function(object, gene.list, datasets, min.cells = 10, min.fraction = .25) {
  if (min.fraction >= 1 || min.fraction <= 0) stop("min.fraction reached limit, please use values > 0 and < 1.0.") else message("Calculating cell-types for each gene...")
  datasets <- if (missing(datasets)) object@datasets else select.datasets(object, datasets)
  if (missing(gene.list)) {
    res <- object@index$geneSupportInCellTypes(object@index$genes(), datasets)
  } else {
    gene.list <- caseCorrect(object, gene.list)
    res <- object@index$geneSupportInCellTypes(gene.list, datasets)
  }

  res.tissue <- res
  names(res.tissue) <- gsub("\\.", "#", names(res.tissue))
  df <- cbind(stack(res), stack(unlist(res.tissue)))
  # df[,4] <- sub("^[^.]+\\.", "", df[,4])
  df[, 1] <- object@index$getCellTypeSupport(sub("^[^.]+\\.", "", df[, 4])) * min.fraction
  if (length(which(df[, 1] < min.cells)) != 0) df[which(df[, 1] < min.cells), 1] <- min.cells
  if (nrow(df) != 0) df <- df[which(df[, 3] > df[, 1]), ] else return(split(rep(0, length(gene.list)), gene.list))
  if (nrow(df) != 0) {
    return(as.list(summary(df[, 2], maxsum = nrow(df))))
  } else {
    return(split(rep(0, length(gene.list)), gene.list))
  }
}

#' @rdname findCellTypeSpecificities
#' @aliases findCellTypeSpecificities
setMethod(
  "findCellTypeSpecificities",
  signature(object = "SCFind"),
  cell.type.specificity
)


#' Find out how many tissues each gene is found
#'
#' @param object the \code{SCFind} object
#' @param gene.list genes to be searched in the gene.index
#' @param min.cells threshold of cell hit of a tissue
#'
#' @name findTissueSpecificities
#' @return the list of number of tissue for each gene
tissue.specificity <- function(object, gene.list, min.cells = 10) {
  if (length(object@datasets) <= 1) stop("Index contains 1 dataset only.") else message("Calculating tissues for each gene...")
  if (missing(gene.list)) {
    res <- object@index$geneSupportInCellTypes(object@index$genes(), object@datasets)
  } else {
    gene.list <- caseCorrect(object, gene.list)
    res <- object@index$geneSupportInCellTypes(gene.list, object@datasets)
  }

  if (length(res) > 0) res.tissue <- res else return(split(rep(0, length(gene.list)), gene.list))
  names(res.tissue) <- gsub("\\.", "#", names(res.tissue))
  df <- cbind(stack(res), stack(unlist(res.tissue)))
  df[, 5] <- gsub("^[^.]*\\.([^.]*)\\..*$", "\\1", df[, 4])
  df <- aggregate(df[, 1], by = list(df[, 5], df[, 2]), FUN = sum)
  df <- df[which(df[, 3] > min.cells), ]

  if (nrow(df) != 0) {
    return(as.list(summary(df[, 2], maxsum = nrow(df))))
  } else {
    return(split(rep(0, length(gene.list)), gene.list))
  }
}

#' @rdname findTissueSpecificities
#' @aliases findTissueSpecificities
setMethod(
  "findTissueSpecificities",
  signature(object = "SCFind"),
  tissue.specificity
)

#' Find the set of genes that are ubiquitously expressed in a query of cell types
#'
#' @param object the \code{SCFind} object
#' @param cell.types a list of cell types for the list to evaluated
#' @param min.recall threshold of minimun recall value
#' @param max.genes threshold of number of genes to be considered for each cell type
#'
#' @importFrom utils txtProgressBar
#' @name findHouseKeepingNodes
#' @return the list of gene that ubiquitously expressed in a query of cell types
#'
house.keeping.nodes <- function(object, cell.types, min.recall = .5, max.genes = 1000) {
  if (min.recall >= 1 || min.recall <= 0) stop("min.recall reached limit, please use values > 0 and < 1.0.")
  if (max.genes > length(object@index$genes())) stop(paste("max.genes exceeded limit, please use values > 0 and < ", length(object@index$genes()))) else message("Searching for house keeping node...")
  df <- cellTypeMarkers(object, cell.types[1], top.k = max.genes, sort.field = "recall")
  house.keeping.genes <- df$nodes[which(df$recall > min.recall)]

  for (i in 2:length(cell.types)) {
    setTxtProgressBar(txtProgressBar(1, length(cell.types), style = 3), i)
    df <- cellTypeMarkers(object, cell.types[i], top.k = max.genes, sort.field = "recall")
    house.keeping.genes <- intersect(house.keeping.genes, df$nodes[which(df$recall > min.recall)])
    if (length(house.keeping.genes) == 0) {
      stop("No house keeping node is found.")
    }
  }
  cat("\n")
  return(house.keeping.genes)
}


#' @rdname findHouseKeepingNodes
#' @aliases findHouseKeepinNodes
setMethod(
  "findHouseKeepingNodes",
  signature(
    object = "SCFind",
    cell.types = "character"
  ),
  house.keeping.nodes
)

#'  Find the signature genes for a cell-type
#'
#' @param object the \code{SCFind} object
#' @param cell.types a list of cell types for the list to evaluated
#' @param max.genes threshold of number of genes to be considered for each cell type
#' @param min.cells threshold of cell hit of a tissue
#' @param max.pval threshold of p-value
#'
#' @importFrom utils setTxtProgressBar
#'
#' @name findNodeSignatures
#' @return the list of gene signatures in a query of cell types
#'
node.signatures <- function(object, cell.types, max.genes = 1000, min.cells = 10, max.pval = 0) {
  message("Searching for node signatures...")
  cell.types.all <- if (missing(cell.types)) object@index$getCellTypes() else cellTypeNames(object)[tolower(cellTypeNames(object)) %in% tolower(cell.types)]
  signatures <- list()
  if (length(cell.types.all) != 0) {
    for (i in 1:length(cell.types.all)) {
      if (i > 1) setTxtProgressBar(txtProgressBar(1, length(cell.types.all), style = 3), i)
      signatures[[cell.types.all[i]]] <- find.signature(object, cell.types.all[i], max.genes = max.genes, min.cells = min.cells, max.pval = max.pval)
    }
    cat("\n")
    return(signatures)
  } else {
    return(message(paste0("Ignored ", toString(cell.types), ". Cell type not found in index.")))
  }
}

#' @rdname findNodeSignatures
#' @aliases findNodeSignatures
setMethod(
  "findNodeSignatures",
  signature(object = "SCFind"),
  node.signatures
)

#'  Look at all other genes and rank them based on the similarity of their expression pattern to the pattern defined by the gene query
#'
#' @param object the \code{SCFind} object
#' @param node.list genes to be searched in the gene.index
#' @param datasets the datasets that will be considered
#' @param top.k how many genes to retrieve
#'
#' @importFrom utils setTxtProgressBar
#' @name findSimilarNodes
#' @return the list of genes and their similarities presented in Jaccard indices
#'
similar.nodes <- function(object, node.list, datasets, top.k = 5) {
  message("Searching for genes with similar pattern...")
  datasets <- if (missing(datasets)) object@datasets else select.datasets(object, datasets)
  node.list <- caseCorrect(object, node.list)
  e <- object@index$findCellTypes(node.list, datasets) # the cells expressing the genes in node.list
  n.e <- length(unlist(e))
  if (n.e > 0) {
    gene.names <- setdiff(object@index$genes(), node.list)
    similarities <- rep(0, length(gene.names))
    ns <- rep(0, length(gene.names))
    ms <- rep(0, length(gene.names))
    for (i in 1:length(gene.names)) {
      setTxtProgressBar(txtProgressBar(1, length(gene.names), style = 3), i)
      f <- object@index$findCellTypes(gene.names[i], datasets) # find expression pattern of other gene
      if (length(f) > 0) {
        m <- rep(0, length(e))
        for (j in 1:length(names(e))) {
          m[j] <- length(intersect(e[[j]], f[[names(e)[j]]]))
        }
        # calculate the Jaccard index for the similarity of the cells expressing the gene
        similarities[i] <- sum(m) / (n.e + length(unlist(f)) - sum(m))
        ns[i] <- length(unlist(f))
        ms[i] <- sum(m)
      }
    }
    cat("\n")
    r <- sort(similarities, index.return = T)
    inds <- tail(r$ix, top.k)
    res <- data.frame("gene" = gene.names[inds], "Jaccard" = similarities[inds], "overlap" = ms[inds], "n" = ns[inds])
    return(res)
  } else {
    message(paste("Cannot find cell expressing", toString(node.list), "in the index."))
    return(c())
  }
}


#' @rdname findSimilarNodes
#' @aliases findSimilarNodes
setMethod(
  "findSimilarNodes",
  signature(
    object = "SCFind",
    node.list = "character"
  ),
  similar.nodes
)
