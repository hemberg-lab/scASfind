
#' merge two elias fano indices
#'
#' @param efdb.root the root index
#' @param efdb the index to be merged
#'
#' @importFrom hash values keys
#'
mergeIndices <- function(efdb.root, efdb) {
  genes.to.merge <- intersect(keys(efdb.root), keys(efdb))
  new.genes <- setdiff(keys(efdb), keys(efdb.root))

  ## print(paste(length(new.genes),length(genes.to.merge)))
  for (gene in new.genes)
  {
    efdb.root[[gene]] <- efdb[[gene]]
  }
  for (gene in genes.to.merge)
  {
    ## maybe we want to merge??
    efdb.root[[gene]][keys(efdb[[gene]])] <- values(efdb[[gene]], simplify = F)
  }
  return(efdb.root)
}


contigency.table <- function(query.results) {
  data <- as.data.frame(lapply(values(query.results), length))
  names(data) <- keys(query.results)
  return(data)
}

caseCorrect <- function(object, gene.list) {
  gene.list <- gene.list[gene.list != ""]

  if (length(gene.list) != 0) {
    gene.corr <- object@index$genes()[match(tolower(gene.list), tolower(object@index$genes()), nomatch = 0)]

    if (length(setdiff(tolower(gene.list), tolower(gene.corr))) != 0) message(paste0("Ignored ", toString(setdiff(gene.list, gene.list[match(tolower(gene.corr), tolower(gene.list), nomatch = 0)])), ". Not found in the index"))

    return(unique(gene.corr))
  } else {
    return(c())
  }
}



#' @importFrom stats aggregate p.adjust phyper setNames
phyper.test <- function(object, result, datasets) {
  df <- query.result.as.dataframe(result)
  datasets <- select.datasets(object, datasets)

  cell.types.df <- aggregate(cell_id ~ cell_type, df, FUN = length)
  colnames(cell.types.df)[colnames(cell.types.df) == "cell_id"] <- "cell_hits"
  cell.types.df$total_cells <- object@index$getCellTypeSupport(cell.types.df$cell_type)
  query.hits <- nrow(df)

  cell.types.df$pval <- p.adjust(1 - phyper(
    cell.types.df$cell_hits, # total observed successes ( query.hits for cell type)
    cell.types.df$total_cells, # total successes ( cell type size )
    sum(cell.types.df$total_cells) - cell.types.df$total_cells, # total failures( total cells excluding cell type)
    query.hits # sample size
  ), n = length(cellTypeNames(object, datasets)))


  return(cell.types.df)
}

query.result.as.dataframe <- function(query.result) {
  if (is.data.frame(query.result)) {
    return(query.result)
  }
  if (length(query.result) == 0) {
    return(data.frame(cell_type = c(), cell_id = c()))
  } else {
    result <- setNames(unlist(query.result, use.names = F), rep(names(query.result), lengths(query.result)))
    return(data.frame(cell_type = names(result), cell_id = result))
  }
}

select.datasets <- function(object, datasets) {
  if (missing(datasets)) {
    ## Select all available datasets
    datasets <- object@datasets
  } else {
    ## datasets should not be a superset of the data
    if (length(setdiff(datasets, object@datasets)) != 0) {
      stop(paste("Dataset", setdiff(datasets, object@datasets), "does not exist in the database"))
    }
  }
  return(datasets)
}



scfind.get.genes.in.db <- function(object) {
  return(object@index$genes())
}

pair.id <- function(cell.list = list()) {
  if (length(cell.list) == 0) {
    return(c())
  } else {
    pair.vec <- stack(cell.list)
    return(paste0(pair.vec$ind, "#", pair.vec$values))
  }
}

find.signature <- function(object, cell.type, max.genes = 1000, min.cells = 10, max.pval = 0) {
  # Use this method to find a gene signature for a cell-type.
  # We do this by ranking genes by recall and then adding genes to the query until we exceed a target p-value threshold or until a minimum number of cells is returned from the query
  df <- cellTypeMarkers(object, cell.type, top.k = max.genes, sort.field = "recall")
  genes <- as.character(df$nodes)
  genes.list <- c()
  thres <- max(c(min.cells, object@index$getCellTypeMeta(cell.type)$total_cells))
  for (j in 1:dim(df)[1]) {
    res <- hyperQueryCellTypes(object, c(genes.list, genes[j]))
    if (dim(res)[1] > 0) {
      ind <- which(res[, 1] == cell.type)
      if (length(ind) == 0) {
        break
      } else {
        if (res[ind, 4] > max.pval | res[ind, 2] < thres) {
          break
        }
      }
    }
    genes.list <- c(genes.list, genes[j])
  }
  return(genes.list)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}



#' reorder correlation matrix
#'
#' @param cormat correlation matrix
#'
#' @import reshape2
#'
reorder_cormat <- function(cormat) {
  # Use correlation between variables as distance
  dd <- as.dist((1 - cormat) / 2)
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
}

#' Get lower triangle of the correlation matrix
#'
#' @param cormat correlation matrix
#'
#' @import reshape2
#'
get_lower_tri <- function(cormat) {
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

#' # Get upper triangle of the correlation matrix
#'
#' @param cormat correlation matrix
#'
#' @import reshape2
#'
get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}

#' This function gets the raw PSI of a node in a cell type
#'
#' @param object the \code{SCFind} object
#' @param gene.list several nodes that we wish to get the raw PSI
#' @param cell.type cell type tp query, can only query one cell type at once
#' @return a dataframe that contains raw psi value in the queried cell type of the gene.list
#'
#'

get.cell.type.raw.psi <- function(object, gene.list, cell.type) {
  if (any(grepl("above|below", object@datasets)) == FALSE) {
    warning("index must contain two datasets with above and below events")
  }
  if (length(cell.type) != 1) {
    warning("can only query one cell type at once")
  }

  ## the reverse process of building index-ready PSI matrix

  ## prepare above and below scaled matrix by querying index

  cell.type.above <- paste("above", cell.type, sep = ".")

  cell.type.below <- paste("below", cell.type, sep = ".")

  cells_psi_above <- as.data.frame(object@index$getCellTypeExpression(cell.type.above))

  scaled_psi_above <- cells_psi_above[which(rownames(cells_psi_above) %in% gene.list), ]

  cells_psi_below <- as.data.frame(object@index$getCellTypeExpression(cell.type.below))

  scaled_psi_below <- cells_psi_below[which(rownames(cells_psi_below) %in% gene.list), ]

  # prepare tissue-mean and gene list for the raw psi matrix

  mean <- object@metadata$stats[gene.list, ]

  basis <- data.frame("gene_id" = gene.list, row.names = gene.list)

  # process above and below to the raw matrix

  dd <- data.frame(matrix(ncol = ncol(scaled_psi_above), nrow = length(gene.list[!(gene.list %in% rownames(scaled_psi_above))])),
    row.names = gene.list[!(gene.list %in% rownames(scaled_psi_above))]
  )
  colnames(dd) <- colnames(scaled_psi_above)

  # scale back

  above <- rbind(dd, scaled_psi_above)[gene.list, ] * 0.01


  dd2 <- data.frame(matrix(ncol = ncol(scaled_psi_below), nrow = length(gene.list[!(gene.list %in% rownames(scaled_psi_below))])),
    row.names = gene.list[!(gene.list %in% rownames(scaled_psi_below))]
  )
  colnames(dd2) <- colnames(scaled_psi_below)

  below <- rbind(dd2, scaled_psi_below)[gene.list, ] * 0.01

  # merge above and below matrix element-wise

  for (i in seq(1, nrow(below))) {
    for (j in seq(1, ncol(below))) {
      if (!is.na(below[i, j]) & (above[i, j] == 0)) {
        above[i, j] <- -1 * below[i, j]
      }
    }
  }
  # identify nodes that are not sufficiently different from mean hence got value 0


  diff_cut <- object@metadata$diff_cut

  diff_check <- diff_cut[gene.list, grepl(cell.type, colnames(diff_cut))]

  for (i in seq(1, nrow(above))) {
    for (j in seq(1, ncol(above))) {
      if ((diff_check[i, j] == 1)) { # as long as diff_cut is 1, mean that the cell have PSI within cut-off range of mean value and was removed

        above[i, j] <- mean[i, "mean"] # as an approximation, set PSI to tissue mean
      } else if ((!is.na(above[i, j])) && (above[i, j] == 0) && (diff_check[i, j] == 0)) {
        above[i, j] <- NA # these cells are truely NA
      } else if ((!is.na(above[i, j])) && (above[i, j] != 0)) { # cells with a sufficiently different PSI, revert the centralization by adding mean

        above[i, j] <- above[i, j] + mean[i, "mean"]
      }
    }
  }

  # due to imprecise retrieval of PSI there could be values slightly beyond limit
  above[above > 1] <- 1
  above[above < 0] <- 0

  return(above)
}
