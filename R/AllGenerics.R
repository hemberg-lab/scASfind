#' The scfind main class object
#'
#' @export
setClass(
  "SCFind",
  representation(
    index = "ANY", # Index is type of Rcpp_EliasFanoDB but this removes the warning
    datasets = "character",
    serialized = "raw",
    metadata = "list"
  )
)

#'
#' @export
setGeneric(
  name = "buildAltSpliceIndex",
  def = function(psival,
                 metadata,
                 dataset.name,
                 column.label = "primary_type",
                 qb = 2) {
    standardGeneric("buildAltSpliceIndex")
  }
)

#'
#' @export
setGeneric(
  name = "addIndexMeta",
  def = function(object, stats, node_list, diff_cut) {
    standardGeneric("addIndexMeta")
  }
)

#'
#' @export
setGeneric(
  name = "hyperQueryCellTypes",
  def = function(object,
                 node.list,
                 datasets) {
    standardGeneric("hyperQueryCellTypes")
  }
)

#'
#' @export
setGeneric(
  name = "nodeDetails",
  def = function(object,
                 node.list) {
    standardGeneric("nodeDetails")
  }
)

#'
#' @export
setGeneric(
  name = "geneNodes",
  def = function(object,
                 gene.list,
                 query.type) {
    standardGeneric("geneNodes")
  }
)

#'
#' @export
setGeneric(
  name = "findNodeSets",
  def = function(object,
                 gene.list,
                 query.type,
                 node.types) {
    standardGeneric("findNodeSets")
  }
)


#'
#' @export
setGeneric(
  name = "getCoordinatedNodes",
  def = function(object,
                 gene.name) {
    standardGeneric("getCoordinatedNodes")
  }
)



#'
#' @export
setGeneric(
  name = "getRawPsi",
  def = function(object,
                 node.list,
                 cell.types,
                 verbose) {
    standardGeneric("getRawPsi")
  }
)


#'
#' @export
setGeneric(
  name = "plotRawPsiCorr",
  def = function(raw_psi,
                 node.list,
                 cell.types) {
    standardGeneric("plotRawPsiCorr")
  }
)


#'
#' @export
setGeneric(
  name = "plotRawPsiHeatmap",
  def = function(raw_psi,
                 node.list,
                 cell.types) {
    standardGeneric("plotRawPsiHeatmap")
  }
)





#' @export
setGeneric(name = "mergeDataset", def = function(object, new.object) {
  standardGeneric("mergeDataset")
})


#' queries cells that contain all the genes from the list
#' @export
setGeneric(name = "findCellTypes", function(object, gene.list, datasets) {
  standardGeneric("findCellTypes")
})


#' return all the gene markers for a specified cell.type
#'
#' @export
#'
setGeneric(name = "cellTypeMarkers", function(object,
                                              cell.types,
                                              background.cell.types,
                                              top.k = 5,
                                              sort.field = "f1") {
  standardGeneric("cellTypeMarkers")
})


#' @export
#'
setGeneric(name = "cellTypeNames", function(object, datasets) {
  standardGeneric("cellTypeNames")
})

#' @export
#'
setGeneric(name = "evaluateMarkers", function(object,
                                              gene.list,
                                              cell.types,
                                              background.cell.types,
                                              sort.field = "f1") {
  standardGeneric("evaluateMarkers")
})



#' Generic to be used instead of readRDS
#' @export
setGeneric(name = "loadObject", function(filename) {
  standardGeneric("loadObject")
})

#' Generic to be used instead of saveRDS
#'
#' @export
setGeneric(name = "saveObject", function(object, file) {
  standardGeneric("saveObject")
})


#' Performs query optimization and return the best candidate gene sets
#'
#' @export

setGeneric(name = "markerNodes", function(object, node.list, datasets, log.message = 0) {
  standardGeneric("markerNodes")
})


#' @export
setGeneric(name = "scfindNodes", function(object) {
  standardGeneric("scfindNodes")
})


#' Find out how many cell-types each gene is found
#'
#' @export
setGeneric(name = "findCellTypeSpecificities", function(object,
                                                        gene.list,
                                                        datasets,
                                                        min.cells = 10,
                                                        min.fraction = .25) {
  standardGeneric("findCellTypeSpecificities")
})

#' Find out how many tissues each gene is found
#'
#' @export
setGeneric(name = "findTissueSpecificities", function(object,
                                                      gene.list,
                                                      min.cells = 10) {
  standardGeneric("findTissueSpecificities")
})

#' Find the set of genes that are ubiquitously expressed in a query of cell types
#'
#' @export
setGeneric(name = "findHouseKeepingNodes", function(object,
                                                    cell.types,
                                                    min.recall = .5,
                                                    max.genes = 1000) {
  standardGeneric("findHouseKeepingNodes")
})

#' Find the signature genes for a cell-type
#'
#' @export
setGeneric(name = "findNodeSignatures", function(object,
                                                 cell.types,
                                                 max.genes = 1000,
                                                 min.cells = 10,
                                                 max.pval = 0) {
  standardGeneric("findNodeSignatures")
})

#' Look at all other genes and rank them based on the similarity of their expression pattern to the pattern defined by the gene query
#'
#' @export
setGeneric(name = "findSimilarNodes", function(object, node.list, datasets, top.k = 5) {
  standardGeneric("findSimilarNodes")
})
