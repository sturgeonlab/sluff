#' Generic to extract pseudotime from CDS object
#'
#' @param x A cell_data_set object.
#' @param reduction_method Reduced dimension to extract pseudotime for.
#'
#' @examples
#'  \donttest{
#'    cell_metadata <- readRDS(system.file('extdata',
#'                             'worm_embryo/worm_embryo_coldata.rds',
#'                             package='monocle3'))
#'    gene_metadata <- readRDS(system.file('extdata',
#'                                         'worm_embryo/worm_embryo_rowdata.rds',
#'                                         package='monocle3'))
#'    expression_matrix <- readRDS(system.file('extdata',
#'                                             'worm_embryo/worm_embryo_expression_matrix.rds',
#'                                              package='monocle3'))
#'
#'    cds <- new_cell_data_set(expression_data=expression_matrix,
#'                             cell_metadata=cell_metadata,
#'                             gene_metadata=gene_metadata)
#'
#'    cds <- preprocess_cds(cds, num_dim=50)
#'    cds <- align_cds(cds, alignment_group = "batch",
#'                     residual_model_formula_str = "~ bg.300.loading + bg.400.loading +
#'                     bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading +
#'                     bg.b02.loading")
#'    cds <- reduce_dimension(cds) 
#'    ciliated_genes <- c("che-1", "hlh-17", "nhr-6", "dmd-6", "ceh-36", "ham-1")
#'    cds <- cluster_cells(cds)
#'    cds <- learn_graph(cds)
#'    cds <- order_cells(cds,root_pr_nodes='Y_27')
#'    ps_tim <- pseudotime(cds)
#'  }
#' 
#' @return Pseudotime values.
#'
#' @export
setGeneric("pseudotime", function(x, reduction_method = "UMAP")
  standardGeneric("pseudotime"))


#' Method to extract pseudotime from CDS object
#' @param x A cell_data_set object.
#' @param reduction_method Reduced dimension to extract clusters for.
#'
#' @return Pseudotime values.
#'
#' @export
setMethod("pseudotime", "cell_data_set",
          function(x, reduction_method = "UMAP") {
            value <- x@principal_graph_aux[[
              reduction_method]]$pseudotime[colnames(exprs(x))]
            if (is.null(value)) {
              stop("No pseudotime calculated for reduction_method = ",
                          reduction_method, ". Please first run ",
                          "order_cells with reduction_method = ",
                          reduction_method, ".")
            }
            return(value)
          })

#' Generic to extract clusters from CDS object
#'
#' @param x A cell_data_set object.
#' @param reduction_method Reduced dimension to extract clusters for.
#'
#' @examples
#'   \donttest{
#'     cell_metadata <- readRDS(system.file('extdata',
#'                                          'worm_embryo/worm_embryo_coldata.rds',
#'                                          package='monocle3'))
#'     gene_metadata <- readRDS(system.file('extdata',
#'                                          'worm_embryo/worm_embryo_rowdata.rds',
#'                                          package='monocle3'))
#'     expression_matrix <- readRDS(system.file('extdata',
#'                                              'worm_embryo/worm_embryo_expression_matrix.rds',
#'                                              package='monocle3'))
#'     cds <- new_cell_data_set(expression_data=expression_matrix,
#'                              cell_metadata=cell_metadata,
#'                              gene_metadata=gene_metadata)
#'
#'     cds <- preprocess_cds(cds)
#'     cds <- reduce_dimension(cds)
#'     cds <- cluster_cells(cds)
#'     clusters_factors <- clusters(cds, "UMAP")
#'   }
#'
#' @return Clusters.
#'
#' @export
setGeneric("clusters", function(x, reduction_method = "UMAP")
  standardGeneric("clusters"))

#' Method to extract clusters from CDS object
#' @param x A cell_data_set object.
#' @param reduction_method Reduced dimension to extract clusters for.
#'
#' @return Clusters.
#'
#' @export
setMethod("clusters", "cell_data_set",
          function(x, reduction_method = "UMAP") {
            value <- x@clusters[[
              reduction_method]]$clusters[colnames(exprs(x))]
            if (is.null(value)) {
              stop("No clusters calculated for reduction_method = ",
                          reduction_method, ". Please first run ",
                          "cluster_cells with reduction_method = ",
                          reduction_method, ".")
            }
            return(value)
          })

#' Generic to extract partitions from CDS object
#'
#' @param x A cell_data_set object.
#' @param reduction_method Reduced dimension to partitions clusters for.
#'
#' @examples
#'   \donttest{
#'     cell_metadata <- readRDS(system.file('extdata',
#'                                          'worm_embryo/worm_embryo_coldata.rds',
#'                                          package='monocle3'))
#'     gene_metadata <- readRDS(system.file('extdata',
#'                                          'worm_embryo/worm_embryo_rowdata.rds',
#'                                          package='monocle3'))
#'     expression_matrix <- readRDS(system.file('extdata',
#'                                              'worm_embryo/worm_embryo_expression_matrix.rds',
#'                                              package='monocle3'))
#'     cds <- new_cell_data_set(expression_data=expression_matrix,
#'                              cell_metadata=cell_metadata,
#'                              gene_metadata=gene_metadata)
#'
#'     cds <- preprocess_cds(cds)
#'     cds <- reduce_dimension(cds)
#'     cds <- cluster_cells(cds)
#'     partitions_factors <- partitions(cds, "UMAP")
#'   }
#'
#' @return Partitions.
#'
#' @export
setGeneric("partitions", function(x, reduction_method = "UMAP")
  standardGeneric("partitions"))

#' Method to extract partitions from CDS object
#' @param x A cell_data_set object.
#' @param reduction_method Reduced dimension to partitions clusters for.
#'
#' @return Partitions.
#'
#' @export
setMethod("partitions", "cell_data_set",
          function(x, reduction_method = "UMAP") {
            value <- x@clusters[[
              reduction_method]]$partitions[colnames(exprs(x))]
            if (is.null(value)) {
              stop("No partitions calculated for reduction_method = ",
                   reduction_method, ". Please first run ",
                   "cluster_cells with reduction_method = ",
                   reduction_method, ".")
            }
            return(value)
          })

#' Generic to extract principal graph from CDS
#' @param x A cell_data_set object.
#'
#' @examples
#'  \donttest{
#'     cell_metadata <- readRDS(system.file('extdata',
#'                                          'worm_embryo/worm_embryo_coldata.rds',
#'                                          package='monocle3'))
#'     gene_metadata <- readRDS(system.file('extdata',
#'                                          'worm_embryo/worm_embryo_rowdata.rds',
#'                                          package='monocle3'))
#'     expression_matrix <- readRDS(system.file('extdata',
#'                                              'worm_embryo/worm_embryo_expression_matrix.rds',
#'                                              package='monocle3'))
#'
#'     cds <- new_cell_data_set(expression_data=expression_matrix,
#'                              cell_metadata=cell_metadata,
#'                              gene_metadata=gene_metadata)
#'
#'     cds <- preprocess_cds(cds)
#'     cds <- align_cds(cds, alignment_group =
#'                      "batch", residual_model_formula_str = "~ bg.300.loading +
#'                       bg.400.loading + bg.500.1.loading + bg.500.2.loading +
#'                       bg.r17.loading + bg.b01.loading + bg.b02.loading")
#'    cds <- reduce_dimension(cds)
#'    ciliated_genes <- c("che-1", "hlh-17", "nhr-6", "dmd-6", "ceh-36", "ham-1")
#'    cds <- cluster_cells(cds)
#'    cds <- learn_graph(cds)
#'    pr_gr <- principal_graph(cds)
#'  }
#'
#' @return Principle graph.
#'
#' @export
setGeneric("principal_graph", function(x) standardGeneric("principal_graph"))

#' Generic to set principal graph to CDS
#' @param x A cell_data_set object.
#' @param value A principal graph object.
#'
#' @return x.
#'
#' @examples
#'  \donttest{
#'     cell_metadata <- readRDS(system.file('extdata',
#'                                          'worm_embryo/worm_embryo_coldata.rds',
#'                                          package='monocle3'))
#'     gene_metadata <- readRDS(system.file('extdata',
#'                                          'worm_embryo/worm_embryo_rowdata.rds',
#'                                          package='monocle3'))
#'     expression_matrix <- readRDS(system.file('extdata',
#'                                              'worm_embryo/worm_embryo_expression_matrix.rds',
#'                                              package='monocle3'))
#'
#'     cds <- new_cell_data_set(expression_data=expression_matrix,
#'                              cell_metadata=cell_metadata,
#'                              gene_metadata=gene_metadata)
#'
#'     cds <- preprocess_cds(cds)
#'     cds <- align_cds(cds, alignment_group =
#'                      "batch", residual_model_formula_str = "~ bg.300.loading +
#'                       bg.400.loading + bg.500.1.loading + bg.500.2.loading +
#'                       bg.r17.loading + bg.b01.loading + bg.b02.loading")
#'    cds <- reduce_dimension(cds)
#'    ciliated_genes <- c("che-1", "hlh-17", "nhr-6", "dmd-6", "ceh-36", "ham-1")
#'    cds <- cluster_cells(cds)
#'    cds <- learn_graph(cds)
#'    pr_gr <- principal_graph(cds)
#'    principal_graph(cds) <- NULL
#'    principal_graph(cds) <- pr_gr
#'  }
#'
#' @export
setGeneric("principal_graph<-", function(x, value)
  standardGeneric("principal_graph<-"))

#' Method to extract principal graph from CDS
#' @param x A cell_data_set object.
#'
#' @return Principle graph.
#'
#' @export
setMethod("principal_graph", "cell_data_set", function(x) {
  value <- x@principal_graph
  return(value)
})

#' Generic to set principal graph to CDS
#' @param x A cell_data_set object.
#' @param value A principal graph object.
#'
#' @examples
#'  \donttest{
#'     cell_metadata <- readRDS(system.file('extdata',
#'                                          'worm_embryo/worm_embryo_coldata.rds',
#'                                          package='monocle3'))
#'     gene_metadata <- readRDS(system.file('extdata',
#'                                          'worm_embryo/worm_embryo_rowdata.rds',
#'                                          package='monocle3'))
#'     expression_matrix <- readRDS(system.file('extdata',
#'                                              'worm_embryo/worm_embryo_expression_matrix.rds',
#'                                              package='monocle3'))
#'
#'     cds <- new_cell_data_set(expression_data=expression_matrix,
#'                              cell_metadata=cell_metadata,
#'                              gene_metadata=gene_metadata)
#'
#'     cds <- preprocess_cds(cds)
#'     cds <- align_cds(cds, alignment_group =
#'                      "batch", residual_model_formula_str = "~ bg.300.loading +
#'                       bg.400.loading + bg.500.1.loading + bg.500.2.loading +
#'                       bg.r17.loading + bg.b01.loading + bg.b02.loading")
#'    cds <- reduce_dimension(cds)
#'    ciliated_genes <- c("che-1", "hlh-17", "nhr-6", "dmd-6", "ceh-36", "ham-1")
#'    cds <- cluster_cells(cds)
#'    cds <- learn_graph(cds)
#'    pr_gr <- principal_graph(cds)
#'    principal_graph(cds) <- NULL
#'    principal_graph(cds) <- pr_gr
#'  }
#'
#' @return x.
#'
#' @export
#' @importClassesFrom S4Vectors List
setReplaceMethod("principal_graph", "cell_data_set", function(x, value) {
  value <- methods::as(value, methods::className("List","S4Vectors"))
  if (is.null(names(value))) {
    names(value) <- character(length(value))
  }
  x@principal_graph <- value
  methods::validObject(x)
  return(x)
})

#' Generic to extract principal graph auxiliary information from CDS
#' @param x A cell_data_set object.
#'
#' @examples
#'  \donttest{
#'     cell_metadata <- readRDS(system.file('extdata',
#'                                          'worm_embryo/worm_embryo_coldata.rds',
#'                                          package='monocle3'))
#'     gene_metadata <- readRDS(system.file('extdata',
#'                                          'worm_embryo/worm_embryo_rowdata.rds',
#'                                          package='monocle3'))
#'     expression_matrix <- readRDS(system.file('extdata',
#'                                              'worm_embryo/worm_embryo_expression_matrix.rds',
#'                                              package='monocle3'))
#'
#'     cds <- new_cell_data_set(expression_data=expression_matrix,
#'                              cell_metadata=cell_metadata,
#'                              gene_metadata=gene_metadata)
#'
#'     cds <- preprocess_cds(cds)
#'     cds <- align_cds(cds, alignment_group =
#'                      "batch", residual_model_formula_str = "~ bg.300.loading +
#'                       bg.400.loading + bg.500.1.loading + bg.500.2.loading +
#'                       bg.r17.loading + bg.b01.loading + bg.b02.loading")
#'    cds <- reduce_dimension(cds)
#'    ciliated_genes <- c("che-1", "hlh-17", "nhr-6", "dmd-6", "ceh-36", "ham-1")
#'    cds <- cluster_cells(cds)
#'    cds <- learn_graph(cds)
#'    pr_gr_aux <- principal_graph_aux(cds)
#'  }
#'
#' @return Principal graph auxiliary information.
#'
#' @export
setGeneric("principal_graph_aux", function(x)
  standardGeneric("principal_graph_aux"))

#' Generic to set principal graph auxiliary information into CDS
#' @param x A cell_data_set object.
#' @param value A S4Vectors::SimpleList of principal graph auxiliary information.
#'
#' @examples
#'  \donttest{
#'     cell_metadata <- readRDS(system.file('extdata',
#'                                          'worm_embryo/worm_embryo_coldata.rds',
#'                                          package='monocle3'))
#'     gene_metadata <- readRDS(system.file('extdata',
#'                                          'worm_embryo/worm_embryo_rowdata.rds',
#'                                          package='monocle3'))
#'     expression_matrix <- readRDS(system.file('extdata',
#'                                              'worm_embryo/worm_embryo_expression_matrix.rds',
#'                                              package='monocle3'))
#'
#'     cds <- new_cell_data_set(expression_data=expression_matrix,
#'                              cell_metadata=cell_metadata,
#'                              gene_metadata=gene_metadata)
#'
#'     cds <- preprocess_cds(cds)
#'     cds <- align_cds(cds, alignment_group =
#'                      "batch", residual_model_formula_str = "~ bg.300.loading +
#'                       bg.400.loading + bg.500.1.loading + bg.500.2.loading +
#'                       bg.r17.loading + bg.b01.loading + bg.b02.loading")
#'    cds <- reduce_dimension(cds)
#'    ciliated_genes <- c("che-1", "hlh-17", "nhr-6", "dmd-6", "ceh-36", "ham-1")
#'    cds <- cluster_cells(cds)
#'    cds <- learn_graph(cds)
#'    pr_gr_aux <- principal_graph_aux(cds)
#'    principal_graph_aux(cds) <- NULL
#'    principal_graph_aux(cds) <- pr_gr_aux
#'  }
#'
#' @return x.
#'
#' @export
setGeneric("principal_graph_aux<-", function(x, value)
  standardGeneric("principal_graph_aux<-"))

#' Method to extract principal graph auxiliary information from CDS
#' @param x A cell_data_set object.
#'
#' @return Principal graph auxiliary information.
#'
#' @export
setMethod("principal_graph_aux", "cell_data_set", function(x) {
  value <- x@principal_graph_aux
  return(value)
})

#' Method to set principal graph auxiliary information into CDS
#' @param x A cell_data_set object.
#' @param value A S4Vectors::SimpleList of principal graph auxiliary information.
#' 
#' @examples
#'  \donttest{
#'     cell_metadata <- readRDS(system.file('extdata',
#'                                          'worm_embryo/worm_embryo_coldata.rds',
#'                                          package='monocle3'))
#'     gene_metadata <- readRDS(system.file('extdata',
#'                                          'worm_embryo/worm_embryo_rowdata.rds',
#'                                          package='monocle3'))
#'     expression_matrix <- readRDS(system.file('extdata',
#'                                              'worm_embryo/worm_embryo_expression_matrix.rds',
#'                                              package='monocle3'))
#'
#'     cds <- new_cell_data_set(expression_data=expression_matrix,
#'                              cell_metadata=cell_metadata,
#'                              gene_metadata=gene_metadata)
#'
#'     cds <- preprocess_cds(cds)
#'     cds <- align_cds(cds, alignment_group =
#'                      "batch", residual_model_formula_str = "~ bg.300.loading +
#'                       bg.400.loading + bg.500.1.loading + bg.500.2.loading +
#'                       bg.r17.loading + bg.b01.loading + bg.b02.loading")
#'    cds <- reduce_dimension(cds)
#'    ciliated_genes <- c("che-1", "hlh-17", "nhr-6", "dmd-6", "ceh-36", "ham-1")
#'    cds <- cluster_cells(cds)
#'    cds <- learn_graph(cds)
#'    pr_gr_aux <- principal_graph_aux(cds)
#'    principal_graph_aux(cds) <- NULL
#'    principal_graph_aux(cds) <- pr_gr_aux
#'  }
#'
#' @return x.
#'
#' @export
#' @importClassesFrom S4Vectors List
setReplaceMethod("principal_graph_aux", "cell_data_set", function(x, value) {
  value <- methods::as(value, methods::className("List","S4Vectors"))
  if (is.null(names(value))) {
    names(value) <- character(length(value))
  }
  x@principal_graph_aux <- value
  methods::validObject(x)
  return(x)
})


# Set of wrappers for easy transition from monocle.

#' Generic to access cds count matrix
#' @param x A cell_data_set object.
#'
#' @examples
#'  \donttest{
#'    cds <- load_a549()
#'    exprs(cds)
#'  }
#'
#' @return Count matrix.
#'
#' @export
setGeneric("exprs", function(x) standardGeneric("exprs"))

#' Method to access cds count matrix
#' @param x A cell_data_set object.
#'
#' @return Count matrix.
#'
#' @export
setMethod("exprs", "cell_data_set", function(x) {
  value <- assays(x)$counts
  return(value)
})

#' Generic to access cds colData table
#' @param x A cell_data_set object.
#'
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'     pData(cds)
#'   }
#'
#' @return colData.
#'
#' @export
setGeneric("pData", function(x) standardGeneric("pData"))

#' Generic to set cds colData table
#' @param x A cell_data_set object.
#' @param value A data frame to set to colData table.
#'
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'     pData(cds)[['row_index']] <- seq(nrow(pData(cds)))
#'   }
#'
#' @return x.
#'
#' @export
setGeneric("pData<-", function(x, value) standardGeneric("pData<-"))

#' Method to access cds colData table
#' @param x A cell_data_set object.
#'
#' @return colData.
#'
#' @export
setMethod("pData", "cell_data_set", function(x) {
  value <- colData(x)
  return(value)
})

#' Method to set cds colData table
#' @param x A cell_data_set object.
#' @param value A data frame to set to colData table.
#'
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'     pData(cds)[['row_index']] <- seq(nrow(pData(cds)))
#'   }
#'
#' @return x.
#'
#' @export
#' @importClassesFrom S4Vectors List
setReplaceMethod("pData", "cell_data_set", function(x, value) {
  colData(x) <- value
  methods::validObject(x)
  return(x)
})

#' Generic to access cds rowData table
#' @param x A cell_data_set object.
#'
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'     fData(cds)
#'   }
#'
#' @return rowData table.
#'
#' @export
setGeneric("fData", function(x) standardGeneric("fData"))

#' Generic to set cds rowData table
#' @param x A cell_data_set object.
#' @param value A data frame to set to colData table.
#'
#' @return x.
#'
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'     fData(cds)[['row_index']] <- seq(nrow(fData(cds)))
#'   }
#'
#' @export
setGeneric("fData<-", function(x, value) standardGeneric("fData<-"))

#' Method to access cds rowData table
#' @param x A cell_data_set object.
#'
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'     fData(cds)
#'   }
#'
#' @return rowData table.
#'
#' @export
setMethod("fData", "cell_data_set", function(x) {
  value <- rowData(x)
  return(value)
})

#' Method to set cds rowData table
#' @param x A cell_data_set object.
#' @param value A data frame to set to colData table.
#'
#' @return x.
#'
#' @export
#' @importClassesFrom S4Vectors List
setReplaceMethod("fData", "cell_data_set", function(x, value) {
  rowData(x) <- value
  methods::validObject(x)
  return(x)
})
