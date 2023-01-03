find_gene_modules_test <- function(cds,
                          reduction_method = c("UMAP"),
                          max_components = 2,
                          umap.metric = "cosine",
                          umap.min_dist = 0.1,
                          umap.n_neighbors = 15L,
                          umap.fast_sgd = FALSE,
                          umap.nn_method = "annoy",
                          k = 20,
                          leiden_iter = 1,
                          partition_qval = 0.05,
                          weight = FALSE,
                          resolution = NULL,
                          random_seed = 0L,
                          cores=1,
                          verbose = FALSE,
                          preprocess_method = c('PCA', 'LSI'),
                          nn_control = list(),
                          ...) {
  method = 'leiden'

  nn_control_default <- get_global_variable('nn_control_annoy_euclidean')
  nn_control <- set_nn_control(mode=3,
                               nn_control=nn_control,
                               nn_control_default=nn_control_default,
                               nn_index=NULL,
                               k=k,
                               verbose=verbose)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(preprocess_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "preprocess_method must be one of 'PCA' or 'LSI'")
  preprocess_method <- match.arg(preprocess_method)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'UMAP', 'PCA' or 'tSNE'")
  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(is.character(reduction_method))
  assertthat::assert_that(assertthat::is.count(k))
  assertthat::assert_that(is.logical(weight))
  assertthat::assert_that(assertthat::is.count(leiden_iter))
  ## TO DO what is resolution?
  assertthat::assert_that(is.numeric(partition_qval))
  assertthat::assert_that(is.logical(verbose))
  assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
                          msg = paste("No dimensionality reduction for",
                                      reduction_method, "calculated.",
                                      "Please run reduce_dimension with",
                                      "reduction_method =", reduction_method,
                                      "before running cluster_cells"))

  # preprocess_mat is gene_loading matrix. The gene_loadings were calculated only
  # for preprocess_method='PCA' in preprocess_cds() but I extend this to 'LSI' and
  # calculate gene_loadings here.
  preprocess_mat <- cds@reduce_dim_aux[[preprocess_method]][['model']]$svd_v %*% diag(cds@reduce_dim_aux[[preprocess_method]][['model']]$svd_sdev)

# Notes:
#   o  the beta vector is in cds@reduce_dim_aux[['Aligned']][['model']][['beta']]
#   o  cds@reduce_dim_aux[['Aligned']][['model']][['beta']] is npc x nfactor, which causes
#      preprocess_mat to have nfactor columns, often one column
#   o  I do not know how to adjust gene_loadings for batch effects
#      so this is disabled for now
#  if (!is.null(cds@reduce_dim_aux[['Aligned']][['model']][['beta']])){
#    preprocess_mat = preprocess_mat %*% (-cds@reduce_dim_aux[['Aligned']][['model']][['beta']])
#  }
  preprocess_mat <- preprocess_mat[intersect(rownames(cds), row.names(preprocess_mat)),]

  # uwot::umap uses a random number generator
  if( random_seed != 0L )
    set.seed( random_seed )

  umap_res = uwot::umap(as.matrix(preprocess_mat),
                        n_components = max_components,
                        metric = umap.metric,
                        min_dist = umap.min_dist,
                        n_neighbors = umap.n_neighbors,
                        fast_sgd = umap.fast_sgd,
                        n_threads=cores,
                        verbose=verbose,
                        nn_method= umap.nn_method,
                        ...)

  row.names(umap_res) <- row.names(preprocess_mat)
  if(ncol(umap_res) < 1) warning('bad loop: ncol(umap_res) < 1')
  colnames(umap_res) <- paste0('dim_', 1:ncol(umap_res))
  reduced_dim_res <- umap_res

  if(verbose)
    message("Running leiden clustering algorithm ...")

  cluster_result <- leiden_clustering(data=reduced_dim_res,
                                      pd=rowData(cds)[row.names(reduced_dim_res),,drop=FALSE],
                                      weight=weight,
                                      nn_index=NULL,
                                      k=k,
                                      nn_control=nn_control,
                                      num_iter=leiden_iter,
                                      resolution_parameter=resolution,
                                      random_seed=random_seed,
                                      verbose=verbose,
                                      ...)

  cluster_graph_res <- compute_partitions(cluster_result$g,
                                          cluster_result$optim_res,
                                          partition_qval, verbose)
  partitions <-
    igraph::components(cluster_graph_res$cluster_g)$membership[
      cluster_result$optim_res$membership]
  names(partitions) <- row.names(reduced_dim_res)
  partitions <- as.factor(partitions)

  gene_module_df <- tibble::tibble(id = row.names(preprocess_mat),
                                   module = factor(
                                     igraph::membership(cluster_result$optim_res)),
                                   supermodule = partitions)
  gene_module_df <- tibble::as_tibble(cbind(gene_module_df, umap_res))

  return(gene_module_df)
}
