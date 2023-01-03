leiden_clustering_test <- function(data,
                              pd,
                              weight=NULL,
                              nn_index=NULL,
                              k=20,
                              nn_control=list(),
                              num_iter=2,
                              resolution_parameter=0.0001,
                              random_seed=NULL,
                              verbose=FALSE, ...) {
  extra_arguments <- list(...)
  if( 'partition_type' %in% names( extra_arguments ) )
    partition_type <- extra_arguments[['partition_type']]
  else
    partition_type <- 'CPMVertexPartition'
  if( 'initial_membership' %in% names( extra_arguments ) )
    initial_membership <- extra_arguments[['initial_membership']]
  else
    initial_membership <- NULL
  if( 'weights' %in% names( extra_arguments ) )
    edge_weights <- extra_arguments[['weights']]
  else
    edge_weights <- NULL
  if( 'node_sizes' %in% names( extra_arguments ) )
    node_sizes <- extra_arguments[['node_sizes']]
  else
    node_sizes <- NULL

  # Check input parameters.
  assertthat::assert_that(assertthat::is.count(k))

  # The following vertex partitions have no resolution parameter.
  if( partition_type %in% c('ModularityVertexPartition','SignificanceVertexPartition','SurpriseVertexPartition') )
  {
    resolution_parameter = NA
  }
  else if( is.null( resolution_parameter ) )
  {
    resolution_parameter = 0.0001
  }
  if( is.null( num_iter ) )
    num_iter = 2
  if( random_seed == 0L )
    random_seed = NULL
  cell_names <- row.names(pd)
  if(!identical(cell_names, row.names(pd)))
    stop("Phenotype and row name from the data don't match")

  graph_result <- cluster_cells_make_graph(data=data,
                                           weight=weight,
                                           cell_names=cell_names,
                                           nn_index,
                                           k=k,
                                           nn_control=nn_control,
                                           verbose=verbose)

  if(verbose)
    message("  Run leiden clustering ...")

  t_start <- Sys.time()

  if(verbose)
  {
    table_results <- data.frame(
      resolution_parameter = double(),
      quality              = double(),
      modularity           = double(),
      significance         = double(),
      number_clusters      = integer() )
  }

  best_modularity <- -1
  best_result <- NULL
  best_resolution_parameter <- 'No resolution'
  # These three vertex partition types have a resolution parameter
  # so scan parameter range, if given.
  if(length(resolution_parameter) < 1) warning("bad loop: length(resolution_parameter) < 1")
  for(i in 1:length(resolution_parameter)) {
    cur_resolution_parameter <- resolution_parameter[i]
    cluster_result <- leidenbase::leiden_find_partition( graph_result[['g']],
                                                         partition_type = partition_type,
                                                         initial_membership = initial_membership,
                                                         edge_weights = edge_weights,
                                                         node_sizes = node_sizes,
                                                         seed = random_seed,
                                                         resolution_parameter = cur_resolution_parameter,
                                                         num_iter = num_iter,
                                                         verbose = verbose )
    quality      <- cluster_result[['quality']]
    modularity   <- cluster_result[['modularity']]
    significance <- cluster_result[['significance']]

    if(verbose)
      table_results <- rbind( table_results, data.frame(
        resolution_parameter = cur_resolution_parameter,
        quality              = quality,
        modularity           = modularity,
        significance         = significance,
        cluster_count        = max(cluster_result[['membership']]) ) )
    if(verbose)
      message('    Current resolution is ', cur_resolution_parameter,
              '; Modularity is ', modularity,
              '; Quality is ', quality,
              '; Significance is ', significance,
              '; Number of clusters is ', max(cluster_result[['membership']]))
    if(modularity > best_modularity) {
      best_result <- cluster_result
      best_resolution_parameter <- cur_resolution_parameter
      best_modularity <- modularity
    }
    if ( is.null( best_result ) ) {
      best_result <- cluster_result
      best_resolution_parameter <- NULL
      best_modularity <- cluster_result[['modularity']]
    }
  }
  t_end <-
    Sys.time()

  if(verbose)
  {
    message('    Done. Run time: ', t_end - t_start, 's\n')
    message('  Clustering statistics')
    selected <- vector( mode='character',
                        length = length( resolution_parameter ) )
    if(length(resolution_parameter) < 1 ) warning("bad loop: length(resolution_parameter) < 1")
    for(irespar in 1:length(resolution_parameter))
    {
      if( identical( table_results[['resolution_parameter']][irespar],
                     best_resolution_parameter ) )
        selected[irespar] <- '*'
      else
        selected[irespar] <- ' '
    }
    print( cbind(' '=' ',table_results, selected ),row.names=FALSE )
    message()
    message('  Cell counts by cluster')
    membership<-best_result[['membership']]
    membership_frequency <- stats::aggregate(data.frame(cell_count = membership),
                                             list(cluster = membership),
                                             length)
    membership_frequency <- cbind(' '=' ',membership_frequency,
                                  cell_fraction=sprintf("%.3f",membership_frequency[['cell_count']]/sum(membership_frequency[['cell_count']])))
    print( membership_frequency,row.names=FALSE)
    message()
    message('  Maximal modularity is ', best_modularity,
            ' for resolution parameter ', best_resolution_parameter)
    message("\n  Run kNN based graph clustering DONE.\n  -Number of clusters: ",
            max(best_result[['membership']]))
  }

  if(igraph::vcount(graph_result[['g']]) < 3000) {
    coord <- NULL
    edge_links <- NULL
  } else {
    coord <- NULL
    edge_links <- NULL
  }

  igraph::V(graph_result[['g']])$names <- as.character(igraph::V(graph_result[['g']]))
  out_result <- list(membership = best_result[['membership']],
                     modularity = best_result[['modularity']] )
  names(out_result$membership) = cell_names

  return(list(g=graph_result[['g']],
              relations=graph_result[['relations']],
              distMatrix=graph_result[['distMatrix']],
              coord=coord,
              edge_links=edge_links,
              optim_res=out_result))
  return(graph_result)
}
