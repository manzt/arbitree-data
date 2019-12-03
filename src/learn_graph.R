learn_graph_ <- function(cds, medioids=NULL, k=NULL, use_partition = TRUE, close_loop = TRUE, learn_graph_control = NULL, verbose = FALSE) {
  reduction_method <- "UMAP"
  
  # Create defaults
  if (!is.null(learn_graph_control)) {
    assertthat::assert_that(methods::is(learn_graph_control, "list"))
    assertthat::assert_that(all(
      names(learn_graph_control) 
      %in% 
      c(
        "euclidean_distance_ratio",
        "geodesic_distance_ratio",
        "minimal_branch_len",
        "orthogonal_proj_tip",
        "prune_graph",
        "scale",
        "ncenter",
        "maxiter",
        "eps",
        "L1.gamma",
        "L1.sigma"
      )),
      msg = "Unknown variable in learn_graph_control"
    )
  }
  euclidean_distance_ratio <- ifelse(is.null(learn_graph_control$euclidean_distance_ratio), 1, learn_graph_control$euclidean_distance_ratio)
  geodesic_distance_ratio <- ifelse(is.null(learn_graph_control$geodesic_distance_ratio), 1/3, learn_graph_control$geodesic_distance_ratio)
  minimal_branch_len <- ifelse(is.null(learn_graph_control$minimal_branch_len), 10, learn_graph_control$minimal_branch_len)
  orthogonal_proj_tip <- ifelse(is.null(learn_graph_control$orthogonal_proj_tip), FALSE, learn_graph_control$orthogonal_proj_tip)
  prune_graph <- ifelse(is.null(learn_graph_control$prune_graph), TRUE, learn_graph_control$prune_graph)
  ncenter <- learn_graph_control$ncenter
  scale <- ifelse(is.null(learn_graph_control$scale), FALSE, learn_graph_control$scale)
  maxiter <- ifelse(is.null(learn_graph_control$maxiter), 10, learn_graph_control$maxiter)
  eps <- ifelse(is.null(learn_graph_control$eps), 1e-5, learn_graph_control$eps)
  L1.gamma <- ifelse(is.null(learn_graph_control$L1.gamma), 0.5, learn_graph_control$L1.gamma)
  L1.sigma <- ifelse(is.null(learn_graph_control$L1.sigma), 0.01, learn_graph_control$L1.sigma)
  
  # Check arguments
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(is.logical(use_partition))
  assertthat::assert_that(is.logical(close_loop))
  assertthat::assert_that(is.logical(verbose))
  assertthat::assert_that(is.logical(orthogonal_proj_tip))
  assertthat::assert_that(is.logical(prune_graph))
  assertthat::assert_that(is.logical(scale))
  assertthat::assert_that(is.numeric(euclidean_distance_ratio))
  assertthat::assert_that(is.numeric(geodesic_distance_ratio))
  assertthat::assert_that(is.numeric(minimal_branch_len))
  if(!is.null(ncenter)) assertthat::assert_that(assertthat::is.count(ncenter))
  assertthat::assert_that(assertthat::is.count(maxiter))
  assertthat::assert_that(is.numeric(eps))
  assertthat::assert_that(is.numeric(L1.sigma))
  assertthat::assert_that(is.numeric(L1.sigma))
  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste("No dimensionality reduction for",
                                      reduction_method, "calculated.",
                                      "Please run reduce_dimensions with",
                                      "reduction_method =", reduction_method,
                                      "and cluster_cells before running",
                                      "learn_graph."))
  assertthat::assert_that(!is.null(cds@clusters[[reduction_method]]),
                          msg = paste("No cell clusters for",
                                      reduction_method, "calculated.",
                                      "Please run cluster_cells with",
                                      "reduction_method =", reduction_method,
                                      "before running learn_graph."))
  
  if (use_partition) {
    partition_list <- cds@clusters[[reduction_method]]$partitions
  } else {
    partition_list <- rep(1, nrow(colData(cds)))
  }
  
  multi_tree_DDRTree_res <- multi_component_RGE(
    cds, 
    medioids = medioids,
    k_init = k,
    scale = scale,
    reduction_method = reduction_method,
    partition_list = partition_list,
    irlba_pca_res = reducedDims(cds)[[reduction_method]],
    max_components = max_components,
    ncenter = ncenter,
    maxiter = maxiter,
    eps = eps,
    L1.gamma = L1.gamma,
    L1.sigma = L1.sigma,
    close_loop = close_loop,
    euclidean_distance_ratio = euclidean_distance_ratio,
    geodesic_distance_ratio = geodesic_distance_ratio,
    prune_graph = prune_graph,
    minimal_branch_len = minimal_branch_len,
    verbose = verbose
  )
  
  rge_res_W <- multi_tree_DDRTree_res$ddrtree_res_W
  rge_res_Z <- multi_tree_DDRTree_res$ddrtree_res_Z
  rge_res_Y <- multi_tree_DDRTree_res$ddrtree_res_Y
  cds <- multi_tree_DDRTree_res$cds
  dp_mst <- multi_tree_DDRTree_res$dp_mst
  
  
  principal_graph(cds)[[reduction_method]] <- dp_mst
  cds@principal_graph_aux[[reduction_method]]$dp_mst <- rge_res_Y
  
  cds <- project2MST(cds, project_point_to_line_segment, orthogonal_proj_tip, verbose, reduction_method, rge_res_Y)
  
  cds
}