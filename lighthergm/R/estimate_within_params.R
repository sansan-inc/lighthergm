#' Estimate a within-block network model.
#' @importFrom parallel mclapply
#' @importFrom magrittr %>%
#' @importFrom ergm ergm
#' @importFrom foreach foreach %do% %dopar%
#' @param formula a within network formula
#' @param network a network object
#' @param z_memb block memberships for each node
#' @param number_cores The number of CPU cores to use.
#' @param verbose A logical or an integer: if this is TRUE/1,
#' the program will print out additional information about the progress of estimation and simulation.
#' @param seeds seed value (integer) for the random number generator
#' @param method_second_step If "MPLE" (the default), then the maximum pseudolikelihood estimator is returned.
#' If "MLE", then an approximate maximum likelihood estimator is returned.
#' @param offset_coef a vector of model parameters to be fixed when estimation.(i.e., not estimated).
#' @param ... Additional arguments, to be passed to lower-level functions
#' @export
estimate_within_params <-
  function(formula,
           network,
           z_memb,
           number_cores = 1,
           verbose = 1,
           seeds = NULL,
           method_second_step = c("MPLE", "MLE"),
           offset_coef = NULL,
           ...) {

    varargs <- list(...)
    # Store block structure in a tibble
    block_structure <-
      tibble::tibble(
        node_id = network::network.vertex.names(network),
        block = z_memb
      )

    # Get number of clusters
    n_cluster <- length(unique(z_memb))

    # Get an edgelist and vertex attributes from the network object using intergraph::asDF
    list_edgelist <- intergraph::asDF(network)

    # Get vertex attributes.
    vertex_attr <- list_edgelist$vertexes %>%
      dplyr::select(-intergraph_id)

    # Extract vertex names from the network
    vertex_name <- list_edgelist$vertexes %>%
      dplyr::select(intergraph_id, vertex.names)

    # Construct an edgelist whose source_id and target_ids correspond to the vertex names of the network.
    # If you use network::as.edgelist instead, the original vertex names are not kept. That causes a problem when you convert an edgelist into a network object.
    edgelist <-
      list_edgelist$edges %>%
      dplyr::select(V1, V2) %>%
      dplyr::left_join(., vertex_name, by = c("V1" = "intergraph_id")) %>%
      dplyr::rename(source_id = vertex.names) %>%
      dplyr::select(V2, source_id) %>%
      dplyr::left_join(., vertex_name, by = c("V2" = "intergraph_id")) %>%
      dplyr::rename(target_id = vertex.names) %>%
      dplyr::select(source_id, target_id) %>%
      tibble::tibble()

    # Get a sparse adjacency matrix for each block. Store them in a list.
    # This computation might be unstable when the network is large. Check this point later.
    # For Windows users
    if (Sys.info()[["sysname"]] == "Windows") {
      # Preparation for parallel computing using foreach
      cluster <- parallel::makeCluster(number_cores, type = "PSOCK")
      doParallel::registerDoParallel(cluster)
      # Start computation
      block_net <- foreach(k = 1:n_cluster) %dopar% {
        # Get a subgraph whose vertices belong to block k.
        subnet <- get_induced_subgraph(block_structure, edgelist, k)
        # Keep vertex ids
        vertex_id <- network::network.vertex.names(subnet)
        # Convert the subgraph into a sparse adjacency matrix.
        sub_net <- as_sparse_adj(subnet)
        # Make a block attribute.
        # The length of this list must be the same with the number of vertices of the subgraph.
        block_attr <- rep(k, length(which(z_memb == k)))
        # Return the objects as a list.
        return(list(net = sub_net, id = vertex_id, block = block_attr))
      }
      parallel::stopCluster(cluster)
    }
    # For non-Windows users
    else {
      block_net <- mclapply(1:n_cluster, function(k) {
        # Get a subgraph whose vertices belong to block k.
        subnet <- get_induced_subgraph(block_structure, edgelist, k)
        # Keep vertex ids
        vertex_id <- network::network.vertex.names(subnet)
        # Convert the subgraph into a sparse adjacency matrix.
        sub_net <- as_sparse_adj(subnet)
        # Make a block attribute.
        # The length of this list must be the same with the number of vertices of the subgraph.
        block_attr <- rep(k, length(which(z_memb == k)))
        # Return the objects as a list.
        list(net = sub_net, id = vertex_id, block = block_attr)
      },
      mc.cores = number_cores
      )
    }

    # Extract info on vertex id.
    vertex_id <- unlist(c(block_net %>%
      purrr::map(function(x) {
        x$id
      })))

    # Extract info on which node belongs to which block.
    block_attr <- unlist(c(block_net %>%
      purrr::map(function(x) {
        x$block
      })))

    # Create a sparse adjacency matrix that only considers within-block connections.
    block_net <- Matrix::bdiag(block_net %>%
      purrr::map(function(x) {
        x$net
      }))

    # Convert the sparse adjacency matrix into a network object, called "block_net".
    edgelist <-
      Matrix::summary(block_net) %>%
      dplyr::select(i, j) %>%
      as.matrix()

    attr(edgelist, "n") <- length(vertex_id)
    attr(edgelist, "vnames") <- vertex_id
    attr(edgelist, "directed") <- FALSE
    attr(edgelist, "bipartite") <- FALSE
    attr(edgelist, "loops") <- FALSE
    attr(edgelist, "class") <- c("matrix", "edgelist")

    block_net <- network::network(edgelist, matrix.type = "edgelist", directed = FALSE)

    # Attach vertex id
    network::network.vertex.names(block_net) <- vertex_id

    # Attach vertex attributes
    df_vertex_attr <-
      tibble::tibble(vertex.names = vertex_id) %>%
      dplyr::left_join(., vertex_attr, by = "vertex.names")

    # Remove non-vertex-attribute columns
    df_vertex_attr <-
      df_vertex_attr %>%
      dplyr::select(-vertex.names)

    ## This part could be written in a much better way without foreach?
    df_vertex_attr_colnames <- colnames(df_vertex_attr)
    foreach(i = 1:ncol(df_vertex_attr)) %do% {
      network::set.vertex.attribute(x = block_net, attrname = df_vertex_attr_colnames[i], value = df_vertex_attr[[i]])
    }

    # Attach block attributes
    network::set.vertex.attribute(block_net, "block", block_attr)

    # Re-arrange the formula in such a way that the LHS is block_net.
    formula_terms <- as.character(formula)[3]

    # If within-block parameters are fixed, need to include `offset` in the formula
    if (!is.null(offset_coef)) {
      second_step_rhs <- as.character(formula)[3]
      second_step_rhs <- unlist(stringr::str_split(string = second_step_rhs, pattern = " \\+ "))
      # Extract offset terms and wrap them by `offset()`
      offset_terms <- second_step_rhs[stringr::str_detect(string = second_step_rhs, pattern = "\"")]
      offset_terms <- glue::glue("offset({offset_terms})")
      # Extract within-block parameters to be estimated
      within_terms <- second_step_rhs[!stringr::str_detect(string = second_step_rhs, pattern = "\"")]
      # Combine all
      formula_terms <- stringr::str_c(c(within_terms, offset_terms), collapse = " + ")
    }

    formula <- as.formula(glue::glue("block_net ~ {formula_terms}"))

    # Estimate the within-block parameters
    # The default estimation method is "MPLE", but you can select "MLE" if you like.
    if(is.null(varargs$control)){
      control <- ergm::control.ergm()
    } else {
      control <- varargs$control
    }

    model_est <- ergm(
      formula = formula,
      constraints = ~ blockdiag("block"),
      estimate = method_second_step,
      offset.coef = offset_coef,
      control = control
    )

    # Remove unnecessary network objects
    model_est$newnetwork <- NULL

    return(model_est)
  }

# ------------------------------------------------------------------------
# -------------- Auxiliary functions -------------------------------------
# ------------------------------------------------------------------------

#' Get a sparse adjacency matrix from a network object
#' @param net a network object
as_sparse_adj <- function(net) {
  n_nodes <- as.integer(net$gal$n)
  net <- network::as.edgelist(net)
  net <- Matrix::sparseMatrix(
    i = net[, 1],
    j = net[, 2],
    x = 1,
    dims = c(n_nodes, n_nodes),
    symmetric = TRUE
  )
  return(net)
}


get_induced_subgraph <- function(block_structure, edgelist, searched_block) {
  # Get the relevant nodes for the searched block
  block_nodes <-
    block_structure %>%
    dplyr::filter(block == searched_block) %>%
    .$node_id

  # Get the number of nodes in the block
  n_nodes_in_block <- length(block_nodes)

  # Get only the edges where nodes on both sides share the searched block membership
  subgraph <-
    edgelist %>%
    dplyr::filter(source_id %in% block_nodes & target_id %in% block_nodes)

  # Keep only the ID attributes
  subgraph <-
    subgraph %>%
    dplyr::select(source_id, target_id) %>%
    as.matrix()

  attr(subgraph, "n") <- n_nodes_in_block
  attr(subgraph, "vnames") <- block_nodes
  attr(subgraph, "directed") <- FALSE
  attr(subgraph, "bipartite") <- FALSE
  attr(subgraph, "loops") <- FALSE
  attr(subgraph, "class") <- c("matrix", "edgelist")

  if (nrow(subgraph) == 0) {
    subgraph <- network::network.initialize(n = n_nodes_in_block, directed = FALSE)
    network::network.vertex.names(subgraph) <- block_nodes
    return(subgraph)
  } else {
    # Return the subnet
    subgraph <-
      network::network(subgraph, directed = FALSE, matrix.type = "edgelist")

    return(subgraph)
  }
}
