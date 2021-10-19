#' Simulate a network
#' @param formula_for_simulation formula for simulating a network
#' @param data_for_simulation a data frame that contains vertex id, block membership, and vertex features.
#' @param colname_vertex_id a column name in the data frame for the vertex id
#' @param seed_edgelist an edgelist used for creating a seed network. It should have the "edgelist" class
#' @param colname_block_membership a column name in the data frame for the block membership
#' @param coef_within_block a vector of within-block parameters. The order of the parameters should match that of the formula.
#' @param coef_between_block a vector of between-block parameters. The order of the parameters should match that of the formula without externality terms.
#' @param ergm_control auxiliary function as user interface for fine-tuning ERGM simulation
#' @param seed seed value (integer) for network simulation.
#' @param directed whether the simulated network is directed
#' @param n_sim number of networks generated
#' @param output Normally character, one of "network" (default), "stats", "edgelist", to determine the output format.
#' @param prevent_duplicate If `TRUE`, the coefficient on nodematch("block") is set to be a very large negative number in drawing between-block links, so that there will be (almost) no within-block links.
#' @param use_fast_between_simulation If `TRUE`, this function uses an effcient way to simulate a between-block network. If the network is very large, you should consider using this option.
#' Note that when you use this, the first element of `coef_between_block` must be the edges parameter.
#' @param list_feature_matrices a list of feature adjacency matrices. If `use_fast_between_simulation`, this must be given.
#' @param verbose If this is TRUE/1, the program will print out additional information about the progress of simulation.
#' @param ... Additional arguments, to be passed to lower-level functions
#'
#' @examples
#' # Load an embedded network object.
#' data(toyNet)
#'
#' # Specify the model that you would like to estimate.
#' model_formula <- toyNet ~ edges + nodematch("x") + nodematch("y") + triangle
#'
#' # Estimate the model
#' hergm_res <-
#'   lighthergm::hergm(
#'     object = model_formula, # The model you would like to estiamte
#'     n_clusters = 4, # The number of blocks
#'     n_em_step_max = 100, # The maximum number of EM algorithm steps
#'     estimate_parameters = TRUE, # Perform parameter estimation after the block recovery step
#'     clustering_with_features = TRUE, # Indicate that clustering must take into account nodematch on characteristics
#'     check_block_membership = TRUE # Keep track of block memberships at each EM iteration
#'   )
#'
#' # Prepare a data frame that contains nodal id and covariates.
#' nodes_data <-
#'   data.frame(
#'     node_id = network::network.vertex.names(toyNet),
#'     block = hergm_res$partition,
#'     x = network::get.vertex.attribute(toyNet, "x"),
#'     y = network::get.vertex.attribute(toyNet, "y")
#'   )
#' # The feature adjacency matrices
#' list_feature_matrices <- lighthergm::get_list_sparse_feature_adjmat(toyNet, model_formula)
#'
#' # Estimated coefficients for the between-community connections
#' coef_between_block <- coef(hergm_res$est_between)
#'
#' # Estimated coefficients for the within-community connections
#' coef_within_block <- coef(hergm_res$est_within)
#'
#' # The MCMC settings
#' sim_ergm_control <- ergm::control.simulate.formula(
#'   MCMC.burnin = 1000000,
#'   MCMC.interval = 100000
#' )
#'
#' # Simulate network stats
#' sim_stats <- lighthergm::simulate_hergm(
#'   formula_for_simulation = model_formula, # Formula for between-blocks
#'   data_for_simulation = nodes_data, # Same as for gof, a dataframe containing nodes attributes
#'   colname_vertex_id = "node_id", # Name of the column containing node IDs
#'   colname_block_membership = "block", # Name of the column containing block IDs
#'   coef_between_block = coef_between_block, # The coefficients for the between connections
#'   coef_within_block = coef_within_block, # The coefficients for the within connections
#'   ergm_control = sim_ergm_control,
#'   n_sim = 100, # Number of simulations to return
#'   output = "stats", # If `stats` a list with network statistics for the between and within connections is returned
#'   use_fast_between_simulation = TRUE, # Simulates between connections by drawing from a logistic distribution. If FALSE, draws between connections by MCMC.
#'   list_feature_matrices = list_feature_matrices
#' )
#' @export
simulate_hergm <- function(formula_for_simulation,
                           data_for_simulation,
                           colname_vertex_id,
                           colname_block_membership,
                           seed_edgelist = NULL,
                           coef_within_block,
                           coef_between_block,
                           ergm_control = ergm::control.simulate.formula(),
                           seed = NULL,
                           directed = FALSE,
                           n_sim = 1,
                           output = "network",
                           prevent_duplicate = TRUE,
                           use_fast_between_simulation = FALSE,
                           list_feature_matrices = NULL,
                           verbose = 0,
                           ...) {

  # The current hergm doesn't support directed networks.
  if (directed) {
    stop("\nThe current hergm doesn't support directed networks. This will be modified in the future.")
  }

  # When use_fast_between_simulation:
  if (use_fast_between_simulation) {
    # Check if list_feature_matrices is given. If not, return an error.
    if (is.null(list_feature_matrices)) {
      stop("\nuse_fast_between_simulation = TRUE, you also need to give list_feature_matrices. \nYou can compute one using get_list_sparse_feature_adjmat().")
    }
  }

  if (verbose > 0) {
    message("Sorting the given dataset by the block ids.")
  }
  # Create a data frame from block memberships, vertex ids, and vertex covariates, sorted by block ids.
  sorted_dataframe <- sort_block_membership(data_for_simulation, colname_vertex_id, colname_block_membership)

  # Remove unnecessary data frame
  rm(data_for_simulation)

  # If seed_edgelist is given
  if (!is.null(seed_edgelist)) {
    # check if its class is "edgelist". It must have an edgelist, the number of vertices, and vertex ids.
    if (!any(class(seed_edgelist) == "edgelist")) {
      stop("\nThe given edgelist must have the 'edgelist' class, i.e., it must have an edgelist, the number of vertices, and vertex names.")
    }
    # Arrange the given edgelist
    list_edgelist <- arrange_edgelist(seed_edgelist, sorted_dataframe)
    seed_edgelist_within <- list_edgelist$edgelist_within
    seed_edgelist_between <- list_edgelist$edgelist_between
  } else {
    seed_edgelist_within <- NULL
    seed_edgelist_between <- NULL
  }

  # Create seed networks for from which networks will be simulated.
  ## For within-block links
  if (verbose > 0) {
    message("Creating a seed within-block network.")
  }
  seed_network_within <-
    generate_seed_network(formula_for_simulation,
      sorted_dataframe,
      edgelist = seed_edgelist_within,
      directed
    )

  # If output == "network":
  if (!is.function(output) && output == "network") {
    temp <- "edgelist"
  } else {
    temp <- output
  }

  # Create formula for simulating within- and between-block networks.
  ## Extract the RHS of the given formula
  formula_rhs <- as.character(formula_for_simulation)[3]
  ## Create a formula for simulation
  ### For within-block network
  formula_for_simulation_within <- as.formula(glue::glue("seed_network_within ~ {formula_rhs}"))
  ### For between-block network
  ### Here, the LHS network object can be arbitrary as long as the RHS is correct.
  formula_for_simulation_between <- separate_formulas(formula_for_simulation_within)[[1]]

  # Send a message about which value of a coefficient is attached to which term.
  # This is to make sure that `coef_within_block` and `coef_between_block` are correctly specified as the user intends.
  attach_terms_to_coef <- function(formula, coef) {
    terms <- as.character(formula)[3]
    terms <- unlist(stringr::str_split(string = terms, pattern = " \\+ "))
    terms <- stringr::str_replace_all(terms, "\"", "'")
    names(coef) <- terms
    return(coef)
  }
  coefs <- purrr::map2(
    list(formula_for_simulation_within, formula_for_simulation_between),
    list(coef_within_block, coef_between_block),
    attach_terms_to_coef
  )

  # Print out the specified coefficients
  if (verbose > 0) {
    message("Specified coefficients for the within-block model:")
    print(coefs[[1]])

    message("Specified coefficients for the between-block model:")
    print(coefs[[2]])
  }

  ############################################################################
  ############# Draw within-block connections ################################
  ############################################################################

  if (verbose > 0) {
    message("Simulating within-block networks.")
  }
  edgelist_within <- draw_within_block_connection(
    seed_network = seed_network_within,
    formula_for_simulation = formula_for_simulation_within,
    coef_within_block = coef_within_block,
    ergm_control = ergm_control,
    output = temp,
    seed = seed,
    n_sim = n_sim,
    verbose = verbose
  )

  ############################################################################
  ############# Draw between-block connections ###############################
  ############################################################################

  if (verbose > 0) {
    message("Simulating between-block networks.")
  }
  edgelist_between <-
    draw_between_block_connection(
      formula_for_simulation = formula_for_simulation_between,
      sorted_dataframe = sorted_dataframe,
      coef_between_block = coef_between_block,
      seed_edgelist_between = seed_edgelist_between,
      list_feature_matrices = list_feature_matrices,
      use_fast_between_simulation = use_fast_between_simulation,
      seed = seed,
      n_sim = n_sim,
      prevent_duplicate = prevent_duplicate,
      verbose = verbose,
      ergm_control = ergm_control,
      output = temp
    )

  if (!is.function(output) && output == "network") {
    if (n_sim == 1) {
      # Combine within- and between-block edges while removing duplicated links.
      edgelist <- combine_within_between_edges(edgelist_within, edgelist_between, use_fast_between_simulation)
      # Create a final network
      output <- generate_network_for_output(formula_for_simulation_within, formula_for_simulation_between, sorted_dataframe, edgelist)
      # Return the output
      if (verbose > 0) {
        message("One entire network has been generated.")
      }
      return(output)
    }
    if (n_sim > 1) {
      # Combine within- and between-block edges while removing duplicated links.
      fn_combine <- function(edgelist_within, edgelist_between) {
        combine_within_between_edges(edgelist_within, edgelist_between, use_fast_between_simulation)
      }
      edgelist <- purrr::map2(edgelist_within, edgelist_between, fn_combine)
      # Create a final network
      fn_final_network <- function(edgelist) {
        return(generate_network_for_output(formula_for_simulation_within, formula_for_simulation_between, sorted_dataframe, edgelist))
      }
      output <- purrr::map(edgelist, fn_final_network)
      # Return the output
      if (verbose > 0) {
        message(glue::glue("{n_sim} entire networks have been generated."))
      }
      return(output)
    }
  }

  if (output == "stats") {
    if (use_fast_between_simulation) {
      output_between <- get_between_stats(edgelist_between, between_formula = formula_for_simulation_between)
    } else {
      output_between <- data.frame(edgelist_between)
    }
    output_within <- data.frame(edgelist_within)
  } else {
    output_between <- edgelist_between
    output_within <- edgelist_within
  }

  # When output != "edgelist", like "stat" to see whether MCMC is converging.
  return(list(within_network = output_within, between_network = output_between))
}

# ----------------------------------------------------------------------------------------------
# Auxiliary functions for simulating networks --------------------------------------------------
# ----------------------------------------------------------------------------------------------

#' Create a data frame from block memberships, vertex ids, and vertex covariates, sorted by block ids.
#' @param data_for_simulation a data frame that contains vertex id, block membership, and vertex features.
#' @param colname_vertex_id a column name in the data frame for the vertex id
#' @param colname_block_membership a column name in the data frame for the block membership.
#' @importFrom magrittr %>%
sort_block_membership <- function(data_for_simulation, colname_vertex_id, colname_block_membership) {
  # Rename column names of vertex id and block membership.
  colname_vertex_id <- dplyr::enquo(colname_vertex_id)
  colname_block_membership <- dplyr::enquo(colname_block_membership)
  output <-
    data_for_simulation %>%
    dplyr::rename(
      vertex_id = !!colname_vertex_id,
      block = !!colname_block_membership
    )

  # Sort block memberships, vertex ids, and vertex features.
  output <-
    output %>%
    dplyr::mutate(block = as.numeric(block)) %>%
    dplyr::arrange(block)

  # If output$vertex_id is factor, convert it into character
  if (is.factor(output$vertex_id)) {
    output$vertex_id <- as.character(output$vertex_id)
  }

  # Return the output
  return(tibble::as_tibble(output))
}


#' Make the given edgelist consistent with the data frame that contains vertex info.
#' @param edgelist an edgelist to be arranged
#' @param sorted_dataframe a data frame sorted by `sort_block_membership`
arrange_edgelist <- function(edgelist, sorted_dataframe) {
  # Number of vertices
  N_node <- attr(edgelist, "n")
  # Vertex names
  vnames <- attr(edgelist, "vnames")
  # Whether the network is directed
  directed <- attr(edgelist, "directed")
  bipartite <- attr(edgelist, "bipartite")
  loops <- attr(edgelist, "loops")
  class <- attr(edgelist, "class")

  # If the order of vnames matches that of sorted_dataframe$vertex_id, we just need to attach block memberships to the edgelist and split it.
  # Otherwise, the following steps are necessary
  if (!all(vnames == sorted_dataframe$vertex_id)) {
    # Convert the given edgelist into a data frame
    df_edgelist <- data.frame(edgelist)
    names(df_edgelist) <- c("tail_old", "head_old")
    # Create a data frame that stores pairs of vnames and serial ids
    df_vnames_id <-
      data.frame(
        id_old = 1:N_node,
        vnames = vnames
      )
    # Attach vnames to the edgelist
    df_edgelist <-
      dplyr::left_join(df_edgelist, df_vnames_id, by = c("tail_old" = "id_old")) %>%
      dplyr::rename(tail_vnames = vnames) %>%
      dplyr::select(tail_old, head_old, tail_vnames) %>%
      dplyr::left_join(., df_vnames_id, by = c("head_old" = "id_old")) %>%
      dplyr::rename(head_vnames = vnames) %>%
      dplyr::select(tail_old, head_old, tail_vnames, head_vnames)

    # Create a data frame that consists of sorted vertex names, their new vertex serial numbers, and block memberships,
    sorted_vnames_block <-
      data.frame(
        vnames_sorted = sorted_dataframe$vertex_id,
        vertex_serial_id = 1:N_node,
        block = sorted_dataframe$block
      )
    # Merge datasets to attach new vertex serial numbers
    output <-
      dplyr::left_join(df_edgelist, sorted_vnames_block, by = c("tail_vnames" = "vnames_sorted")) %>%
      dplyr::rename(
        tail = vertex_serial_id,
        tail_block = block
      ) %>%
      dplyr::select(tail, head_old, head_vnames, tail_block) %>%
      dplyr::left_join(., sorted_vnames_block, by = c("head_vnames" = "vnames_sorted")) %>%
      dplyr::rename(
        head = vertex_serial_id,
        head_block = block
      ) %>%
      dplyr::select(tail, head, tail_block, head_block) %>%
      dplyr::arrange(tail_block, head_block)
  } else {
    # When the order of vnames matches that of sorted_dataframe$vertex_id:
    # Convert the given edgelist into a data frame
    df_edgelist <- data.frame(edgelist)
    names(df_edgelist) <- c("tail", "head")
    df_block <-
      data.frame(
        id = 1:N_node,
        block = sorted_dataframe$block
      )
    # Attach block memberships to the edgelist
    output <-
      dplyr::left_join(df_edgelist, df_block, by = c("tail" = "id")) %>%
      dplyr::rename(tail_block = block) %>%
      dplyr::select(tail, head, tail_block) %>%
      dplyr::left_join(., df_block, by = c("head" = "id")) %>%
      dplyr::rename(head_block = block) %>%
      dplyr::select(tail, head, tail_block, head_block)
  }

  # Get within- and between-block edgelists.
  ## Within-block edgelist
  edgelist_within <-
    output %>%
    dplyr::filter(tail_block == head_block) %>%
    as.matrix()
  attr(edgelist_within, "n") <- N_node
  attr(edgelist_within, "vnames") <- sorted_dataframe$vertex_id
  attr(edgelist_within, "directed") <- directed
  attr(edgelist_within, "bipartite") <- bipartite
  attr(edgelist_within, "loops") <- loops
  attr(edgelist_within, "class") <- class

  ## Between-block edgelist
  edgelist_between <-
    output %>%
    dplyr::filter(tail_block != head_block) %>%
    as.matrix()
  attr(edgelist_between, "n") <- N_node
  attr(edgelist_between, "vnames") <- sorted_dataframe$vertex_id
  attr(edgelist_between, "directed") <- directed
  attr(edgelist_between, "bipartite") <- bipartite
  attr(edgelist_between, "loops") <- loops
  attr(edgelist_between, "class") <- class

  # Return the output
  return(list(
    edgelist_within = edgelist_within,
    edgelist_between = edgelist_between
  ))
}


#' Create a seed network from which a network will be simulated.
#' @param formula_for_simulation formula for simulating a network
#' @param sorted_dataframe a data frame generated by `sort_block_membership`
#' @param edgelist an edgelist from which a seed network is generated. The class of the edgelist should be "edgelist", i.e. it should contain as attributes the number of nodes and vertex names.
#' @param directed a boolean of whether the network is directed.
#' @importFrom foreach foreach %do%
generate_seed_network <- function(formula_for_simulation, sorted_dataframe, edgelist = NULL, directed) {
  # Number of nodes
  N_node <- nrow(sorted_dataframe)

  # Initialize a network object
  if (is.null(edgelist)) {
    # If no edgelist is given:
    g <- network::network.initialize(n = N_node, directed = directed)
  } else {
    # If an edgelist is given:
    directed <- attr(edgelist, "directed")
    g <- network::network(edgelist, matrix.type = "edgelist", directed = directed)
  }

  # Get variable names from formula (extract strings sandwiched by double quotes)
  list_varname <- as.character(formula_for_simulation)[3]
  list_varname <- unlist(stringr::str_extract_all(string = list_varname, pattern = '"[^"]*"'))
  list_varname <- stringr::str_remove_all(string = list_varname, pattern = '\"')

  # Attach vertex features
  foreach(i = 1:length(list_varname)) %do% {
    feature <-
      dplyr::pull(sorted_dataframe, list_varname[i])
    # Network objects can't accept factors.
    if (class(feature) == "factor") {
      feature <- as.character(feature)
    }
    network::set.vertex.attribute(g, attrname = list_varname[i], value = feature)
  }

  # Attach block info
  block <- dplyr::pull(sorted_dataframe, "block")
  # Network objects can't accept factors.
  if (class(block) == "factor") {
    block <- as.character(block)
  }
  network::set.vertex.attribute(g, attrname = "block", value = block)

  # Attach vertex id
  # When the seed network is generated from an edgelist, this step is not necessary.
  if (is.null(edgelist)) {
    vertex_id <- dplyr::pull(sorted_dataframe, "vertex_id")
    # Network objects can't accept factors.
    if (class(vertex_id) == "factor") {
      vertex_id <- as.character(vertex_id)
    }
    network::set.vertex.attribute(g, attrname = "vertex.names", value = vertex_id)
  }

  # Return the seed network
  return(g)
}


#' Draw within-block connections
#' @param seed_network a seed network from which a network will be simulated.
#' @param formula_for_simulation formula for simulating a network
#' @param coef_within_block a vector of within-block parameters. The order of the parameters should match that of the formula.
#' @param ergm_control auxiliary function as user interface for fine-tuning ERGM simulation
#' @param output Normally character, one of "network" (default), "stats", "edgelist", to determine the output format.
#' @param seed seed value (integer) for the random number generator.
#' @param n_sim Number of networks to be randomly drawn from the given distribution on the set of all networks.
#' @param verbose If this is TRUE/1, the program will print out additionalinformation about the progress of simulation.
#' @param ... Additional arguments, to be passed to lower-level functions
draw_within_block_connection <- function(seed_network,
                                         formula_for_simulation,
                                         coef_within_block,
                                         ergm_control,
                                         output,
                                         seed,
                                         n_sim,
                                         verbose,
                                         ...) {
  # Extract the RHS of the formula
  formula_rhs <- as.character(formula_for_simulation)[3]
  # Create a formula for simulation
  sim_formula <- as.formula(glue::glue("seed_network ~ {formula_rhs}"))
  # Simulate within-block links
  within_conn <- ergm::simulate_formula(
    object = sim_formula,
    nsim = n_sim,
    coef = coef_within_block,
    constraints = ~ blockdiag("block"),
    seed = seed,
    control = ergm_control,
    output = output,
    verbose = verbose
  )
  # If output == "edgelist", attach the class c("edgelist", "matrix")
  if (output == "edgelist") {
    attr(within_conn, "vnames") <- network::network.vertex.names(seed_network)
    class(within_conn) <- c("matrix_edgelist", "edgelist", class(within_conn))
  }
  # Return the output
  return(within_conn)
}



#' Draw between-block connections. There may be some edges that appear both in within- and between-block links.
#' The overlapped edges will be removed after this step.
#' @param formula_for_simulation formula for simulating a between-block network
#' @param sorted_dataframe a data frame generated by `sort_block_membership`
#' @param coef_between_block a vector of between-block parameters. The order of the parameters should match that of the formula.
#' @param seed_edgelist_between a seed edgelist from which a between-block network is simulated.
#' @param use_fast_between_simulation If `TRUE`, this function uses an effcient way to simulate a between-block network.
#' If the network is very large, you should consider using this option.
#' Note that when you use this, the first element of `coef_between_block` must be the edges parameter.
#' @param list_feature_matrices a list of feature adjacency matrices. This is used when `use_fast_between_simulation`.
#' @param seed seed value (integer) for the random number generator.
#' @param n_sim number of networks generated.
#' @param prevent_duplicate If `TRUE`, the coefficient on nodematch("block") is set to be a very large negative number in drawing between-block links,
#' so that there will be (almost) no within-block links.
#' @param verbose If this is TRUE/1, the program will print out additionalinformation about the progress of simulation.
#' @param ergm_control auxiliary function as user interface for fine-tuning ERGM simulation
#' @param output Normally character, one of "network" (default), "stats", "edgelist", to determine the output format.
#' @param ... Additional arguments, to be passed to lower-level functions
draw_between_block_connection <- function(formula_for_simulation,
                                          sorted_dataframe,
                                          coef_between_block,
                                          seed_edgelist_between = NULL,
                                          use_fast_between_simulation = FALSE,
                                          list_feature_matrices = NULL,
                                          seed = NULL,
                                          n_sim = 1,
                                          prevent_duplicate = TRUE,
                                          verbose = 0,
                                          ergm_control = ergm::control.simulate.formula(),
                                          output = "edgelist",
                                          ...) {
  # If use_fast_between_simulation, simulate between-block networks using the cpp function `simulate_between_network()`.
  if (use_fast_between_simulation) {
    # Set seed
    set.seed(seed)
    # Number of nodes
    N_node <- length(sorted_dataframe$vertex_id)

    sim_between <- function() {
      # Simulate one between-block network. The output format is a sparse matrix.
      between_conn <- simulate_between_network(
        numOfVertices = N_node,
        list_feature_adjmat = list_feature_matrices,
        coef_between = coef_between_block,
        block_membership = sorted_dataframe$block,
        directed = FALSE
      )
      # Convert the sparse matrix into an edgelist.
      between_conn <-
        as.data.frame(Matrix::summary(between_conn)) %>%
        dplyr::select(i, j) %>%
        dplyr::rename(X1 = i, X2 = j)

      between_conn <-
        as.matrix(between_conn)
      attr(between_conn, "n") <- N_node
      attr(between_conn, "vnames") <- sorted_dataframe$vertex_id
      attr(between_conn, "directed") <- FALSE
      attr(between_conn, "bipartite") <- FALSE
      attr(between_conn, "loops") <- FALSE
      attr(between_conn, "class") <- c("edgelist", "matrix")
      # Return the output
      return(between_conn)
    }

    # If you simulate just one between-block network:
    if (n_sim == 1) {
      output <- sim_between()
      if (verbose > 0) {
        message("Simulated one between-block network.")
      }
      # Return the output
      return(output)
    }

    # If you simulate more than one between-block network:
    if (n_sim > 1) {
      output <-
        foreach(i = 1:n_sim) %do% {
          net <- sim_between()
          if (verbose > 0) {
            message(glue::glue("Finished between-block network simulation {i} of {n_sim}."))
          }
          return(net)
        }
      if (verbose > 0) {
        message("Finished simulating between-block networks.")
      }
      return(output)
    }
  }

  # If use_fast_between_simulation = FALSE, simulate between-block networks using ergm's simulate().
  if (!use_fast_between_simulation) {
    # Generate a between-block seed network.
    seed_network_between <-
      generate_seed_network(formula_for_simulation,
        sorted_dataframe,
        edgelist = seed_edgelist_between,
        directed = FALSE
      )
    # Extract the RHS of the formula
    formula_rhs <- as.character(formula_for_simulation)[3]
    # If prevent_duplicate, set the probability that within-block links are generated to be almost zero.
    if (prevent_duplicate) {
      sim_formula <- as.formula(glue::glue("seed_network_between ~ {formula_rhs} + nodematch(\"block\")"))
      coef_between_block <- c(coef_between_block, -1000)
    } else {
      sim_formula <- as.formula(glue::glue("seed_network_between ~ {formula_rhs}"))
    }

    # Simulate within-block links
    between_conn <- ergm::simulate_formula(
      object = sim_formula,
      nsim = n_sim,
      coef = coef_between_block,
      seed = seed,
      control = ergm_control,
      output = output,
      verbose = verbose
    )
    # Return the output
    return(between_conn)
  }
}


#' Combine within- and between-block edges while removing duplicated links.
#' @param edgelist_within an within-block edgelist
#' @param edgelist_between a between-block edgelist (Potentially, there are edges that also appear in the within-block edgelist)
#' @param use_fast_between_simulation If `TRUE`, this function uses an effcient way to simulate a between-block network. If the network is very large, you should consider using this option.
combine_within_between_edges <- function(edgelist_within, edgelist_between, use_fast_between_simulation) {
  # Store network information necessary to reconstruct a network object later.
  n <- attr(edgelist_within, "n")
  vnames <- attr(edgelist_within, "vnames")
  directed <- attr(edgelist_within, "directed")
  bipartite <- attr(edgelist_within, "directed")
  loops <- attr(edgelist_within, "loops")

  # Convert the edgelists into data frames
  edgelist_within <- data.frame(edgelist_within)
  colnames(edgelist_within) <- c("tail", "head")
  edgelist_between <- data.frame(edgelist_between)
  colnames(edgelist_between) <- c("tail", "head")

  # From df_between, remove between-block edges that also appear in df_within.
  # To do so, we use dplyr::anti_join.
  # If use_fast_between_simulation, skip this step.
  if (!use_fast_between_simulation) {
    edgelist_between <-
      dplyr::anti_join(edgelist_between, edgelist_within, by = c("tail", "head"))
  }

  # Bind within- and between-block edges
  output <-
    rbind(edgelist_within, edgelist_between) %>%
    dplyr::arrange(tail) %>%
    as.matrix()

  # Attach network info
  attr(output, "n") <- n
  attr(output, "vnames") <- vnames
  attr(output, "directed") <- directed
  attr(output, "bipartite") <- bipartite
  attr(output, "loops") <- loops
  class(output) <- c("matrix_edgelist", "edgelist", class(output))

  # Return the output
  return(output)
}


#' Create a final network
#' @param formula_for_simulation_within formula for simulating a within network
#' @param formula_for_simulation_between formula for simulating a between network
#' @param sorted_dataframe a data frame generated by `sort_block_membership`
#' @param edgelist an edgelist that contain both within- and between-block edges without duplication
#' @importFrom foreach foreach %do%
generate_network_for_output <- function(formula_for_simulation_within, formula_for_simulation_between, sorted_dataframe, edgelist) {
  directed <- attr(edgelist, "directed")
  # Initialize a network object from the edgelist
  g <- network::network(edgelist, directed = directed, matrix.type = "edgelist")

  # Get variable names from formula (extract strings sandwiched by double quotes)
  ## From the within-block network formula
  list_varname_within <- as.character(formula_for_simulation_within)[3]
  list_varname_within <- unlist(stringr::str_extract_all(string = list_varname_within, pattern = '"[^"]*"'))
  list_varname_within <- stringr::str_remove_all(string = list_varname_within, pattern = '\"')
  ## From the between-block network formula
  list_varname_between <- as.character(formula_for_simulation_between)[3]
  list_varname_between <- unlist(stringr::str_extract_all(string = list_varname_between, pattern = '"[^"]*"'))
  list_varname_between <- stringr::str_remove_all(string = list_varname_between, pattern = '\"')
  # Combine them.
  list_varname <- unique(c(list_varname_within, list_varname_between))

  # Attach vertex features
  foreach(i = 1:length(list_varname)) %do% {
    feature <-
      dplyr::pull(sorted_dataframe, list_varname[i])
    network::set.vertex.attribute(g, attrname = list_varname[i], value = feature)
  }
  # Attach block info
  network::set.vertex.attribute(g, attrname = "block", value = sorted_dataframe$block)
  # Return the network
  return(g)
}



#' Converts an edgelist into a matrix of sufficient network statistics
#' @param net the net to extract the covariates from
#' @param edgelist the edgelist
#' @param between_formula the formula for the between connections
#' @return a matrix of sufficient network statistics
edgelist_to_stats <- function(net, edgelist, between_formula) {
  g_copy <- network::network(edgelist, matrix.type = "edgelist", directed = FALSE, loops = FALSE)
  attrs <- network::list.vertex.attributes(net)
  for (attr in attrs) {
    network::set.vertex.attribute(g_copy, attr, network::get.vertex.attribute(net, attr))
  }
  form <- as.formula(paste("g_copy", "~", as.character(between_formula)[3]))
  summary(form)
}

#' Converts a list of edgelists into a data frame of network statistics
#' @param edgelists the list of edgelists
#' @param between_formula the formula for the between connections
#' @return a data frame of sufficient network statistics
get_between_stats <- function(edgelists, between_formula) {
  net <- get(as.character(between_formula)[2], envir = environment(between_formula))
  between_stats <-
    edgelists %>%
    purrr::map(function(el) {
      edgelist_to_stats(net, el, between_formula)
    }) %>%
    purrr::reduce(function(a, b) {
      rbind(a, b)
    })
  rownames(between_stats) <- NULL
  data.frame(between_stats)
}


#' Obtains network statistics based on MCMC simulations including only the
#' within-blocks connections.
#' @param formula_for_simulations formula for simulating a network
#' @param data_for_simulation a data frame that contains vertex id, block membership, and vertex features.
#' @param colname_vertex_id a column name in the data frame for the vertex ids
#' @param colname_block_membership a column name in the data frame for the block membership
#' @param seed_edgelist an edgelist used for creating a seed network. It should have the "edgelist" class
#' @param coef_within_block a vector of within-block parameters. The order of the parameters should match that of the formula.
#' @param output The desired output of the simulation (any of `stats`, `network` or `edgelist`). Defaults to `stats`
#' @param ergm_control auxiliary function as user interface for fine-tuning ERGM simulation
#' @param seed seed value (integer) for network simulation.
#' @param n_sim number of networks generated
#' @param verbose If this is TRUE/1, the program will print out additional information about the progress of simulation.
#' @param ... arguments to be passed to low level functions
#' @importFrom magrittr %>% %<>%
#' @export
simulate_hergm_within <- function(formula_for_simulation,
                                  data_for_simulation,
                                  colname_vertex_id,
                                  colname_block_membership,
                                  coef_within_block,
                                  seed_edgelist = NULL,
                                  output = 'stats',
                                  ergm_control = ergm::control.simulate.formula(),
                                  seed = NULL,
                                  n_sim = 1,
                                  verbose = 0,
                                  ...) {

  if (verbose > 0) {
    message("Sorting the given dataset by the block ids.")
  }
  # Create a data frame from block memberships, vertex ids, and vertex covariates, sorted by block ids.
  sorted_dataframe <- sort_block_membership(data_for_simulation, colname_vertex_id, colname_block_membership)
  rm(data_for_simulation)

  # Create seed networks from which networks will be simulated.
  if (verbose > 0) {
    message("Creating a seed within-block network.")
  }

  # Seed edgelist
  if (!is.null(seed_edgelist)) {
    # Validate the `seed_edgelist` object.
    if (!any(class(seed_edgelist) == "edgelist")) {
      stop("\nThe given edgelist must have the 'edgelist' class, i.e., it must have an edgelist, the number of vertices, and vertex names.")
    }
    # Arrange the given edgelist
    seed_edgelist_within <- arrange_edgelist(seed_edgelist, sorted_dataframe)$edgelist_within
  } else {
    seed_edgelist_within <- NULL
  }

  seed_network_within <-
    generate_seed_network(formula_for_simulation, sorted_dataframe, directed = FALSE)

  # Create formula for simulating within-block networks.
  ## Extract the RHS of the given formula
  formula_rhs <- as.character(formula_for_simulation)[3]
  ## Create a formula for simulation
  ### For within-block network
  formula_for_simulation_within <- as.formula(glue::glue("seed_network_within ~ {formula_rhs}"))

  # Send a message about which value of a coefficient is attached to which term.
  # This is to make sure that `coef_within_block` are correctly specified as the user intends.
  names(coef_within_block) <- statnet.common::list_rhs.formula(formula_for_simulation_within) %>%
    as.character %>%
    stringr::str_replace_all("\"", "'")

  # Simulate connections
  if (verbose > 0) {
    message("Simulating within-block networks.")
  }
  simulation_within <- draw_within_block_connection(
    seed_network = seed_network_within,
    formula_for_simulation = formula_for_simulation_within,
    coef_within_block = coef_within_block,
    ergm_control = ergm_control,
    output = output,
    seed = seed,
    n_sim = n_sim,
    verbose = verbose
  )

  if(output == 'stats'){
    simulation_within %<>% as.data.frame
  }

  return(simulation_within)
}
