#' Hierarchical exponential-family random graph models (HERGMs) with local dependence
#' @description
#' The function hergm estimates and simulates three classes of hierarchical exponential-family random graph models.
#' @useDynLib lighthergm
#' @importFrom Rcpp sourceCpp
#' @param object A formula or `lighthergm` class object. A `lighthergm` is returned by `hergm()`.
#' When you pass a `lighthergm` class object to `hergm()`, you can restart the EM step.
#' @param n_clusters The number of blocks. This must be specified by the user.
#' When you pass a "lighthergm" class object to `hergm()`, you don't have to specify this argument.
#' @param n_cores The number of CPU cores to use.
#' @param block_membership The pre-specified block memberships for each node.
#' If `NULL`, the latent community structure is estimated, assuming that the number of communities is `n_clusters`.
#' @param estimate_parameters If `TRUE`, both clustering and parameter estimation are implemented.
#' If `FALSE`, only clustering is executed.
#' @param verbose A logical or an integer: if this is TRUE/1,
#' the program will print out additional information about the progress of estimation and simulation.
#' A higher value yields lower level information.
#' @param n_em_step_max The maximum number of EM iterations.
#' Currently, no early stopping criteria is introduced. Thus `n_em_step_max` EM iterations are exactly implemented.
#' @param initialization_method Cluster initialization method.
#' If `1` (the default), `igraph`'s infomap is implemented.
#' If `2`, the initial clusters are randomly uniformally selected.
#' If `3`, spectral clustering is conducted.
#' @param use_infomap_python If `TRUE`, the cluster initialization is implemented using Pythons' infomap.
#' When using this, make sure that Python's infomap is installed and callable from `system().
#' @param seed_infomap seed value (integer) for Python's infomap.
#' @param weight_for_initialization weight value used for cluster initialization. The higher this value, the more weight is put on the initialized alpha.
#' @param seeds seed value (integer) for the random number generator
#' @param initialized_cluster_data initialized cluster data from which the EM iterations begin.
#' This can be either a vector of block affiliations of each node or initialized cluster data by Python's infomap (given by .clu format).
#' @param method_second_step If "MPLE" (the default), then the maximum pseudolikelihood estimator is implemented when estimating the within-block network model.
#' If "MLE", then an approximate maximum likelihood estimator is conducted.
#' @param clustering_with_features If `TRUE`, clustering is implemented using the discrete covariates specified in the formula.
#' @param list_multiplied_feature_matrices a list of multiplied feature adjacency matarices necessary for EM step.
#' If `NULL`, `hergm()` automatically calculates. Or you can calculate by `compute_multiplied_feature_matrices()`.
#' @param fix_covariate_parameter If `TRUE`, when estimating the within-block network model,
#' parameters for covariates are fixed at the estimated of the between-block network model.
#' @param compute_pi If `TRUE`, this function keeps track of pi matrices at each EM iteration.
#' If the network is large, we strongly recommend to set to be `FALSE`.
#' @param check_alpha_update If `TRUE`, this function keeps track of alpha matrices at each EM iteration.
#' If the network is large, we strongly recommend to set to be `FALSE`.
#' @param check_block_membership If TRUE, this function keeps track of estimated block memberships at each EM iteration.
#' @param cache a `cachem` cache object used to store intermediate calculations such as eigenvector decomposition results.
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
#' @export
hergm <- function(object,
                  n_clusters,
                  n_cores = 1,
                  block_membership = NULL,
                  estimate_parameters = TRUE,
                  verbose = 0,
                  n_em_step_max = 100,
                  initialization_method = 1,
                  use_infomap_python = FALSE,
                  seed_infomap = NULL,
                  weight_for_initialization = 1000,
                  seeds = NULL,
                  initialized_cluster_data = NULL,
                  method_second_step = "MPLE",
                  clustering_with_features = TRUE,
                  list_multiplied_feature_matrices = NULL,
                  fix_covariate_parameter = FALSE,
                  compute_pi = FALSE,
                  check_alpha_update = FALSE,
                  check_block_membership = FALSE,
                  cache = NULL,
                  ...) {
  ###################################################################################
  ###### Preparations for estimation ################################################
  ###################################################################################
  # When the given object is formula:
  if ("formula" %in% class(object)) {
    # If n_cluster is missing and block_membership and initialized_cluster_data are NULL, stop the process.
    if (missing(n_clusters) & is.null(block_membership) & is.null(initialized_cluster_data)) {
      stop("\nThe argument 'n_clusters' is missing. Please specify the number of clusters.")
    }

    # Get the formula
    formula <- object

    # When use_infomap_python = TRUE, check if it is installed.
    if (use_infomap_python) {
      tryCatch(
        {
          system("infomap --version")
        },
        warning = function(e) {
          stop("\nPython's infomap is not installed or your default Python version is 2.x. Please install infomap by
           \nsystem('pip3 install infomap')
           \nIf already installed, make sure that the default Python version is 3.x. You can check the default version by
           \nsystem('python --version')")
        }
      )
    }
    # If the initialized cluster data is given by .clu format:
    if (!is.null(initialized_cluster_data) && all(stringr::str_detect(initialized_cluster_data, ".clu"))) {
      if (verbose > 0) {
        message(paste("Reading initial clusters data from: ", initialized_cluster_data))
      }
      initialized_cluster_data <- readr::read_delim(initialized_cluster_data, delim = " ", skip = 9, col_names = c("node_id", "block", "flow"), col_types = "iid")
      initialized_cluster_data <- as.numeric(dplyr::arrange(initialized_cluster_data, by_group = node_id)$block)
    }

    # Update the number of clusters if `initialized_cluster_data` or `block_membership` is given:
    if (!is.null(initialized_cluster_data)) {
      n_clusters <- length(unique(initialized_cluster_data))
    }
    else if (!is.null(block_membership)) {
      n_clusters <- length(unique(block_membership))
    }

    # Get network object from formula
    network <- hergm.getnetwork(formula, n_clusters)

    # The current hergm doesn't support directed networks.
    if (network$gal$directed) {
      stop("\nThe current hergm doesn't support directed networks. This will be modified in the future.")
    }

    # If list_multiplied_feature_matrices is not NULL, check if the order of features names is the same with that of formula.
    if (!is.null(list_multiplied_feature_matrices)) {
      # If clustering_with_features = FALSE, stop the process.
      if (!clustering_with_features) {
        stop("\nSet clustering_with_features = TRUE since you are going to use vertex features for cluster estimation.")
      }
    }

    # Convert vertex.names into string
    if (!is.character(network::network.vertex.names(network))) {
      network::network.vertex.names(network) <- as.character(network::network.vertex.names(network))
    }

    # An N x K matrix which stores the probability that node i belongs to block k (i = 1,..., N, k = 1,..., K)
    # If all_block_memberships_fixed == TRUE, this matrix remains NULL.
    # sbm_pi <- NULL
  }

  EM_restart_object <- NULL
  # When the given object has a class of "lighthergm":
  if ("lighthergm" %in% class(object)) {
    # Inherit the following arguments from the previous estimation.
    # These arguments are so important that they won't be replaced by the arguments given by the user.
    network <- object$checkpoint$network
    formula <- object$checkpoint$formula
    n_clusters <- object$checkpoint$n_clusters
    clustering_with_features <- object$checkpoint$clustering_with_features
    list_multiplied_feature_matrices <- object$checkpoint$list_multiplied_feature_matrices

    # Inherit the following arguments from the previous estimation if not given by the user.
    # If given by the user, discard the inherited argument and use the given one.
    vec_arguments <-
      c(
        "n_cores",
        "estimate_parameters",
        "verbose",
        "n_em_step_max",
        "seeds",
        "method_second_step",
        "fix_covariate_parameter",
        "compute_pi",
        "check_alpha_update",
        "check_block_membership"
      )

    message("Arguments:")
    message(glue::glue("Number of clusters = {n_clusters}"))
    for (i in 1:length(vec_arguments)) {
      if (do.call(missing, list(vec_arguments[[i]])) == TRUE) {
        assign(vec_arguments[[i]], object$checkpoint[[vec_arguments[[i]]]])
      }
      message(glue::glue("{vec_arguments[[i]]} = {get(vec_arguments[[i]])}"))
    }

    # Prepare an object to restart the EM with.
    EM_restart_object <-
      list(
        alpha = object$alpha,
        list_alpha = object$EM_list_alpha,
        list_z = object$EM_list_z,
        z_memb_init = object$z_memb_final_before_kmeans,
        change_in_alpha = object$EM_change_in_alpha,
        lower_bound = object$EM_lower_bound,
        counter_e_step = object$counter_e_step,
        adjacency_matrix = object$adjacency_matrix
      )
  }

  # If the block_memberships of each node are specified in the variable 'block_membership' as integers between 1 and n_clusters,
  # the specified block memberships are fixed.
  all_block_memberships_fixed <- ifelse(is.null(block_membership), FALSE, TRUE)
  # Make sure that if all_block_memberships_fixed == TRUE, estimate_parameters must also be TRUE.
  estimate_parameters <- ifelse(all_block_memberships_fixed, TRUE, estimate_parameters)

  ###################################################################################
  ###### First step: Estimating block memberships ###################################
  ###################################################################################

  # Estimate the memberships if they are not specified.
  if (!all_block_memberships_fixed) {
    set.seed(seeds[1])
    # Estimate the block memberships
    answer <- EM_wrapper(
      network = network,
      formula = formula,
      n_clusters = n_clusters,
      n_em_step_max = n_em_step_max,
      min_size = 2,
      initialization_method = initialization_method,
      use_infomap_python = use_infomap_python,
      seed_infomap = seed_infomap,
      initialized_cluster_data = initialized_cluster_data,
      clustering_with_features = clustering_with_features,
      list_multiplied_feature_matrices = list_multiplied_feature_matrices,
      verbose = verbose,
      weight_for_initialization = weight_for_initialization,
      compute_pi = compute_pi,
      check_alpha_update = check_alpha_update,
      check_block_membership = check_block_membership,
      EM_restart_object = EM_restart_object,
      cache = cache
    )

    block_membership <- answer$z_memb_final
    initial_block <- answer$z_memb_init
    membership_before_kmeans <- answer$z_memb_final_before_kmeans
    sbm_pi <- answer$Pi
    EM_list_alpha <- answer$list_alpha
    EM_list_z <- answer$list_z
    EM_change_in_alpha <- answer$change_in_alpha
    EM_lower_bound <- answer$lower_bound
    alpha <- answer$alpha
    counter_e_step <- answer$counter_e_step
    adjacency_matrix <- answer$adjacency_matrix
    list_multiplied_feature_matrices <- answer$list_multiplied_feature_matrices
  }
  else {
    if (verbose > 0) {
      message("\nSkipping Steps 1 and 2: z specified")
    }
    param_EM_wrapper <- NULL
    initial_block <- NULL
    membership_before_kmeans <- NULL
    sbm_pi <- NULL
    EM_list_z <- NULL
    EM_list_alpha <- NULL
    EM_change_in_alpha <- NULL
    EM_lower_bound <- NULL
    alpha <- NULL
    counter_e_step <- NULL
    adjacency_matrix <- NULL
    list_multiplied_feature_matrices <- NULL
  }

  ####################################################################################################
  ###### Second step: Estimating between-block parameters ############################################
  ####################################################################################################

  if (estimate_parameters) {
    if (verbose > 0) {
      message("\nEstimate between-block parameters")
    }
    ## Estimate between-block parameters
    est_between <- estimate_between_param(
      formula = formula,
      network = network,
      block = block_membership
    )
  } else {
    # If you don't estimate any parameters...
    est_between <- NULL
  }

  ####################################################################################################
  ###### Third step: Estimating within-block parameters ##############################################
  ####################################################################################################

  # When you estimate with-block parameters, then:
  if (estimate_parameters) {
    if (verbose > 0) {
      message("\n\nStep 3: Estimate parameters conditional on z")
    }
    # If fixing within-block feature parameters:
    if (fix_covariate_parameter) {
      offset_coef <- coef(est_between)
      offset_coef <- offset_coef[-1]
    } else {
      offset_coef <- NULL
    }
    # Estimate within-block parameters
    est_within <-
      estimate_within_params(
        formula,
        network,
        z_memb = block_membership,
        n_cores,
        verbose,
        seeds = NULL,
        method_second_step,
        offset_coef
      )


    ####################################################################################################
    ###### Store the results ###########################################################################
    ####################################################################################################

    estimation_status <- ifelse(est_within$failure, "Estimation failed", "Estimated")
  } else {
    # estimate_parameters = FALSE:
    message("\n")
    estimation_status <- "Not estimated"
    est_within <- NULL
  }

  # Store the given arguments. These will be inherited for the next EM.
  checkpoint <- list(
    network = network,
    formula = formula,
    n_clusters = n_clusters,
    clustering_with_features = clustering_with_features,
    list_multiplied_feature_matrices = list_multiplied_feature_matrices,
    n_cores = n_cores,
    estimate_parameters = estimate_parameters,
    verbose = verbose,
    n_em_step_max = n_em_step_max,
    seeds = seeds,
    method_second_step = method_second_step,
    fix_covariate_parameter = fix_covariate_parameter,
    compute_pi = compute_pi,
    check_alpha_update = check_alpha_update,
    check_block_membership = check_block_membership
  )

  # Store the results in a list
  output <- list(
    partition = block_membership,
    initial_block = initial_block,
    sbm_pi = sbm_pi,
    EM_list_z = EM_list_z,
    EM_list_alpha = EM_list_alpha,
    EM_change_in_alpha = EM_change_in_alpha,
    EM_lower_bound = EM_lower_bound,
    alpha = alpha,
    counter_e_step = counter_e_step,
    adjacency_matrix = adjacency_matrix,
    estimation_status = estimation_status,
    est_within = est_within,
    est_between = est_between,
    checkpoint = checkpoint,
    membership_before_kmeans = membership_before_kmeans
  )

  # Return the output
  return(structure(output, class = "lighthergm"))
}
