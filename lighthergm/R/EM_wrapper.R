# Wrapper function that conducts clustering.
EM_wrapper <-
  function(network,
           formula,
           n_clusters,
           n_em_step_max,
           min_size = 2,
           initialization_method = 1,
           use_infomap_python = FALSE,
           seed_infomap = NULL,
           initialized_cluster_data = NULL,
           clustering_with_features = FALSE,
           list_multiplied_feature_matrices = NULL,
           verbose = 0,
           weight_for_initialization = 1000,
           compute_pi = FALSE,
           check_alpha_update = FALSE,
           check_block_membership = FALSE,
           EM_restart_object = NULL,
           minAlpha = 1e-6,
           ...) {
    n_nodes <- as.integer(network$gal$n)
    is_undirected <- !network$gal$directed

    eigenvectors_sparse_cached <- memoise::memoise(eigenvectors_sparse)

    # If restart_object is NULL, initialize memberships.
    if (is.null(EM_restart_object)) {
      # Step 0: Initialization
      n_nodes <- as.integer(network$gal$n)
      is_undirected <- !network$gal$directed

      if (verbose > 0) message("Converting network to edgelist...")
      edgelist <- network::as.edgelist(network)

      if (verbose > 0) message("Converting edgelist to sparse matrix...")
      g <- Matrix::sparseMatrix(i = edgelist[, 1], j = edgelist[, 2], x = 1, dims = c(n_nodes, n_nodes), symmetric = is_undirected)

      ## Calculate network statistics needed for variational approximation of p(Z|X)
      if (verbose > 0) message("\nStep 1: Initialize z")

      # Step 1: Get initial estimate of Z memberships
      # If you already have an initialized clustering result, start the EM iterations from the given clustering result.
      if (!is.null(initialized_cluster_data)) {
        z_memb_init <- factor(initialized_cluster_data, levels = 1:n_clusters)
        z_memb <- check_clusters_sparse(z_memb_init, g, n_clusters, eigenvectors_sparse_cached, min_size, verbose = verbose)
      }

      # Start cluster initialization
      else if (initialization_method == 1) # Infomap
        {
          if (verbose > 0) message("Using Infomap to initialize clustering step...")

          if (use_infomap_python) {
            # Use python's infomap, which is much faster than R's infomap.
            # reticulate::py_config()
            temp_dir <- tempdir()
            temp <- tempfile(fileext = c(".csv"))
            temp_clu <- stringr::str_replace(string = temp[1], pattern = ".csv", replacement = ".clu")
            readr::write_delim(as.data.frame(edgelist), file = temp[1], col_names = FALSE)
            # Implement infomap
            # If the seed for infomap is specified, implement infomap using the seed.
            if (is.numeric(seed_infomap)) {
              system(glue::glue("infomap --flow-model undirected --preferred-number-of-modules {n_clusters} --clu {temp[1]} {temp_dir[1]} --seed {seed_infomap}"))
            } else {
              system(glue::glue("infomap --flow-model undirected --preferred-number-of-modules {n_clusters} --clu {temp[1]} {temp_dir[1]}"))
            }
            # Read the initialized cluster data
            b <- readr::read_delim(temp_clu[1], delim = " ", skip = 9, col_names = c("node_id", "block", "flow"), col_types = "iid")
            b <- as.numeric(dplyr::arrange(b, by_group = node_id)$block)
            z_memb_init <- factor(b, levels = 1:n_clusters)
            z_memb <- check_clusters_sparse(z_memb_init, g, n_clusters, eigenvectors_sparse_cached, min_size, verbose = verbose)
            # Remove tempfiles just in case
            invisible(file.remove(temp))
            invisible(file.remove(temp_clu))
          }
          else {
            # Use igraph's infomap
            set.seed(seed_infomap)
            b <- igraph::cluster_infomap(intergraph::asIgraph(network))
            z_memb_init <- factor(b$membership, levels = 1:n_clusters)
            z_memb <- check_clusters_sparse(z_memb_init, g, n_clusters, eigenvectors_sparse_cached, min_size, verbose = verbose)
          }
        }
      else if (initialization_method == 2) {
        if (verbose > 0) message("Initializing with uniformly distributed clusters...")
        z_memb <- z_memb_init <- factor(sample.int(n_clusters, size = n_nodes, replace = TRUE))
      }
      else # Spectral clustering
      {
        z_memb <- spec_clust_sparse(g, n_clusters, eigenvectors_sparse_cached)
        z_memb_init <- z_memb
      }

      # Step 2: Find A(Z=z) ~ P(Z=z|X=x)
      if (verbose > 0) {
        message(paste("\n\nStep 2: Find variational approximation A(Z=z) ~ P(Z=z|X=x)", sep = ""))
      }
      block_sizes <- table(z_memb)
      memb_inds <- as.numeric(names(block_sizes))
      n_blocks <- length(block_sizes)

      # Initialize alpha
      if (verbose > 0) message("Initializing alpha...")
      alpha <- matrix(1, n_nodes, n_clusters)
      if (verbose > 0) message("Replacing zeros in alpha...")
      for (i in which(!is.na(z_memb))) {
        alpha[i, z_memb[i]] <- weight_for_initialization
      }
      if (verbose > 0) message("Normalizing alpha...")
      alpha <- alpha / Matrix::rowSums(alpha)

      # Keep track of alpha, z, changes in alpha, and the lower bound
      list_alpha <- list()
      list_z <- list()
      change_in_alpha <- c()
      lower_bound <- c()
      if (check_alpha_update) {
        list_alpha[[1]] <- alpha
      }
      if (check_block_membership) {
        list_z[[1]] <- z_memb
      }
      # Set the counter of EM iteration
      counter_e_step <- 0
    }
    # If restart_object is given, restore the values.
    else {
      if (verbose > 0) message("Making alpha matrix dense...")
      alpha <- as.matrix(EM_restart_object$alpha)
      alpha[alpha == 0] <- minAlpha
      list_alpha <- EM_restart_object$list_alpha
      z_memb_init <- EM_restart_object$z_memb_init
      list_z <- EM_restart_object$list_z
      change_in_alpha <- EM_restart_object$change_in_alpha
      lower_bound <- EM_restart_object$lower_bound
      counter_e_step <- EM_restart_object$counter_e_step
      g <- EM_restart_object$adjacency_matrix
      n_blocks <- n_clusters
    }

    # Get a sparse feature adjacency matrix if clustering with features
    if (clustering_with_features) {
      # If a list of feature matrices is given from outside, skip this step.
      if (!is.null(list_multiplied_feature_matrices)) {
        if (verbose > 0) {
          message("Skipped forming features matrix")
        }
      } else {
        if (verbose > 0) {
          message("Forming features matrix")
        }
        list_feature_matrices <- get_list_sparse_feature_adjmat(network, formula)
        if (verbose > 0) {
          message("Got features matrix")
        }
      }
    }

    # Get multiplied feature matrices if clustering with features
    if (clustering_with_features) {
      if (is.null(list_multiplied_feature_matrices)) {
        if (verbose > 0) {
          message("Getting multiplied matrices used for E step")
        }
        list_multiplied_feature_matrices <- get_elementwise_multiplied_matrices(g, list_feature_matrices)
        if (verbose > 0) {
          message("Got multiplied matrices")
        }
        rm(list_feature_matrices)
      } else {
        message("Skipped getting multiplied matrices")
      }
    }

    # Compute gamma (parameter for multinomial distribution)
    gamma <- colMeans(alpha)


    # If clustering with features:
    if (verbose > 0) {
      message(paste("EM algorithm started at ", counter_e_step))
      EM_start <- Sys.time()
      message(paste("Started at ", EM_start))
    }

    # This is for internal use
    inner_counter <- 0

    if (clustering_with_features) {
      # EM iterations start
      repeat {
        inner_counter <- inner_counter + 1
        counter_e_step <- counter_e_step + 1

        if (verbose > 0) {
          message(paste("EM algorithm iteration: ", counter_e_step))
          iter_start <- Sys.time()
          message(glue::glue("EM algorithm iteration {counter_e_step} started at {iter_start}"))
        }

        # If you would like to see the change in \alpha during EM iterations, keep \alpha_{t-1} before updating.
        if (check_alpha_update) {
          alpha_prev <- rlang::duplicate(alpha)
        }

        # Update alpha
        alpha_LB <-
          run_EM_with_features(
            as.integer(n_nodes),
            as.integer(n_clusters),
            as.double(gamma),
            list_multiplied_feature_matrices,
            alpha,
            verbose
          )

        alpha <- alpha_LB[[1]]
        gamma <- colMeans(alpha)

        # Keep track of the lower bound:
        lower_bound <- c(lower_bound, alpha_LB[[2]])

        # If you would like to see the change in \alpha during EM iterations, print max_{i, k} |\alpha_{t}(i, k) - \alpha_{t-1}(i, k)|.
        if (check_alpha_update) {
          list_alpha[[counter_e_step + 1]] <- alpha
          change_in_alpha <- c(change_in_alpha, max(abs(alpha - alpha_prev)))
        }

        # If you would like to keep track on block memberships over EM iterations:
        if (check_block_membership) {
          z <- factor(apply(alpha, 1, which.max), levels = 1:n_clusters)
          list_z[[counter_e_step + 1]] <- z
        }

        if (verbose > 0) {
          iter_end <- Sys.time()
          message(glue::glue("EM algorithm iteration {counter_e_step} finished at {iter_end}"))
          message(glue::glue("EM iterations {counter_e_step} took {signif(difftime(iter_end, iter_start, units = 'mins'), digits = 3)} minutes."))
        }

        if (inner_counter %% 5) {
          gc()
        }

        if (inner_counter >= n_em_step_max) {
          break
        }
      }
      # End of EM iteration
      if (verbose > 0) {
        EM_end <- Sys.time()
        message("Finished the EM cycle at ", EM_end)
        diff <- EM_end - EM_start
        message(glue::glue("{inner_counter} EM iterations took {signif(difftime(EM_end, EM_start, units = 'mins'), digits = 3)} minutes.\n"))
        message(glue::glue("Total EM iterations: {counter_e_step}"))
      }

      # Compute pi if necessary
      Pi <- NULL
      if (compute_pi) {
        Pi <- compute_pi_with_features(n_nodes, n_clusters, list_multiplied_feature_matrices, alpha)
      }
    } else { # If clustering without features:
      list_multiplied_feature_matrices <- NULL
      # EM iterations start
      repeat {
        counter_e_step <- counter_e_step + 1
        inner_counter <- inner_counter + 1

        if (verbose > 0) {
          message(paste("EM algorithm iteration: ", counter_e_step))
          iter_start <- Sys.time()
          message(glue::glue("EM algorithm iteration {counter_e_step} started at {iter_start}"))
        }

        # If you would like to see the change in \alpha during EM iterations, keep \alpha_{t-1} before updating.
        if (check_alpha_update) {
          alpha_prev <- rlang::duplicate(alpha)
        }

        alpha_LB <-
          run_EM_without_features(
            n_nodes,
            n_blocks,
            gamma,
            alpha,
            g,
            verbose
          )

        # \alpha and \gamma
        alpha <- alpha_LB[[1]]
        gamma <- colMeans(alpha)

        # Keep track of the lower bound:
        lower_bound <- c(lower_bound, alpha_LB[[2]])

        # If you would like to see the change in \alpha during EM iterations, print max_{i, k} |\alpha_{t}(i, k) - \alpha_{t-1}(i, k)|.
        if (check_alpha_update) {
          list_alpha[[counter_e_step + 1]] <- alpha
          change_in_alpha <- c(change_in_alpha, max(abs(alpha - alpha_prev)))
        }

        # If you would like to keep track on block memberships over EM iterations:
        if (check_block_membership) {
          z <- factor(apply(alpha, 1, which.max), levels = 1:n_clusters)
          list_z[[counter_e_step + 1]] <- z
        }

        if (verbose > 0) {
          iter_end <- Sys.time()
          message(glue::glue("EM algorithm iteration {counter_e_step} finished at {iter_end}"))
          message(glue::glue("EM iterations {counter_e_step} took {signif(difftime(iter_end, iter_start, units = 'mins'), digits = 3)} minutes."))
        }

        if (inner_counter %% 5) {
          gc()
        }

        if (inner_counter >= n_em_step_max) {
          break
        }
      }

      # End of EM iteration
      if (verbose > 0) {
        EM_end <- Sys.time()
        message("Finished the EM cycle at ", EM_end)
        diff <- EM_end - EM_start
        message(glue::glue("{inner_counter} EM iterations took {signif(difftime(EM_end, EM_start, units = 'mins'), digits = 3)} minutes.\n"))
        message(glue::glue("Total EM iterations: {counter_e_step}"))
      }

      # Compute pi if necessary
      Pi <- NULL
      if (compute_pi) {
        Pi <- (t(alpha) %*% g %*% alpha) / compute_sumTaus(n_nodes, n_blocks, alpha)
      }
    }

    # Estimate the block membership
    z_memb <-
      factor(apply(alpha, 1, which.max), levels = 1:n_clusters)
    # Keep the final clustering result before check_cluster_sparse()
    z_memb_final_before_kmeans <- z_memb

    if (verbose > 0) message("Making alpha matrix sparse...")
    alpha[alpha <= minAlpha] <- 0
    alpha <- as(alpha, "dgCMatrix")

    # Check bad clusters
    z_memb_final <-
      factor(check_clusters_sparse(z_memb, g, n_clusters, eigenvectors_sparse_cached, min_size, verbose = verbose),
        levels = 1:n_clusters
      )

    # Return the output
    return(
      list(
        Pi = Pi,
        z_memb_init = z_memb_init,
        list_z = list_z,
        z_memb_final_before_kmeans = z_memb_final_before_kmeans,
        z_memb_final = z_memb_final,
        list_alpha = list_alpha,
        change_in_alpha = change_in_alpha,
        lower_bound = lower_bound,
        counter_e_step = counter_e_step,
        adjacency_matrix = g,
        alpha = alpha,
        list_multiplied_feature_matrices = list_multiplied_feature_matrices
      )
    )
  }

# ------------------------------------------------------------------------
# -------------- Auxiliary functions -------------------------------------
# ------------------------------------------------------------------------

#' a function that gets a network object from a formula
#' @param form a formula that contains a network object in the left-hand side
#' @param n_clusters the maximum number of blocks
#' @param loopswarning if TRUE, this function warns that the network contains loops.
hergm.getnetwork <- function(form, n_clusters, loopswarning = TRUE) {
  if ((dc <- data.class(form)) != "formula") {
    stop(paste("Invalid formula of class ", dc))
  }
  trms <- terms(form)
  if (trms[[1]] != "~") {
    stop("Formula must be of the form 'network ~ model'.")
  }

  nw.env <- environment(form)
  if (!exists(x = paste(trms[[2]]), envir = nw.env)) {
    stop(paste("The network in the formula '", capture.output(print(form)), "' cannot be found.", sep = ""))
  }
  nw <- try(
    {
      tmp <- eval(trms[[2]], envir = nw.env)
      if (network::is.network(tmp)) {
        return(tmp)
      } else {
        return(network::as.network(tmp))
      }
    },
    silent = TRUE
  )
  if (inherits(nw, "try-error")) {
    stop("Invalid network. Is the left-hand-side of the formula correct?")
  }
  if (loopswarning) {
    e <- as.edgelist(nw)
    if (any(e[, 1] == e[, 2])) {
      warning("This network contains loops")
    } else if (has.loops(nw)) {
      warning("This network is allowed to contain loops")
    }
  }
  nw$terms <- 1
  if (is.null(n_clusters)) {
    nw$n_clusters <- nw$gal$n
  } else {
    nw$n_clusters <- n_clusters
  }
  # print("nw$n_clusters")
  # print(nw$n_clusters)
  model <- ergm::ergm_model(form, nw, drop = FALSE, expanded = TRUE)
  Clist <- ergm::ergm.Cprepare(nw, model)
  nw$terms <- Clist$nterms
  return(nw)
}

permute_tau <- function(tau, labels) {
  new_tau <- tau
  for (i in 1:length(labels)) {
    # temp[z_memb == i] <- labels[as.numeric(names(labels))==i]
    new_tau[, i] <- tau[, labels[i]]
  }
  new_tau
}


#' function for spectral clustering
#' @param network a sparse adjacency matrix
#' @param n_clusters number of specified clusters
spec_clust_sparse <- function(network, n_clusters, eigenvectors_fn) {
  n <- nrow(network)
  n_vec <- ceiling(sqrt(n))
  message("Calculating eigenvectors...")
  b <- eigenvectors_fn(network, n_vec)
  message("K-means clustering...")
  c <- kmeans(
    b,
    centers = n_clusters,
    nstart = 100,
    iter.max = 20,
    algorithm = "Hartigan-Wong"
  )
  c.ind <- c$cluster
  z_memb <- factor(c.ind, levels = 1:n_clusters)
  z_memb
}


check_clusters_sparse <-
  function(z_memb, network, n_clusters, eigenvectors_fn, min_size = 2, verbose = verbose) {
    n <- nrow(network)
    n_vec <- ceiling(sqrt(n))

    if (verbose > 0) {
      message("Eigenvalue decomposition")
    }

    b <- eigenvectors_fn(network, n_vec)
    z_memb <- factor(z_memb, levels = 1:n_clusters)

    if (verbose > 0) {
      message("Checking for bad clusters...")
    }

    counter_bad_cluster <- 0
    length_prev_bad_clusters <- NULL
    repeat {
      z_memb_temp <- as.vector(z_memb)
      z_memb <- factor(z_memb_temp, levels = 1:n_clusters)
      bad_clusters <- which(table(z_memb) < min_size)

      if (verbose > 0) {
        message(paste(c("Remaining bad clusters:", length(bad_clusters))))
      }

      # If the number of bad clusters doesn't change over one iteration:
      if (!is.null(length_prev_bad_clusters)) {
        if (length(bad_clusters) == length_prev_bad_clusters) {
          counter_bad_cluster <- counter_bad_cluster + 1
        } else {
          counter_bad_cluster <- 0
        }
      }
      length_prev_bad_clusters <- length(bad_clusters)

      # If removing bad clusters seems to fail, stop the process.
      if (counter_bad_cluster >= 100) {
        stop("\nRemoving bad clusters failed. The number of pre-specified blocks might be too high.")
      }

      if (length(bad_clusters) == 0) {
        break
      }

      bad_cluster <- which(table(z_memb) < 2)[1]
      donor <- which.max(table(z_memb))
      indeces <- c(bad_clusters, donor)
      temp <-
        kmeans(
          b[z_memb %in% indeces, ],
          centers = 2,
          nstart = 100,
          iter.max = 20,
          algorithm = "Hartigan-Wong"
        )
      for (i in 1:length(indeces)) {
        z_memb_temp[z_memb %in% indeces][temp$cluster == i] <- indeces[i]
      }
      z_memb <- factor(z_memb_temp, levels = 1:n_clusters)
    }

    if (verbose > 0) {
      message("Done checking clusters")
    }

    return(z_memb)
  }
