test_that("computing pi_d1x0 works", {
  rm(list = ls())
  # Number of nodes
  N <- 12
  # Number of clusters
  K <- 3

  # Create an adjacency matrix
  edgelist <-
    tibble::tibble(
      tail = 1:N,
      head = 1:N
    ) %>%
    tidyr::expand(tail, head) %>%
    dplyr::filter(tail < head) %>%
    dplyr::mutate(connect = rep(0:1, nrow(.) / 2)) %>%
    dplyr::filter(connect == 1)

  net <- network::network(edgelist, matrix.type = "edgelist", directed = FALSE)
  adj <- network::as.matrix.network.adjacency(net)
  adj <- as(adj, "dgCMatrix")

  # Create feature matrices
  x <- as.integer(unlist(purrr::rbernoulli(n = N)))
  S <- Matrix::sparseMatrix(i = {}, j = {}, dims = c(N, N))
  S <- as(S, "dgCMatrix")
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        s_ij <- ifelse(x[i] == x[j], 1, 0)
        S[i, j] <- s_ij
      }
    }
  }

  y <- as.integer(unlist(purrr::rbernoulli(n = N)))
  V <- Matrix::sparseMatrix(i = {}, j = {}, dims = c(N, N))
  V <- as(V, "dgCMatrix")
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        v_ij <- ifelse(y[i] == y[j], 1, 0)
        V[i, j] <- v_ij
      }
    }
  }

  z <- as.integer(unlist(purrr::rbernoulli(n = N)))
  W <- Matrix::sparseMatrix(i = {}, j = {}, dims = c(N, N))
  W <- as(W, "dgCMatrix")
  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        w_ij <- ifelse(z[i] == z[j], 1, 0)
        W[i, j] <- w_ij
      }
    }
  }



  # Create a N x K matrix whose (i, k) element represents the probability that node i belongs to block k.
  tau <-
    matrix(c(
      0.2, 0.5, 0.3,
      0.4, 0.4, 0.2,
      0.1, 0.4, 0.5,
      0.4, 0.4, 0.2,
      0.1, 0.1, 0.8,
      0.05, 0.05, 0.9,
      0.8, 0.1, 0.1,
      0.3, 0.4, 0.3,
      0.1, 0.8, 0.1,
      0.5, 0.4, 0.1,
      0.3, 0.3, 0.4,
      0.8, 0.1, 0.1
    ),
    nrow = K, ncol = N
    )
  tau <- t(tau)

  # Compute the true denominator
  one <- matrix(1, nrow = N, ncol = N)
  mat <- (one - S) * (one - V) * (one - W)
  diag(mat) <- 0
  denom_for_pi0_true <- t(tau) %*% mat %*% tau
  denom_for_pi0_true <- as.matrix(denom_for_pi0_true)

  # Compute the true denominator in a naive way
  denom_for_pi0_naive <- matrix(0, nrow = K, ncol = K)
  for (k in 1:K) {
    for (l in 1:K) {
      for (i in 1:N) {
        for (j in 1:N) {
          if (i != j & S[i, j] == 0 & V[i, j] == 0 & W[i, j] == 0) {
            denom_for_pi0_naive[k, l] <- denom_for_pi0_naive[k, l] + tau[i, k] * tau[j, l]
          }
        }
      }
    }
  }

  # Check if they are the same. This verifies that the formula is correct.
  expect_equal(denom_for_pi0_true, denom_for_pi0_naive, check.attributes = FALSE, tolerance = 1e-10)

  # Compute the denominator using the c++ function
  denom <- get_matrix_for_denominator(N, list(S, V, W))
  denom_for_pi0 <- compute_denominator_for_pi_d1x0(N, K, denom, tau, verbose = 0)

  # Check if the computed matrix is correct.
  expect_equal(denom_for_pi0, denom_for_pi0_naive, check.attributes = FALSE, tolerance = 1e-10)


  # Compute true pi1 in a naive way
  pi1 <- matrix(0, nrow = K, ncol = K)
  for (k in 1:K) {
    for (l in 1:K) {
      for (i in 1:N) {
        for (j in 1:N) {
          if (i != j & adj[i, j] == 1 & S[i, j] == 0 & V[i, j] == 0 & W[i, j] == 0) {
            pi1[k, l] <- pi1[k, l] + tau[i, k] * tau[j, l]
          }
        }
      }
    }
  }
  pi1_true <- pi1 / denom_for_pi0_naive

  # Remove extremely small values in pi1
  minPi <- 1e-4
  for (k in 1:K) {
    for (l in 1:K) {
      if (pi1_true[k, l] < minPi) {
        pi1_true[k, l] <- minPi
      }
    }
  }

  # Compute pi0 using the Rcpp function
  list_multiplied_adjmat <- get_elementwise_multiplied_matrices(adj, list(S, V, W))
  list_multiplied_adjmat[[1]] <- denom
  pi1 <- compute_pi_d1x0(N, K, list_multiplied_adjmat, tau, verbose = 0)

  # Check if the computed conditional probability is correct.
  expect_equal(pi1, pi1_true, check.attributes = FALSE, tolerance = 1e-10)
})
