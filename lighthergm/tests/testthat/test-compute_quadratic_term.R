test_that("quadratic term calculation without features works", {
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
    dplyr::mutate(connect = as.integer(unlist(purrr::rbernoulli(n = nrow(.), p = 0.5)))) %>%
    dplyr::filter(connect == 1)

  net <- network::network(edgelist, matrix.type = "edgelist", directed = FALSE)
  adj <- network::as.matrix.network.adjacency(net)
  adj <- as(adj, "dgCMatrix")

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

  # Create a K x K matrix whose (k, l) element represents Pr(D_ij = 1 | Z_i = k, Z_j = l).
  sumTaus <- compute_sumTaus(N, K, tau)
  pi <- (t(tau) %*% adj %*% tau) / sumTaus

  # Compute gamma (parameter of multinomial distribution)
  alpha <- colSums(tau)

  # Compute the true quadratic term in a naive way
  A <- matrix(0, nrow = N, ncol = K)

  for (i in 1:N) {
    for (k in 1:K) {
      for (j in 1:N) {
        if (i != j) {
          for (l in 1:K) {
            pi_ij <- pi
            # When D_ij = 0, we must use 1 - pi.
            if (adj[i, j] == 0) {
              pi_ij <- 1 - pi
            }
            a_ij <- tau[j, l] * log(pi_ij[k, l])
            A[i, k] <- A[i, k] + a_ij
          }
        }
      }
    }
  }

  A <- 1 - A / 2

  # Divide A by alpha_{ik}
  A <- A / tau

  A_cpp <- compute_quadratic_term(N, K, alpha, tau, adj, LB = 0)

  # Check if computation works as expected
  expect_equal(A, A_cpp, check.attributes = FALSE, tolerance = 1e-10)

  # Check if Michael's formula is correct
  # Compute the first term
  A_true <- 0

  for (i in 1:N) {
    for (j in i:N) {
      if (i != j) {
        for (k in 1:K) {
          for (l in 1:K) {
            pi_ij <- pi
            # When D_ij = 0, we must use 1 - pi.
            if (adj[i, j] == 0) {
              pi_ij <- 1 - pi
            }
            A_true <- A_true + tau[i, k]^2 * tau[j, l] * log(pi_ij[k, l]) / (2 * tau[i, k]) + tau[j, l]^2 * tau[i, k] * log(pi_ij[k, l]) / (2 * tau[j, l])
          }
        }
      }
    }
  }

  A <- 0

  for (i in 1:N) {
    for (k in 1:K) {
      for (j in 1:N) {
        if (i != j) {
          for (l in 1:K) {
            pi_ij <- pi
            # When D_ij = 0, we must use 1 - pi.
            if (adj[i, j] == 0) {
              pi_ij <- 1 - pi
            }
            A <- A + tau[i, k]^2 * tau[j, l] * log(pi_ij[k, l]) / (2 * tau[i, k])
          }
        }
      }
    }
  }

  # Check if they are the same
  expect_equal(A, A_true, check.attributes = FALSE, tolerance = 1e-10)
})
