test_that("linear term computation without features works", {
  # Number of nodes
  N <- 12
  # Number of clusters
  K <- 3

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

  # Compute gamma (parameter of multinomial distribution)
  alpha <- colSums(tau)

  # Compute the true linear term in a naive way
  s_true <- matrix(0, nrow = N, ncol = K)

  for (i in 1:N) {
    for (k in 1:K) {
      s_ik <- 1 + log(alpha[k]) - log(tau[i, k])
      s_true[i, k] <- s_true[i, k] + s_ik
    }
  }

  s <- compute_linear_term(N, K, alpha, tau, 0)
  expect_equal(s, s_true, check.attributes = FALSE, tolerance = 1e-10)
})
