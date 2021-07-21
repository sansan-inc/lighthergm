test_that("element-wise matrix multiplication works", {
  # Prepare matrices
  N <- 50
  G <- matrix(as.integer(unlist(purrr::rbernoulli(n = N * N, p = 0.2))), nrow = N, ncol = N)
  diag(G) <- 0
  G <- as(G, "sparseMatrix")
  G <- Matrix::forceSymmetric(G)

  S <- matrix(as.integer(unlist(purrr::rbernoulli(n = N * N, p = 0.2))), nrow = N, ncol = N)
  diag(S) <- 0
  S <- as(S, "sparseMatrix")
  S <- Matrix::forceSymmetric(S)

  V <- matrix(as.integer(unlist(purrr::rbernoulli(n = N * N, p = 0.2))), nrow = N, ncol = N)
  diag(V) <- 0
  V <- as(V, "sparseMatrix")
  V <- Matrix::forceSymmetric(V)

  # Prepare true results
  output_true <- list()
  # The first element of the list is filled with - (S + V) + (S % V), which will be used to compute pi_d0x0.
  output_true[[1]] <- -(S + V) + (S * V)
  # N = 1: (1, 0, 0)
  output_true[[2]] <- Matrix::drop0(G * (1 - S) * (1 - V))
  # N = 2: (0, 1, 0)
  output_true[[3]] <- Matrix::drop0((1 - G) * S * (1 - V))
  # N = 3: (1, 1, 0)
  output_true[[4]] <- Matrix::drop0(G * S * (1 - V))
  # N = 4: (0, 0, 1)
  output_true[[5]] <- Matrix::drop0((1 - G) * (1 - S) * V)
  # N = 5: (1, 0, 1)
  output_true[[6]] <- Matrix::drop0(G * (1 - S) * V)
  # N = 6: (0, 1, 1)
  output_true[[7]] <- Matrix::drop0((1 - G) * S * V)
  # N = 7: (1, 1, 1)
  output_true[[8]] <- Matrix::drop0(G * S * V)

  # Compute element-wise mutiplied matrices
  output <- get_elementwise_multiplied_matrices(G, list(S, V))

  for (i in 1:8) {
    expect_equal(output[[i]], output_true[[i]], check.attributes = FALSE, )
  }
})
