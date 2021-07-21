test_that("computing Yule's phi-coefficient works", {
  # Number of nodes
  N <- 1000
  # Number of clusters
  K <- 100

  # An extreme case
  z1 <- rep(1:K, each = N / K)
  z_star <- z1

  # Check if it equals one.
  expect_equal(compute_yule_coef(z_star, z1), 1)

  # Another case
  z2 <- sample(1:K, size = N, replace = TRUE)

  # Compute Yule's phi-coefficient naively.
  n00 <- 0
  n01 <- 0
  n10 <- 0
  n11 <- 0

  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      if (z_star[i] == z_star[j] & z2[i] == z2[j]) {
        n11 <- n11 + 1
      }
      if (z_star[i] == z_star[j] & z2[i] != z2[j]) {
        n10 <- n10 + 1
      }
      if (z_star[i] != z_star[j] & z2[i] == z2[j]) {
        n01 <- n01 + 1
      }
      if (z_star[i] != z_star[j] & z2[i] != z2[j]) {
        n00 <- n00 + 1
      }
    }
  }

  phi_naive <- (n00 * n11 - n01 * n10) / sqrt((n00 + n01) * (n10 + n11) * (n00 + n10) * (n01 + n11))

  # Check if it works
  expect_equal(compute_yule_coef(z_star, z2), phi_naive, tolerance = 1e-10)
})


test_that("Removing missing values works", {
  # Number of nodes
  N <- 1000
  # Number of clusters
  K <- 100

  # An extreme case
  z1 <- rep(1:K, each = N / K)
  z_star <- z1
  z1[1:4] <- NA

  # Check if it equals one.
  compute_yule_coef(z_star, z1)
  expect_silent(compute_yule_coef(z_star, z1))

  z_star[50:56] <- NA
  expect_silent(compute_yule_coef(z_star, z1))
})
