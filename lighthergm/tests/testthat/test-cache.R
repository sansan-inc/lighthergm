set.seed(334)

# Simulate a random network for testing
simulate_network <- function(){
  N <- 1000
  K <- 50
  memb <- rep(1:K, each = N / K)
  x <- sample(1:10, size = N, replace = TRUE)
  y <- sample(1:10, size = N, replace = TRUE)
  list_within_params <- c(-1, 1, 1, 0.5)
  list_between_params <- c(-3.5, 0.5, 0.5)

  formula <- g ~ edges + nodematch("x") + nodematch("y") + triangle

  vertex_id <- 1:N

  df <- tibble::tibble(
    id = vertex_id,
    memb = memb,
    x = x,
    y = y
  )

  simulate_hergm(
      formula_for_simulation = formula,
      data_for_simulation = df,
      colname_vertex_id = "id",
      colname_block_membership = "memb",
      coef_within_block = list_within_params,
      coef_between_block = list_between_params,
      ergm_control = ergm::control.simulate.formula(MCMC.burnin = 1000000, MCMC.interval = 1000),
      seed = 1,
      n_sim = 1,
      directed = FALSE,
      output = "network"
    )
}

# Function that checks that the number of files in the directory is the expected number
check_files <- function(dir, expected_number){
  expect_equal(length(list.files(dir, '*.rds')), expected_number)
}

cleanup <- function(){
  do.call(file.remove, list(list.files(tempdir(), '*.rds', full.names = TRUE)))
}

test_that('Estimation with a disk cache stores data in the correct directory', {
  on.exit(cleanup())

  # Use this directory for caching in disk
  dir <- tempdir()

  # Simulate a network
  g_1 <- simulate_network()

  # There should be no cached files in the directory
  check_files(dir, 0)

  # Perform estimation
  lighthergm::hergm(
    object = g_1 ~ edges + nodematch("x"),
    n_clusters = 20,
    n_em_step_max = 3,
    initialization_method = 1,
    clustering_with_features = TRUE,
    verbose=2,
    cache = cachem::cache_disk(dir)
  )

  # The estimation should have stored one RDS object in the cache directory.
  check_files(dir, 1)

  lighthergm::hergm(
    object = g_1 ~ edges + nodematch("x"),
    n_clusters = 20,
    n_em_step_max = 3,
    initialization_method = 1,
    clustering_with_features = TRUE,
    verbose=2,
    cache = cachem::cache_disk(dir)
  )

  # Running again the estimation on the same network should reuse the previously stored RDS object
  # and not store a new one.
  check_files(dir, 1)

  # Generate a different network
  g_2 <- simulate_network()

  # Perform estimation on the new network.
  lighthergm::hergm(
    object = g_2 ~ edges + nodematch("x"),
    n_clusters = 20,
    n_em_step_max = 3,
    initialization_method = 1,
    clustering_with_features = TRUE,
    verbose=2,
    cache = cachem::cache_disk(dir)
  )

  # The network has changed, so the previously cached RDS is not reused and a new cache file is generated.
  check_files(dir, 2)

})
