test_that("starting EM iterations and parameter estimation from a given vector of block memberships works", {
  set.seed(334)
  # Simulate a network to work with in this unit test.
  # Number of nodes
  N <- 1000
  # Number of blocks
  K <- 50
  # Block memberships (same block size)
  memb <- rep(1:K, each = N / K)
  # Covariates
  x <- sample(1:10, size = N, replace = TRUE)
  y <- sample(1:10, size = N, replace = TRUE)

  # Within-block parameters: edges, nodematch("x"), nodematch("y"), triangle
  list_within_params <- c(-1, 1, 1, 0.5)
  # Between-block parameters: edges, nodematch("x"), nodematch("y")
  list_between_params <- c(-3.5, 0.5, 0.5)

  formula <- g ~ edges + nodematch("x") + nodematch("y") + triangle

  vertex_id <- 1:N

  df <- tibble::tibble(
    id = vertex_id,
    memb = memb,
    x = x,
    y = y
  )

  g_sim <-
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

  # Conduct clustering
  cluster_with_feature <-
    lighthergm::hergm(g_sim ~ edges + nodematch("x") + nodematch("y") + triangles,
      n_clusters = K,
      estimate_parameters = FALSE,
      verbose = 0,
      n_em_step_max = 3,
      initialization_method = 3,
      infomap_python = FALSE,
      clustering_with_features = TRUE,
      check_alpha_update = TRUE,
      compute_pi = TRUE,
      check_lower_bound = TRUE,
      check_block_membership = TRUE,
      seeds = 334
    )

  # Check if starting from the previously estimated block memberships works.
  expect_error(result <-
    lighthergm::hergm(g_sim ~ edges + nodematch("x") + nodematch("y") + triangles,
      initialized_cluster_data = cluster_with_feature$partition,
      n_em_step_max = 2,
      estimate_parameters = FALSE
    ), NA)

  # Check if starting from block memberships initialized Python's infomap works.
  expect_error(result2 <-
    lighthergm::hergm(g_sim ~ edges + nodematch("x") + nodematch("y") + triangles,
      initialized_cluster_data = system.file("extdata", "initialized_cluster_data_by_infomap.clu", package = "lighthergm"),
      n_em_step_max = 1,
      estimate_parameters = FALSE,
      verbose = 1
    ), NA)

  # Check if starting paramter estimation from a given vector of block memberships works.
  expect_error(result3 <-
    lighthergm::hergm(g_sim ~ edges + nodematch("x") + nodematch("y") + triangles,
      block_membership = result$partition,
      verbose = 1
    ), NA)

  # Check if not specifying n_clusters when initialized_cluster_data and block_membership are null yields an error.
  expect_error(result4 <-
    lighthergm::hergm(g_sim ~ edges + nodematch("x") + nodematch("y") + triangles,
      verbose = 1
    ))
})
