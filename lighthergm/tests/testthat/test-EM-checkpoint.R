test_that("Setting a checkpoint for EM iterations works", {
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

  ############# 1. Clustering with features ##############################
  # Conduct clustering at once
  initial_weight <- 1000

  cluster_with_feature <-
    lighthergm::hergm(g_sim ~ edges + nodematch("x") + nodematch("y") + triangles,
      n_clusters = K,
      estimate_parameters = TRUE,
      verbose = 0,
      n_em_step_max = 10,
      initialization_method = 3,
      infomap_python = FALSE,
      clustering_with_features = TRUE,
      check_alpha_update = TRUE,
      compute_pi = TRUE,
      check_lower_bound = TRUE,
      check_block_membership = TRUE,
      weight_for_initialization = initial_weight,
      seeds = 334
    )

  # Conduct clustering in two steps
  first_step <-
    lighthergm::hergm(
      object = g_sim ~ edges + nodematch("x") + nodematch("y") + triangles,
      n_clusters = K,
      estimate_parameters = TRUE,
      verbose = 0,
      n_em_step_max = 7,
      initialization_method = 3,
      infomap_python = FALSE,
      clustering_with_features = TRUE,
      check_block_membership = TRUE,
      weight_for_initialization = initial_weight,
      check_alpha_update = TRUE,
      seeds = 334
    )

  second_step <-
    lighthergm::hergm(
      object = first_step,
      n_em_step_max = 3
    )

  # Check if the calculated lower bounds are identical (both length and values)
  expect_equal(cluster_with_feature$EM_lower_bound, second_step$EM_lower_bound)

  # Check if block memberships are identical over iterations (both length and Yule's coefficient)
  expect_true(all(unlist(purrr::map2(cluster_with_feature$EM_list_z, second_step$EM_list_z, compute_yule_coef) == rep(1, length(second_step$EM_list_z)))))

  # The partition should be the same at the end of the second estimation as the one obtained after running 10 iterations of the EM algorithm
  expect_equal(compute_yule_coef(cluster_with_feature$partition, second_step$partition), 1)

  # Check if alphas are identical over iterations
  for (i in 1:length(cluster_with_feature$EM_list_alpha)) {
    expect_equal(cluster_with_feature$EM_list_alpha[[i]], second_step$EM_list_alpha[[i]], check.attribute = FALSE, tolerance = 1e-2)
  }
  # Check if the alpha after the second checkpoint is the same as the alpha when performing estimation with 10 EM iterations without a checkpoint.
  expect_equal(cluster_with_feature$alpha, second_step$alpha)

  # Check if estimated coefficients with and without a checkpoint are identical.
  expect_equal(coef(cluster_with_feature$est_between), coef(second_step$est_between), tolerance = 1e-10)
  expect_equal(coef(cluster_with_feature$est_within), coef(second_step$est_within), tolerance = 1e-10)


  ############# 2. Clustering without features ##############################
  # Conduct clustering at once
  initial_weight <- 1000

  cluster_without_feature <-
    lighthergm::hergm(g_sim ~ edges + nodematch("x") + nodematch("y") + triangles,
      n_clusters = K,
      estimate_parameters = TRUE,
      verbose = 0,
      n_em_step_max = 10,
      initialization_method = 3,
      infomap_python = FALSE,
      clustering_with_features = FALSE,
      check_alpha_update = TRUE,
      compute_pi = TRUE,
      check_lower_bound = TRUE,
      check_block_membership = TRUE,
      weight_for_initialization = initial_weight,
      seeds = 334
    )

  # Conduct clustering in two steps
  first_step_without_feature <-
    lighthergm::hergm(g_sim ~ edges + nodematch("x") + nodematch("y") + triangles,
      n_clusters = K,
      estimate_parameters = TRUE,
      verbose = 0,
      n_em_step_max = 3,
      initialization_method = 3,
      infomap_python = FALSE,
      clustering_with_features = FALSE,
      check_alpha_update = TRUE,
      compute_pi = TRUE,
      check_lower_bound = TRUE,
      check_block_membership = TRUE,
      weight_for_initialization = initial_weight,
      seeds = 334
    )

  second_step_without_feature <-
    lighthergm::hergm(
      object = first_step_without_feature,
      n_em_step_max = 7
    )

  # Check if the calculated lower bounds are identical (both length and values)
  expect_equal(cluster_without_feature$EM_lower_bound, second_step_without_feature$EM_lower_bound)

  # Check if block memberships are identical over iterations (both length and Yule's coefficient)
  expect_true(all(unlist(purrr::map2(cluster_without_feature$EM_list_z, second_step_without_feature$EM_list_z, compute_yule_coef) ==
    rep(1, length(second_step_without_feature$EM_list_z)))))

  # The partition should be the same at the end of the second estimation as the one obtained after running 10 iterations of the EM algorithm
  expect_equal(compute_yule_coef(cluster_without_feature$partition, second_step_without_feature$partition), 1)

  # Check if alphas are identical over iterations
  for (i in 1:length(cluster_without_feature$EM_list_alpha)) {
    expect_equal(cluster_without_feature$EM_list_alpha[[i]], second_step_without_feature$EM_list_alpha[[i]], check.attribute = FALSE, tolerance = 1e-2)
  }
  # Check if the alpha after the second checkpoint is the same as the alpha when performing estimation with 10 EM iterations without a checkpoint.
  expect_equal(cluster_without_feature$alpha, second_step_without_feature$alpha)

  # Check if estimated coefficients with and without a checkpoint are identical.
  expect_equal(coef(cluster_without_feature$est_between), coef(second_step_without_feature$est_between), tolerance = 1e-10)
  expect_equal(coef(cluster_without_feature$est_within), coef(second_step_without_feature$est_within), tolerance = 1e-10)
})
