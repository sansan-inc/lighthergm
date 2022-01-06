test_that("estimating within-block parameters works", {
  #################################################################
  ### Network preparation #########################################
  #################################################################
  # Create an adjacency matrix
  adj <- c(
    c(0, 1, 0, 1, 1, 0),
    c(1, 0, 1, 0, 0, 1),
    c(0, 1, 0, 1, 1, 0),
    c(1, 0, 1, 0, 1, 1),
    c(1, 0, 1, 1, 0, 1),
    c(0, 1, 0, 1, 1, 0)
  )
  adj <- matrix(data = adj, nrow = 6, ncol = 6)
  rownames(adj) <- as.character(1001:1006)
  colnames(adj) <- as.character(1001:1006)

  # Vertex attribute
  x <- c(1, 0, 0, 1, 1, 0)

  # Block
  block <- c(1, 2, 3, 1, 3, 2)

  g <- network::network(adj, matrix.type = "adjacency")
  network::set.vertex.attribute(g, attrname = "x", value = x)

  #################################################################
  ### Estimate the within-block parameters ########################
  #################################################################
  suppressWarnings(est <- estimate_within_params(
    formula = g ~ edges + nodematch("x"),
    network = g,
    z_memb = block,
    parallel = FALSE,
    verbose = 0,
    initial_estimate = NULL,
    seeds = NULL,
    method_second_step = "MPLE"
  ))


  # Get the network used for estimation
  g_est <- est$network

  # Get the adjacency matrix for the network
  adj_est <- network::as.matrix.network.adjacency(g_est)

  #################################################################
  ### Test if the created and true networks are the same  #########
  #################################################################

  # Prepare the answer
  adj_ans <- c(
    c(0, 1, 0, 0, 0, 0),
    c(1, 0, 0, 0, 0, 0),
    c(0, 0, 0, 1, 0, 0),
    c(0, 0, 1, 0, 0, 0),
    c(0, 0, 0, 0, 0, 1),
    c(0, 0, 0, 0, 1, 0)
  )

  # Order of vertex id after diagonization: 1001, 1004, 1002, 1006, 1003, 1005
  # Order of original x: c(1, 0, 0, 1, 1, 0)
  vertex_id_ans <- as.character(c(1001, 1004, 1002, 1006, 1003, 1005))
  x_ans <- c(1, 1, 0, 0, 0, 1)

  # Test
  expect_equal(adj_est, adj_ans, check.attributes = FALSE)
  expect_equal(network::get.vertex.attribute(g_est, "x"), x_ans)
  expect_equal(network::network.vertex.names(g_est), vertex_id_ans)
})

test_that("estimating within-block parameters works with non-consecutive block names", {
  adj <- c(
    c(0, 1, 0, 1, 1, 0),
    c(1, 0, 1, 0, 0, 1),
    c(0, 1, 0, 1, 1, 0),
    c(1, 0, 1, 0, 1, 1),
    c(1, 0, 1, 1, 0, 1),
    c(0, 1, 0, 1, 1, 0)
  )
  adj <- matrix(data = adj, nrow = 6, ncol = 6)
  rownames(adj) <- as.character(1001:1006)
  colnames(adj) <- as.character(1001:1006)

  x <- c(1, 0, 0, 1, 1, 0)

  # Use non-consecutive block names
  block <- c(50, 70, 95, 50, 95, 70)

  g <- network::network(adj, matrix.type = "adjacency")
  network::set.vertex.attribute(g, attrname = "x", value = x)

  suppressWarnings(est <- estimate_within_params(
    formula = g ~ edges + nodematch("x"),
    network = g,
    z_memb = block,
    parallel = FALSE,
    verbose = 0,
    initial_estimate = NULL,
    seeds = NULL,
    method_second_step = "MPLE"
  ))


  # Get the network used for estimation
  g_est <- est$network

  # Get the adjacency matrix for the network
  adj_est <- network::as.matrix.network.adjacency(g_est)

  adj_ans <- c(
    c(0, 1, 0, 0, 0, 0),
    c(1, 0, 0, 0, 0, 0),
    c(0, 0, 0, 1, 0, 0),
    c(0, 0, 1, 0, 0, 0),
    c(0, 0, 0, 0, 0, 1),
    c(0, 0, 0, 0, 1, 0)
  )

  vertex_id_ans <- as.character(c(1001, 1004, 1002, 1006, 1003, 1005))
  x_ans <- c(1, 1, 0, 0, 0, 1)

  # Check that the network is the same
  expect_equal(adj_est, adj_ans, check.attributes = FALSE)
  expect_equal(network::get.vertex.attribute(g_est, "x"), x_ans)
  expect_equal(network::network.vertex.names(g_est), vertex_id_ans)

  # Check that the blocks are assigned to the right nodes
  g_nodes_data <- data.frame(
    id = network::network.vertex.names(g),
    block = block
  ) %>% dplyr::arrange(id)

  est_g_nodes_data <- data.frame(
    id = network::network.vertex.names(g_est),
    block = as.double(network::get.vertex.attribute(g_est, 'block'))
  ) %>% dplyr::arrange(id)

  expect_equal(g_nodes_data$id, est_g_nodes_data$id)
  expect_equal(g_nodes_data$block, est_g_nodes_data$block)
})

test_that("control.ergm settings can be passed to the within estimation from hergm function", {
  # Define some settings for testing purposes
  test_burnin <- 9797
  test_interval <- 3434
  test_method <- 'Stepping'

  hergm_formula <- g ~ edges + triangle + nodematch("x")

  n_nodes <- 100
  n_clusters <- 2

  nodes_data <- tibble::tibble(
    node_id = 1:n_nodes,
    x = sample(1:2, size = n_nodes, replace = T),
    block = sample(1:n_clusters, size = n_nodes, replace = T)
  )

  g <- network::network.initialize(n = n_nodes)
  network::set.vertex.attribute(g, "x", nodes_data$x)
  list_feature_matrices <- lighthergm::get_list_sparse_feature_adjmat(g, hergm_formula)

  coef_between_block <- c(-3, 1)
  coef_within_block <- c(-2, 0.1, 0.5)

  sim_ergm_control <- ergm::control.simulate.formula(
    MCMC.burnin = 4000000,
    MCMC.interval = 200000
  )

  g <- lighthergm::simulate_hergm(
    formula_for_simulation = hergm_formula,
    data_for_simulation = nodes_data,
    colname_vertex_id = "node_id",
    colname_block_membership = "block",
    coef_between_block = coef_between_block,
    coef_within_block = coef_within_block,
    ergm_control = sim_ergm_control,
    fast_between_simulation = TRUE,
    list_feature_matrices = list_feature_matrices
  )

  hergm_res <- lighthergm::hergm(
    g ~ edges + nodematch("x") + triangle,
    n_clusters = n_clusters,
    n_em_step_max = 10,
    estimate_parameters = T,
    clustering_with_features = T,
    method_second_step = 'MLE',
    control = ergm::control.ergm(
      MCMC.burnin = test_burnin,
      MCMC.interval = test_interval,
      main.method = test_method
    )
  )

  used_control <- hergm_res$est_within$control

  expect_equal(used_control$MCMC.burnin, test_burnin)
  expect_equal(used_control$MCMC.interval, test_interval)
  expect_equal(used_control$main.method, test_method)
})
