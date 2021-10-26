test_that("combining within- and between-block edgelists while removing duplicated edges works", {
  edgelist_within <- matrix(c(1, 1, 2, 3, 2, 3, 4, 5), nrow = 4, ncol = 2)
  edgelist_between <- matrix(c(1, 2, 3, 4, 3, 5), nrow = 3, ncol = 2)

  # Create a true edgelist
  true_edgelist <-
    data.frame(
      tail = c(1, 1, 2, 3, 1, 2, 3),
      head = c(2, 3, 4, 5, 4, 3, 5)
    ) %>%
    dplyr::distinct(tail, head) %>%
    dplyr::arrange(tail) %>%
    as.matrix()

  # Get edgelist from R function
  edgelist <- combine_within_between_edges(edgelist_within, edgelist_between, use_fast_between_simulation = FALSE)

  # Check if it works
  expect_equal(true_edgelist, unclass(edgelist), check.attributes = FALSE)
})

test_that("correctly attaching vertex ids, block memberships, and vertex features to the simulated network works", {
  set.seed(1)
  # Prepare ingredients for simulating a network
  N <- 1000
  K <- 10

  list_within_params <- c(-3, 1, 1, 0.76, 0.08)
  list_between_params <- c(-5, 2, 2)
  formula <- g ~ edges + nodematch("x") + nodematch("y") + triangle + kstar(2)

  memb <- sample(1:K, size = N, replace = TRUE)
  vertex_id <- as.character(11:(11 + N - 1))

  x <- sample(1:20, size = N, replace = TRUE)
  y <- sample(LETTERS, size = N, replace = TRUE)


  df <- tibble::tibble(
    id = vertex_id,
    memb = memb,
    x = x,
    y = y
  )

  # Simulate a network
  g_sim <-
    simulate_hergm(
      formula_for_simulation = formula,
      data_for_simulation = df,
      colname_vertex_id = "id",
      colname_block_membership = "memb",
      coef_within_block = list_within_params,
      coef_between_block = list_between_params,
      ergm_control = ergm::control.simulate.formula(),
      seed_for_within = 1,
      seed_for_between = 1,
      n_sim = 1,
      directed = FALSE
    )

  # Prepare a dataframe sorted by block
  df_sorted <-
    df %>%
    dplyr::arrange(memb)

  # Check if vertex ids, block memberships, and vertex features are stored correctly in the network object.
  expect_equal(df_sorted$id, network::get.vertex.attribute(g_sim, "vertex.names"), check.attributes = FALSE)
  expect_equal(df_sorted$memb, network::get.vertex.attribute(g_sim, "block"), check.attributes = FALSE)
  expect_equal(df_sorted$x, network::get.vertex.attribute(g_sim, "x"), check.attributes = FALSE)
  expect_equal(df_sorted$y, network::get.vertex.attribute(g_sim, "y"), check.attributes = FALSE)
})

test_that("simulating a network from a given edgelist works", {
  set.seed(1)
  # Prepare ingredients for simulating a network
  N <- 1000
  K <- 10

  list_within_params <- c(-3, 1, 1, 0.76, 0.08)
  list_between_params <- c(-5, 2, 2)
  formula <- g ~ edges + nodematch("x") + nodematch("y") + triangle + kstar(2)

  memb <- sample(1:K, size = N, replace = TRUE)
  vertex_id <- as.character(11:(11 + N - 1))

  x <- sample(1:20, size = N, replace = TRUE)
  y <- sample(LETTERS, size = N, replace = TRUE)


  df <- tibble::tibble(
    id = vertex_id,
    memb = memb,
    x = x,
    y = y
  )

  # Simulate a network
  g_sim <-
    simulate_hergm(
      formula_for_simulation = formula,
      data_for_simulation = df,
      colname_vertex_id = "id",
      colname_block_membership = "memb",
      coef_within_block = list_within_params,
      coef_between_block = list_between_params,
      ergm_control = ergm::control.simulate.formula(),
      seed_for_within = 1,
      seed_for_between = 1,
      n_sim = 1,
      directed = FALSE
    )

  # Convert to an edgelist
  g_sim <- network::as.edgelist(g_sim)

  # Simulate a network from a given edgelist
  g2 <-
    simulate_hergm(
      formula_for_simulation = formula,
      data_for_simulation = df,
      colname_vertex_id = "id",
      colname_block_membership = "memb",
      coef_within_block = list_within_params,
      coef_between_block = list_between_params,
      ergm_control = ergm::control.simulate.formula(),
      seed_for_within = 1,
      seed_for_between = 1,
      seed_edgelist = g_sim,
      n_sim = 1,
      directed = FALSE
    )

  # Check if the network is correctly generated
  expect_match(class(g2), "network")
})

test_that("generating multiple networks works", {
  set.seed(1)
  # Prepare ingredients for simulating a network
  N <- 1000
  K <- 10

  list_within_params <- c(-3, 1, 1, 0.76, 0.08)
  list_between_params <- c(-5, 2, 2)
  formula <- g ~ edges + nodematch("x") + nodematch("y") + triangle + kstar(2)

  memb <- sample(1:K, size = N, replace = TRUE)
  vertex_id <- as.character(11:(11 + N - 1))

  x <- sample(1:20, size = N, replace = TRUE)
  y <- sample(LETTERS, size = N, replace = TRUE)


  df <- tibble::tibble(
    id = vertex_id,
    memb = memb,
    x = x,
    y = y
  )

  # Simulate a network
  g_sim <-
    simulate_hergm(
      formula_for_simulation = formula,
      data_for_simulation = df,
      colname_vertex_id = "id",
      colname_block_membership = "memb",
      coef_within_block = list_within_params,
      coef_between_block = list_between_params,
      ergm_control = ergm::control.simulate.formula(),
      seed_for_within = 1,
      seed_for_between = 1,
      n_sim = 3,
      directed = FALSE
    )

  expect_equal(length(g_sim), 3)
  expect_match(class(g_sim[[1]]), "network")
  expect_match(class(g_sim[[2]]), "network")
  expect_match(class(g_sim[[3]]), "network")
})

test_that("generating stats and avoiding generating within-block links while simulating a between-block network works", {
  set.seed(1)
  # Prepare ingredients for simulating a network
  N <- 1000
  K <- 10

  list_within_params <- c(-3, 1, 1, 0.76, 0.08)
  list_between_params <- c(-5, 2, 2)
  formula <- g ~ edges + nodematch("x") + nodematch("y") + triangle + kstar(2)

  memb <- sample(1:K, size = N, replace = TRUE)
  vertex_id <- as.character(11:(11 + N - 1))

  x <- sample(1:20, size = N, replace = TRUE)
  y <- sample(LETTERS, size = N, replace = TRUE)


  df <- tibble::tibble(
    id = vertex_id,
    memb = memb,
    x = x,
    y = y
  )

  # Simulate a network
  g_sim <-
    simulate_hergm(
      formula_for_simulation = formula,
      data_for_simulation = df,
      colname_vertex_id = "id",
      colname_block_membership = "memb",
      coef_within_block = list_within_params,
      coef_between_block = list_between_params,
      ergm_control = ergm::control.simulate.formula(MCMC.burnin = 1000, MCMC.interval = 1000),
      seed_for_within = 1,
      seed_for_between = 1,
      n_sim = 10,
      directed = FALSE,
      output = "stats",
      prevent_duplicate = TRUE
    )

  # Check whether network stats are generated.
  expect_equal(c("edges", "nodematch.x", "nodematch.y", "triangle", "kstar2"), colnames(g_sim[[1]]))
  expect_equal(c("edges", "nodematch.x", "nodematch.y", "nodematch.block"), colnames(g_sim[[2]]))

  # Check whther no within-block links are generated in the between-block network.
  expect_equal(g_sim[[2]][, 4], rep(0, 10))
})


test_that("simulating between-block networks using the cpp function works", {
  set.seed(1)
  # Prepare ingredients for simulating a network
  N <- 1000
  K <- 10

  list_within_params <- c(-3, 1, 1, 0.76, 0.08)
  list_between_params <- c(-5, 2, 2)
  formula <- g ~ edges + nodematch("x") + nodematch("y") + triangle + kstar(2)

  memb <- sample(1:K, size = N, replace = TRUE)
  vertex_id <- as.character(11:(11 + N - 1))

  x <- sample(1:20, size = N, replace = TRUE)
  y <- sample(LETTERS, size = N, replace = TRUE)

  # Prepare a list of feature adjacency matrices
  X <- get_sparse_feature_adjmat(x)
  Y <- get_sparse_feature_adjmat_from_string(y)
  list_feature_adjmat <- list(X, Y)

  df <- tibble::tibble(
    id = vertex_id,
    memb = memb,
    x = x,
    y = y
  )

  # When use_fast_between_simulation = TRUE but list_feature_matrices is not given, it should yield an error.
  expect_error(g_sim <-
    simulate_hergm(
      formula_for_simulation = formula,
      data_for_simulation = df,
      colname_vertex_id = "id",
      colname_block_membership = "memb",
      coef_within_block = list_within_params,
      coef_between_block = list_between_params,
      ergm_control = ergm::control.simulate.formula(MCMC.burnin = 1000, MCMC.interval = 1000),
      seed_for_within = 1,
      seed_for_between = 1,
      n_sim = 1,
      directed = FALSE,
      output = "edgelist",
      use_fast_between_simulation = TRUE
    ))

  # When simulating only one network:
  g_sim <-
    simulate_hergm(
      formula_for_simulation = formula,
      data_for_simulation = df,
      colname_vertex_id = "id",
      colname_block_membership = "memb",
      coef_within_block = list_within_params,
      coef_between_block = list_between_params,
      ergm_control = ergm::control.simulate.formula(MCMC.burnin = 1000, MCMC.interval = 1000),
      seed_for_within = 1,
      seed_for_between = 1,
      n_sim = 1,
      directed = FALSE,
      output = "network",
      use_fast_between_simulation = TRUE,
      list_feature_matrices = list_feature_adjmat,
      verbose = 1
    )

  # Prepare a dataframe sorted by block
  df_sorted <-
    df %>%
    dplyr::arrange(memb)

  # Check if vertex ids, block memberships, and vertex features are stored correctly in the network object.
  expect_equal(df_sorted$id, network::get.vertex.attribute(g_sim, "vertex.names"), check.attributes = FALSE)
  expect_equal(df_sorted$memb, network::get.vertex.attribute(g_sim, "block"), check.attributes = FALSE)
  expect_equal(df_sorted$x, network::get.vertex.attribute(g_sim, "x"), check.attributes = FALSE)
  expect_equal(df_sorted$y, network::get.vertex.attribute(g_sim, "y"), check.attributes = FALSE)


  # Check if multiple networks can be generated.
  g_sim <-
    simulate_hergm(
      formula_for_simulation = formula,
      data_for_simulation = df,
      colname_vertex_id = "id",
      colname_block_membership = "memb",
      coef_within_block = list_within_params,
      coef_between_block = list_between_params,
      ergm_control = ergm::control.simulate.formula(MCMC.burnin = 1000, MCMC.interval = 1000),
      seed_for_within = 1,
      seed_for_between = NULL,
      n_sim = 3,
      directed = FALSE,
      output = "network",
      use_fast_between_simulation = TRUE,
      list_feature_matrices = list_feature_adjmat,
      verbose = 1
    )

  expect_equal(length(g_sim), 3)
  expect_match(class(g_sim[[1]]), "network")
  expect_match(class(g_sim[[2]]), "network")
  expect_match(class(g_sim[[3]]), "network")

  # Check if the function works when output = "stats" and use_fast_between_simulation = TRUE.
  expect_error(
    g_sim <-
      simulate_hergm(
        formula_for_simulation = formula,
        data_for_simulation = df,
        colname_vertex_id = "id",
        colname_block_membership = "memb",
        coef_within_block = list_within_params,
        coef_between_block = list_between_params,
        ergm_control = ergm::control.simulate.formula(MCMC.burnin = 10000, MCMC.interval = 10000),
        seed_for_within = 1,
        seed_for_between = NULL,
        n_sim = 3,
        directed = FALSE,
        output = "stats",
        use_fast_between_simulation = TRUE,
        list_feature_matrices = list_feature_adjmat
      ),
    NA
  )
})

test_that("Yielding the same between-block network using the cpp function when the seed is fixed works", {
  # Number of nodes
  N <- 1000

  # Number of blocks
  K <- 50

  # Block membership
  block <- rep(1:K, each = N / K)

  # Parameters
  param_edges <- -1
  covar_param <- c(2, 2)
  coef_between <- c(param_edges, covar_param)

  # Covariates
  x <- sample(1:10, size = N, replace = TRUE)
  y <- sample(1:10, size = N, replace = TRUE)
  z <- sample(1:10, size = N, replace = TRUE)

  # Initialize a sparse adjacency matrix
  G <- Matrix::sparseMatrix(i = {}, j = {}, dims = c(N, N))

  # Prepare feature adjacency matrices
  X <- get_sparse_feature_adjmat(x)
  Y <- get_sparse_feature_adjmat(y)
  Z <- get_sparse_feature_adjmat(z)
  list_feature_adjmat <- list(X, Y, Z)

  # Simulate three between-block networks.
  set.seed(334)
  adj1 <- simulate_between_network(N, list_feature_adjmat, coef_between, block, directed = FALSE)
  set.seed(334)
  adj2 <- simulate_between_network(N, list_feature_adjmat, coef_between, block, directed = FALSE)
  set.seed(334)
  adj3 <- simulate_between_network(N, list_feature_adjmat, coef_between, block, directed = FALSE)

  # Check if all the adjacency matrices are identical.
  expect_equal(adj1, adj2, check.attributes = FALSE)
  expect_equal(adj1, adj3, check.attributes = FALSE)
  expect_equal(adj2, adj3, check.attributes = FALSE)

  # If seed is not given, different adjacency matrices must be generated.
  # Simulate three between-block networks.
  set.seed(NULL)
  adj1 <- simulate_between_network(N, list_feature_adjmat, coef_between, block, directed = FALSE)
  adj2 <- simulate_between_network(N, list_feature_adjmat, coef_between, block, directed = FALSE)
  adj3 <- simulate_between_network(N, list_feature_adjmat, coef_between, block, directed = FALSE)

  # Check if they are all different.
  expect_false(all(adj1 == adj2))
  expect_false(all(adj1 == adj3))
  expect_false(all(adj2 == adj3))
})


test_that("Simulating networks using formula without externality terms works", {
  K <- 10
  N <- 200
  formula_1 <- g ~ edges + nodematch("x")
  x <- sample(c(1, 2), size = N, replace = TRUE)
  block <- sample(1:K, size = N, replace = TRUE)

  nodes_data <-
    tibble::tibble(x = x, node_id = as.character(1:N), block = block)

  coef_between <- c(-6.5, 0)
  coef_within <- c(-4.5, 5)

  expect_error(
    sim_net_1 <- lighthergm::simulate_hergm(
      formula_for_simulation = formula_1,
      data_for_simulation = nodes_data,
      colname_vertex_id = 'node_id',
      colname_block_membership = 'block',
      coef_between_block = coef_between,
      coef_within_block = coef_within,
      n_sim = 1
    ),
    NA
  )
})

test_that("Simulating networks using formula without covariates works", {
  K <- 10
  N <- 200
  formula_1 <- g ~ edges
  x <- sample(c(1, 2), size = N, replace = TRUE)
  block <- sample(1:K, size = N, replace = TRUE)

  nodes_data <-
    tibble::tibble(x = x, node_id = as.character(1:N), block = block)

  coef_between <- c(-6.5)
  coef_within <- c(-4.5)

  expect_error(
    sim_net_1 <- lighthergm::simulate_hergm(
      formula_for_simulation = formula_1,
      data_for_simulation = nodes_data,
      colname_vertex_id = 'node_id',
      colname_block_membership = 'block',
      coef_between_block = coef_between,
      coef_within_block = coef_within,
      n_sim = 1
    ),
    NA
  )
})
