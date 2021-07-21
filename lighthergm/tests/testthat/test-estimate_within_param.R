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
