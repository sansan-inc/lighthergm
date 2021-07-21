test_that("getting a subgraph from a network object works", {
  # Which node belongs to which block
  df <-
    tibble::tibble(
      block = c(1, 1, 1, 1, 2, 2, 2, 2),
      node_id = c("A", "B", "C", "D", "E", "F", "G", "H")
    )

  # Edgelist
  edgelist <-
    tibble::tibble(
      source_id = c("A", "A", "A", "B"),
      target_id = c("B", "C", "E", "F")
    )

  # When not all nodes are isolated in a block
  subgraph1 <- get_induced_subgraph(block_structure = df, edgelist = edgelist, searched_block = 1)
  adj_true <- matrix(c(0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), nrow = 4, ncol = 4)
  expect_equal(c("A", "B", "C", "D"), network::network.vertex.names(subgraph1))
  expect_equal(adj_true, network::as.matrix.network.adjacency(subgraph1), check.attributes = FALSE)

  # When all nodes are isolated in a block
  subgraph2 <- get_induced_subgraph(block_structure = df, edgelist = edgelist, searched_block = 2)
  adj_true <- matrix(0, nrow = 4, ncol = 4)
  expect_equal(c("E", "F", "G", "H"), network::network.vertex.names(subgraph2))
  expect_equal(adj_true, network::as.matrix.network.adjacency(subgraph2), check.attributes = FALSE)
})
