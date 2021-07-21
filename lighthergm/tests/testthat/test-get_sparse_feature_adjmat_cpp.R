test_that("creating a feature adjacency matrix works", {
  library(magrittr)

  get_sparse_feature_adjacency_matrix <- function(feature) {
    S <- Matrix::sparseMatrix(i = {}, j = {}, dims = c(length(feature), length(feature)))
    S <- as(S, "dgCMatrix")
    for (i in 1:length(feature)) {
      for (j in 1:length(feature)) {
        if (i != j) {
          if (feature[i] == feature[j]) {
            S[i, j] <- 1
          }
        }
      }
    }
    return(S)
  }

  # Number of nodes
  N <- 100

  # Features
  x <- sample(x = c(1:100), size = N, replace = TRUE)
  y <- sample(x = c(1:30), size = N, replace = TRUE)
  z <- sample(x = c(1:50), size = N, replace = TRUE)
  w <- sample(x = c(LETTERS, letters), size = N, replace = TRUE)

  # Create an edgelist
  edgelist <-
    tibble::tibble(tail = 1:N, head = 1:N) %>%
    tidyr::expand(tail, head) %>%
    dplyr::filter(tail < head) %>%
    dplyr::mutate(connect = unlist(as.integer(purrr::rbernoulli(n = nrow(.), p = 0.005)))) %>%
    dplyr::filter(connect == 1) %>%
    dplyr::select(tail, head)

  # Create a network object
  g <- network::network.initialize(n = N, directed = FALSE)
  network::add.edges(x = g, tail = edgelist$tail, head = edgelist$head)
  network::set.vertex.attribute(x = g, attrname = "x", value = x)
  network::set.vertex.attribute(x = g, attrname = "y", value = y)
  network::set.vertex.attribute(x = g, attrname = "z", value = z)
  network::set.vertex.attribute(x = g, attrname = "w", value = w)

  # Create a formula
  form <- g ~ edges + triangles + nodematch("x") + nodematch("y") + nodematch("z") + nodematch("w")

  # True list of feature adjacency metrices
  list_adjmat_true <- list(
    get_sparse_feature_adjacency_matrix(x),
    get_sparse_feature_adjacency_matrix(y),
    get_sparse_feature_adjacency_matrix(z),
    get_sparse_feature_adjacency_matrix(w)
  )

  # Create a list
  list_adjmat <- get_list_sparse_feature_adjmat(network = g, formula = form)

  # Check if it works
  for (i in 1:4) {
    expect_equal(list_adjmat[[i]], list_adjmat_true[[i]], check.attributes = FALSE, tolerance = 1e-10)
  }
})
