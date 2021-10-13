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

  # Obtain the stats
  within_sim_stats <- lighthergm::simulate_hergm_within(
    formula_for_simulation = formula,
    data_for_simulation = df,
    colname_vertex_id = 'id',
    colname_block_membership = 'memb',
    coef_within_block = list_within_params,
    ergm_control = ergm::control.simulate.formula(),
    seed = 1,
    n_sim = 3
  )

  expected_terms <- statnet.common::list_rhs.formula(formula)

  expect_equal(nrow(within_sim_stats), 3)
  expect_equal(length(expected_terms), length(names(within_sim_stats)))

})
