get_dummy_net <- function(n_nodes, n_clusters, em_iters = 10) {
  hergm_formula <- g ~ edges + triangle + nodematch("x")

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
    n_em_step_max = em_iters,
    estimate_parameters = T,
    clustering_with_features = T
  )

  list(
    hergm_res = hergm_res,
    g = g,
    nodes_data = nodes_data,
    K = n_clusters,
    list_feature_matrices = list_feature_matrices,
    vertex_id_var = "node_id",
    block_id_var = "block",
    ergm_control = sim_ergm_control
  )
}


test_that("Returned GOF dataframe has the correct fields", {
  sim <- get_dummy_net(50, 2)
  g <- sim$g

  test_gof_res <- lighthergm::gof_lighthergm(
    g,
    list_feature_matrices = sim$list_feature_matrices,
    data_for_simulation = sim$nodes_data,
    colname_vertex_id = sim$vertex_id_var,
    colname_block_membership = sim$block_id_var,
    lighthergm_results = sim$hergm_res,
    ergm_control = sim$ergm_control,
    n_sim = 3
  )

  for (stat_type in c("original", "simulated")) {
    stats <- test_gof_res[[stat_type]]
    expect_false(is.null(stats))
    for (stat in c("network_stats", "degree_dist", "esp_dist")) {
      expect_false(is.null(stats[[stat]]))
    }
    expect_true(is.null(stats[["geodesic_dist"]]))
  }
})

test_that("GOF network stats have the right fields and terms", {
  sim <- get_dummy_net(50, 2)
  g <- sim$g

  test_gof_res <- lighthergm::gof_lighthergm(
    g,
    list_feature_matrices = sim$list_feature_matrices,
    data_for_simulation = sim$nodes_data,
    colname_vertex_id = sim$vertex_id_var,
    colname_block_membership = sim$block_id_var,
    lighthergm_results = sim$hergm_res,
    ergm_control = sim$ergm_control,
    n_sim = 3
  )

  expected_terms <- ergm::ergm_model(sim$hergm_res$est_within$formula)$terms %>%
    purrr::map(function(t) {
      `$`(t, name)
    })
  unlist

  for (stat_type in c("original", "simulated")) {
    stat_type_df <- test_gof_res[[stat_type]]
    actual_terms <- colnames(stat_type_df$network_stats)

    actual_terms[stringr::str_detect(actual_terms, "n_sim", negate = TRUE)] %>%
      setdiff(c("value", "stat")) %>%
      length() %>%
      expect_equal(0)

    stat_type_df$network_stats$stat %>%
      unique() %>%
      stringr::str_replace("[.].*", "") %>%
      setdiff(expected_terms) %>%
      length() %>%
      expect_equal(0)
  }
})

test_that("GOF degree stats have the right fields and terms", {
  sim <- get_dummy_net(50, 2)
  g <- sim$g

  test_gof_res <- lighthergm::gof_lighthergm(
    g,
    list_feature_matrices = sim$list_feature_matrices,
    data_for_simulation = sim$nodes_data,
    colname_vertex_id = sim$vertex_id_var,
    colname_block_membership = sim$block_id_var,
    lighthergm_results = sim$hergm_res,
    ergm_control = sim$ergm_control,
    n_sim = 3
  )

  for (stat_type in c("original", "simulated")) {
    stat_type_df <- test_gof_res[[stat_type]]
    actual_terms <- colnames(stat_type_df$degree_dist)

    actual_terms[stringr::str_detect(actual_terms, "n_sim", negate = TRUE)] %>%
      setdiff(c("degree", "share")) %>%
      length() %>%
      expect_equal(0)

    expect_lte(max(stat_type_df$degree_dist$degree), g$gal$n)
    expect(
      min(stat_type_df$degree_dist$share) >= 0 && max(stat_type_df$degree_dist$share) <= 1,
      failure_message = "Some degree shares are out of bounds"
    )
  }
})

test_that("GOF esp stats have the right fields and terms", {
  sim <- get_dummy_net(50, 2)
  g <- sim$g

  test_gof_res <- lighthergm::gof_lighthergm(
    g,
    list_feature_matrices = sim$list_feature_matrices,
    data_for_simulation = sim$nodes_data,
    colname_vertex_id = sim$vertex_id_var,
    colname_block_membership = sim$block_id_var,
    lighthergm_results = sim$hergm_res,
    ergm_control = sim$ergm_control,
    n_sim = 3
  )

  for (stat_type in c("original", "simulated")) {
    stat_type_df <- test_gof_res[[stat_type]]
    actual_terms <- colnames(stat_type_df$esp_dist)
    actual_terms[stringr::str_detect(actual_terms, "n_sim", negate = TRUE)] %>%
      setdiff(c("label", "esp")) %>%
      length() %>%
      expect_equal(0)

    expect_lte(max(stat_type_df$esp_dist$label), min(g$gal$n, 10))
    expect(
      min(stat_type_df$esp_dist$esp) >= 0 && max(stat_type_df$esp_dist$esp) <= (g$gal$n^2),
      failure_message = "Some esp counts are out of bounds."
    )
  }
})

test_that("GOF geodesic distance is returned when requested", {
  sim <- get_dummy_net(50, 2)
  g <- sim$g

  test_gof_res <- lighthergm::gof_lighthergm(
    g,
    list_feature_matrices = sim$list_feature_matrices,
    data_for_simulation = sim$nodes_data,
    colname_vertex_id = sim$vertex_id_var,
    colname_block_membership = sim$block_id_var,
    lighthergm_results = sim$hergm_res,
    ergm_control = sim$ergm_control,
    n_sim = 3,
    compute_geodesic_distance = TRUE
  )

  for (stat_type in c("original", "simulated")) {
    stat_type_df <- test_gof_res[[stat_type]]
    actual_terms <- colnames(stat_type_df$geodesic_dist)

    actual_terms[stringr::str_detect(actual_terms, "n_sim", negate = TRUE)] %>%
      setdiff(c("dist", "pairs")) %>%
      length() %>%
      expect_equal(0)

    # Some of the distances will be Inf, and that's ok (that's how ergm returns them).
    non_inf <- stat_type_df$geodesic_dist$dist[!is.infinite(stat_type_df$geodesic_dist$dist)]
    expect_lte(max(non_inf), g$gal$n)
    expect(
      (min(stat_type_df$geodesic_dist$pairs) >= 0) && (max(stat_type_df$geodesic_dist$pairs) <= (g$gal$n^2)),
      failure_message = "Some geodesic distance pairs are out of bounds."
    )
  }
})

test_that("Return GOF statistics including only within-block connections", {
  sim <- get_dummy_net(50, 2)
  g <- sim$g

  test_gof_res <- lighthergm::gof_lighthergm(
    g,
    list_feature_matrices = sim$list_feature_matrices,
    data_for_simulation = sim$nodes_data,
    colname_vertex_id = sim$vertex_id_var,
    colname_block_membership = sim$block_id_var,
    lighthergm_results = sim$hergm_res,
    type = 'within',
    ergm_control = sim$ergm_control,
    n_sim = 3
  )

  # check that the network stats belong to the within-block sub network only
  edgelist <- network::as.edgelist(g) %>% as.data.frame
  colnames(edgelist) <- c('src', 'dst')
  nodes_with_blocks <- data.frame(id = 1:length(network::network.vertex.names(g)), block=network::get.vertex.attribute(g, 'block'))
  actual_within_conns <- edgelist %>%
    dplyr::left_join(nodes_with_blocks, by = c('src' = 'id')) %>%
    dplyr::left_join(nodes_with_blocks, by = c('dst' = 'id'), suffix=c('.src', '.dst')) %>%
    dplyr::filter(block.src == block.dst) %>%
    nrow

  within_conns_from_gof <- (test_gof_res$original$network_stats %>% dplyr::filter(stat == 'edges'))[, 2]

  expect_equal(within_conns_from_gof, actual_within_conns)

  for (stat_type in c("original", "simulated")) {
    stats <- test_gof_res[[stat_type]]
    expect_false(is.null(stats))
    for (stat in c("network_stats", "degree_dist", "esp_dist")) {
      expect_false(is.null(stats[[stat]]))
    }
    expect_true(is.null(stats[["geodesic_dist"]]))
  }
})

test_that("Within-connections GOF can be started from the observed network", {
  sim <- get_dummy_net(100, 4)
  g <- sim$g

  ergm_control <- ergm::control.simulate.formula(
    MCMC.burnin = 0,
    MCMC.interval = 1
  )

  test_gof_res <- lighthergm::gof_lighthergm(
    g,
    list_feature_matrices = sim$list_feature_matrices,
    data_for_simulation = sim$nodes_data,
    colname_vertex_id = sim$vertex_id_var,
    colname_block_membership = sim$block_id_var,
    lighthergm_results = sim$hergm_res,
    type = 'within',
    ergm_control = ergm_control,
    n_sim = 2,
    start_from_observed = TRUE
  )

  first_simulation_stats <-test_gof_res$simulated$network_stats %>%
    dplyr::filter(n_sim == 1) %>%
    dplyr::select(-n_sim)

  original_network_stats <- test_gof_res$original$network_stats
  expect_equal(original_network_stats, first_simulation_stats)
})

test_that("Full GOF can be started from the observed network", {
  sim <- get_dummy_net(100, 4)
  g <- sim$g

  ergm_control <- ergm::control.simulate.formula(
    MCMC.burnin = 0,
    MCMC.interval = 1
  )

  test_gof_res <- lighthergm::gof_lighthergm(
    g,
    list_feature_matrices = sim$list_feature_matrices,
    data_for_simulation = sim$nodes_data,
    colname_vertex_id = sim$vertex_id_var,
    colname_block_membership = sim$block_id_var,
    lighthergm_results = sim$hergm_res,
    type = 'full',
    ergm_control = ergm_control,
    n_sim = 2,
    start_from_observed = TRUE
  )

  first_simulation_stats <-test_gof_res$simulated$network_stats %>%
    dplyr::filter(n_sim == 1)

  # If it starts from the observed network, the stats should not be zero
  expect_true(all(first_simulation_stats['value'] > 0))
})
