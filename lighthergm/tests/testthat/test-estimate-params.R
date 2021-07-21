test_that("estimating between-block parameters by logit works", {
  set.seed(334)
  # Prepare data
  edgelist <-
    tibble::tribble(
      ~head, ~tail,
      1, 9,
      2, 6,
      2, 7,
      2, 9,
      3, 5,
      3, 9,
      4, 7,
      4, 11,
      4, 15,
      5, 11,
      5, 15,
      7, 8,
      7, 16,
      9, 13,
      9, 14,
      9, 16,
      10, 14,
      11, 15,
      13, 15,
      13, 16
    )
  edgelist <-
    as.matrix(edgelist)
  attr(edgelist, "n") <- 16
  attr(edgelist, "vnames") <-
    c(
      "Acciaiuoli", "Albizzi", "Barbadori", "Bischeri", "Castellani", "Ginori",
      "Guadagni", "Lamberteschi", "Medici", "Pazzi", "Peruzzi", "Pucci", "Ridolfi",
      "Salviati", "Strozzi", "Tornabuoni"
    )
  attr(edgelist, "directed") <- FALSE
  attr(edgelist, "bipartite") <- FALSE
  attr(edgelist, "loops") <- FALSE
  attr(edgelist, "class") <- c("edgelist", "matrix")

  g <- network::network(edgelist, matrix.type = "edgelist", directed = FALSE)

  x1 <- as.integer(unlist(purrr::rbernoulli(n = g$gal$n)))
  network::set.vertex.attribute(x = g, attrname = "x1", value = x1)

  # Cluster
  z_memb <- rep(1:4, each = 4)
  network::set.vertex.attribute(x = g, attrname = "block", value = z_memb)

  # Create dataset for test
  g_link <- intergraph::asDF(g)$edges
  g_attr <- intergraph::asDF(g)$vertexes

  df_g <-
    tibble::tibble(
      head = 1:g$gal$n,
      tail = 1:g$gal$n
    ) %>%
    tidyr::expand(tail, head) %>%
    dplyr::filter(tail < head) %>%
    dplyr::left_join(., g_attr, by = c("tail" = "intergraph_id")) %>%
    dplyr::left_join(., g_attr, by = c("head" = "intergraph_id")) %>%
    dplyr::mutate(
      nodematch.x1 = ifelse(x1.x == x1.y, 1, 0),
      same_block = ifelse(block.x == block.y, 1, 0)
    ) %>%
    dplyr::select(tail, head, nodematch.x1:same_block) %>%
    dplyr::left_join(., g_link, by = c("tail" = "V1", "head" = "V2")) %>%
    dplyr::mutate(connected = ifelse(is.na(na), 0, 1)) %>%
    dplyr::select(-na)

  # Estimate the model
  formula <- g ~ edges + nodematch("x1") + triangle + kstar(2)
  est_between <- estimate_between_param(
    formula = formula,
    network = g,
    block = z_memb
  )

  # Check if between-block connections are all zero.
  g_logit <- est_between$network
  edgelist <- intergraph::asDF(g_logit)$edges

  true_edgelist <-
    df_g %>%
    dplyr::filter(same_block == 0 & connected == 1) %>%
    dplyr::select(tail, head) %>%
    dplyr::arrange(tail, head)

  # Does it work!!!?
  expect_equal(edgelist$V1, true_edgelist$tail)
  expect_equal(edgelist$V2, true_edgelist$head)

  # Check if estimates for between-block parameters are the same.
  param_est <- stats::coef(est_between)
  logit_true <- glm(
    formula = connected ~ nodematch.x1,
    data = df_g %>% dplyr::mutate(connected = ifelse(same_block == 1, 0, connected)),
    family = "binomial"
  )

  param_est_true <- stats::coef(logit_true)

  # Does it work?
  expect_equal(param_est, param_est_true, check.attributes = FALSE, tolerance = 1e-7)

  # Check if within-block parameter estiamtion works
  expect_error(estimate_within_params(
    formula = formula,
    network = g,
    z_memb = z_memb,
    parallel = FALSE,
    verbose = 0,
    initial_estimate = NULL,
    seeds = NULL,
    method_second_step = "MPLE"
  ), NA)
})
