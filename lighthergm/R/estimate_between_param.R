#' Estimate between-block parameters by logit
#' @importFrom magrittr %>%
#' @importFrom ergm ergm
#' @param formula formula for estimating between-block parameters
#' @param network network object
#' @param block a vector that represents which node belongs to which node
#' @export
estimate_between_param <- function(formula, network, block) {
  # Create a data frame that stores node-block correspondences.
  df_block <-
    tibble::tibble(
      intergraph_id = 1:length(block),
      block = block
    )

  # Create an edgelist for within-block connections
  within_block_link <-
    intergraph::asDF(network)$edges %>%
    dplyr::rename(
      tail = V1,
      head = V2
    ) %>%
    dplyr::left_join(., df_block, by = c("tail" = "intergraph_id")) %>%
    dplyr::select(tail, head, block) %>%
    dplyr::rename(block_tail = block) %>%
    dplyr::left_join(., df_block, by = c("head" = "intergraph_id")) %>%
    dplyr::select(tail, head, block_tail, block) %>%
    dplyr::rename(block_head = block) %>%
    dplyr::filter(block_tail == block_head)

  # Get edge id for within-block links
  within_block_link_eid <-
    unlist(network::get.dyads.eids(
      x = network,
      tails = within_block_link$tail,
      heads = within_block_link$head
    ))

  # Create a network for logit estimation
  g_logit <- network

  # Delete within-block links
  network::delete.edges(x = g_logit, eid = within_block_link_eid)

  # Create a formula that contains only dyad-independent terms. i.e. exclude externality terms like triangle.
  terms <- ergm::ergm_model(formula)$terms
  varnames <-
    statnet.common::list_rhs.formula(formula) %>%
    as.character()
  dep_terms <-
    terms %>% purrr::map(function(t) {
      dep <- t$dependence
      is_dep <- is.null(dep) || dep
    }) %>% unlist()
  between_rhs <- varnames[!dep_terms]
  between_formula <- as.formula(paste("g_logit ~ ", paste(between_rhs, collapse = " + ")))

  # Estimate logit
  between_logit <- ergm(
    formula = between_formula,
    estimate = "MPLE"
  )

  # Remove unnecessary network objects
  between_logit$newnetwork <- NULL

  # Return the output
  return(between_logit)
}
