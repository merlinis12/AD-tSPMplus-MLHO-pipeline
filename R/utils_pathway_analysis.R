## =============================================================================
## utils_pathway_analysis.R
## Build a directed graph of selected transitions, enumerate multi-step
## pathways, and test their prevalence / association with the outcome.
## =============================================================================

#' Split "A->B" transition labels into a from/to edge list
#'
#' @param transitions character vector of "A->B" labels
#' @return data.frame(from, to)
transitions_to_edgelist <- function(transitions) {
  tibble::tibble(
    transition = stringr::str_replace_all(transitions, "\\s+", ""),
    from = stringr::str_extract(transitions, ".*(?=->)"),
    to   = stringr::str_replace(transitions, "^.*->", "")
  )
}

#' Build a directed igraph object from an edge list and enumerate all simple
#' paths between every pair of nodes up to `max_len` edges.
#'
#' @param edges data.frame(from, to)
#' @param max_len maximum path length (number of nodes)
#' @return list(graph = igraph object,
#'              paths = tibble(path, nodes, length))
enumerate_pathways <- function(edges, max_len = 7) {
  g <- igraph::graph_from_data_frame(edges, directed = TRUE)

  node_names <- igraph::V(g)$name
  all_paths <- list()
  for (src in node_names) {
    for (tgt in node_names) {
      if (src == tgt) next
      found <- igraph::all_simple_paths(g, from = src, to = tgt, cutoff = max_len)
      for (p in found) {
        nodes <- node_names[p]
        all_paths[[paste(nodes, collapse = "->")]] <- nodes
      }
    }
  }

  paths_df <- tibble::tibble(path = names(all_paths), nodes = all_paths) %>%
    dplyr::mutate(length = purrr::map_int(nodes, length)) %>%
    dplyr::arrange(dplyr::desc(length))

  list(graph = g, paths = paths_df)
}

#' Check whether `path_nodes` occurs as a consecutive subsequence of a
#' patient's ordered event vector.
path_occurs <- function(events, path_nodes) {
  n <- length(path_nodes); m <- length(events)
  if (m < n) return(FALSE)
  for (i in seq_len(m - n + 1)) {
    if (all(events[i:(i + n - 1)] == path_nodes)) return(TRUE)
  }
  FALSE
}

#' Compute prevalence and outcome association (Fisher's exact test / OR) for
#' every enumerated pathway.
#'
#' @param paths_df output of enumerate_pathways()$paths
#' @param patients tibble with one row per patient, containing:
#'   `events` (list-column of ordered event labels) and `outcome_flag` (0/1)
compute_pathway_statistics <- function(paths_df, patients) {

  path_prevalence <- function(path_nodes) {
    occurs <- purrr::map_lgl(patients$events, path_occurs, path_nodes = path_nodes)
    tab <- table(occurs, patients$outcome_flag)

    a <- if (nrow(tab) >= 2) tab[2, 1] else 0  # occurs & outcome = 0
    b <- if (nrow(tab) >= 2) tab[2, 2] else 0  # occurs & outcome = 1
    c <- tab[1, 1]                              # !occurs & outcome = 0
    d <- tab[1, 2]                              # !occurs & outcome = 1

    ft <- tryCatch(
      fisher.test(matrix(c(a, b, c, d), nrow = 2)),
      error = function(e) list(p.value = NA, estimate = NA, conf.int = c(NA, NA))
    )

    tibble::tibble(
      n_patients = sum(occurs),
      prevalence = mean(occurs),
      occ_and_outcome0 = a, occ_and_outcome1 = b,
      noocc_and_outcome0 = c, noocc_and_outcome1 = d,
      or = suppressWarnings(as.numeric(ft$estimate)),
      pval = ft$p.value,
      conf_low  = suppressWarnings(ft$conf.int[1]),
      conf_high = suppressWarnings(ft$conf.int[2])
    )
  }

  paths_df %>%
    dplyr::mutate(stats = purrr::map(nodes, path_prevalence)) %>%
    tidyr::unnest(stats)
}
