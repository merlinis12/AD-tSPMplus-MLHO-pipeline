## =============================================================================
## 03_pathway_network_analysis.R
##
## Step 3 of the pipeline: build a directed graph from the transitions
## selected in step 2, enumerate multi-step candidate pathways, and test
## their prevalence and association with the outcome in the full cohort.
##
## Inputs: outputs of 01_/02_ scripts
## Outputs:
##   - network_basic.svg              (graph of selected transitions)
##   - enumerated_paths.csv
##   - pathway_stats.csv              (prevalence + odds ratio per pathway)
##   - top_paths_by_or.svg / top_paths_by_prevalence.svg
## =============================================================================

source("R/00_setup.R")
source("config.R")
source("R/utils_pathway_analysis.R")

selected_features <- read.csv(file.path(config$output_dir, "selected_features.csv"))
events   <- readRDS(file.path(config$output_dir, "event_table_numeric_id.rds"))
patients <- readRDS(file.path(config$output_dir, "patient_metadata_numeric_id.rds"))

## ---- build graph of selected transitions -----------------------------------
edges <- transitions_to_edgelist(selected_features$feature_name)

svg(file.path(config$output_dir, "network_basic.svg"))
g_preview <- igraph::graph_from_data_frame(edges, directed = TRUE)
plot(g_preview, vertex.size = 30, vertex.label.cex = 0.8, edge.arrow.size = 0.6,
    main = "Selected transitions network")
dev.off()

## ---- enumerate multi-step pathways -----------------------------------------
pathway_result <- enumerate_pathways(edges, max_len = config$max_path_length)
write.csv(
  pathway_result$paths %>% dplyr::mutate(nodes = purrr::map_chr(nodes, paste, collapse = ";")),
  file.path(config$output_dir, "enumerated_paths.csv"), row.names = FALSE
)

svg(file.path(config$output_dir, "network_pathways.svg"), width = 15, height = 7)
ggraph::ggraph(pathway_result$graph, layout = "sugiyama") +
  ggraph::geom_edge_link(arrow = grid::arrow(length = grid::unit(3, "mm")), alpha = 0.5) +
  ggraph::geom_node_label(ggplot2::aes(label = name)) +
  ggplot2::theme_minimal()
dev.off()

## ---- per-patient ordered event sequences (for pathway occurrence checks) ---
patient_events <- events %>%
  dplyr::arrange(.data[[config$id_col]], start_date) %>%
  dplyr::group_by(.data[[config$id_col]]) %>%
  dplyr::summarise(events = list(DESCRIPTION), .groups = "drop") %>%
  dplyr::left_join(
    patients %>% dplyr::select(dplyr::all_of(c(config$id_col, config$outcome_col))),
    by = config$id_col
  ) %>%
  dplyr::rename(outcome_flag = dplyr::all_of(config$outcome_col))

## ---- prevalence + odds-ratio statistics for each pathway -------------------
pathway_stats <- compute_pathway_statistics(pathway_result$paths, patient_events)
write.csv(pathway_stats %>% dplyr::select(-nodes),
         file.path(config$output_dir, "pathway_stats.csv"), row.names = FALSE)

top_by_or <- pathway_stats %>% dplyr::arrange(dplyr::desc(abs(or))) %>% dplyr::slice(1:15)
p1 <- ggplot2::ggplot(top_by_or, ggplot2::aes(x = reorder(path, or), y = or)) +
  ggplot2::geom_col() + ggplot2::coord_flip() +
  ggplot2::labs(x = "Pathway", y = "Odds ratio (occurrence -> outcome)",
               title = "Top pathways by odds ratio") +
  ggplot2::theme_minimal()
ggplot2::ggsave(file.path(config$output_dir, "top_paths_by_or.svg"), p1, width = 10, height = 6)

top_by_prev <- pathway_stats %>% dplyr::arrange(dplyr::desc(prevalence)) %>% dplyr::slice(1:15)
p2 <- ggplot2::ggplot(top_by_prev, ggplot2::aes(x = reorder(path, prevalence), y = prevalence)) +
  ggplot2::geom_col() + ggplot2::coord_flip() +
  ggplot2::labs(x = "Pathway", y = "Prevalence (fraction of patients)",
               title = "Top pathways by prevalence") +
  ggplot2::theme_minimal()
ggplot2::ggsave(file.path(config$output_dir, "top_paths_by_prevalence.svg"), p2, width = 10, height = 6)
