## =============================================================================
## 04_patient_trajectory_analysis.R
##
## Step 4 of the pipeline: reconstruct each patient's ordered trajectory
## through the selected event categories, summarize the most common
## trajectories leading to each outcome, visualize them (Sankey diagram,
## frequency bar charts, sequence-state plots), and stratify by demographic /
## genetic-risk subgroup.
##
## Inputs: outputs of 01_/03_ scripts, plus a mapping from fine-grained event
##   DESCRIPTION to a small number of broader domains (edit `domain_map`
##   below to match your own event categories).
## Outputs:
##   - patient_sequences.csv
##   - sankey_trajectories.html
##   - top_trajectories_by_outcome.svg
##   - trajectory_frequency_by_group.svg
##   - sequence_state_distribution.svg (TraMineR)
## =============================================================================

source("R/00_setup.R")
source("config.R")
source("R/utils_trajectory_analysis.R")
source("R/utils_plotting.R")

events   <- readRDS(file.path(config$output_dir, "event_table_numeric_id.rds"))
patients <- readRDS(file.path(config$output_dir, "patient_metadata_numeric_id.rds"))
selected_features <- read.csv(file.path(config$output_dir, "selected_features.csv"))

## Optional: restrict events to those involved in the selected transitions,
## and/or to a specific subgroup (e.g., only patients with the outcome).
selected_descriptions <- unique(unlist(strsplit(selected_features$feature_name, "->")))
events <- dplyr::filter(events, DESCRIPTION %in% selected_descriptions)

## ---- group fine-grained events into broader clinical domains --------------
## Edit this mapping to reflect your own event categories.
domain_map <- c(
  "Cardiovascular_event_1" = "Cardiovascular",
  "Cardiovascular_event_2" = "Cardiovascular",
  "Metabolic_event_1"      = "Metabolic",
  "Endocrine_event_1"      = "Endocrine",
  "Neuropsychiatric_event_1" = "Neuropsychiatric",
  "Immune_event_1"         = "Immune/Anti-inflammatory"
)
events$DESCRIPTION <- dplyr::recode(events$DESCRIPTION, !!!domain_map, .default = "Other")

## ---- build ordered per-patient trajectories --------------------------------
patient_index <- patients %>%
  dplyr::rename(index_date = dplyr::all_of(config$outcome_date_col))

df_ordered <- build_patient_trajectories(
  events, patient_index,
  id_col = config$id_col, event_date_col = "start_date",
  index_date_col = "index_date", outcome_col = config$outcome_col,
  outcome_labels = c("Outcome", "Outcome-free"),
  covariate_cols = c("sex", "genetic_risk")   # edit to match your covariates
)

patient_sequences <- summarize_patient_sequences(
  df_ordered, id_col = config$id_col, event_date_col = "start_date",
  covariate_cols = c("sex", "genetic_risk")
)
write.csv(patient_sequences, file.path(config$output_dir, "patient_sequences.csv"),
         row.names = FALSE)

## ---- Sankey diagram of stepwise transitions --------------------------------
sankey_tables <- build_sankey_tables(df_ordered, id_col = config$id_col)
sankey <- networkD3::sankeyNetwork(
  Links = sankey_tables$links, Nodes = sankey_tables$nodes,
  Source = "IDsource", Target = "IDtarget", Value = "value", NodeID = "name",
  fontSize = 12, nodeWidth = 20
)
networkD3::saveNetwork(sankey, file.path(config$output_dir, "sankey_trajectories.html"))

## ---- most frequent trajectories, overall and by outcome --------------------
sequence_counts <- patient_sequences %>% dplyr::count(sequence, outcome, sort = TRUE)
write.csv(sequence_counts, file.path(config$output_dir, "sequence_counts.csv"), row.names = FALSE)

top_trajectories <- sequence_counts %>% dplyr::slice_max(n, n = 15)
p1 <- ggplot2::ggplot(top_trajectories, ggplot2::aes(x = reorder(sequence, n), y = n, fill = outcome)) +
  ggplot2::geom_col() + ggplot2::coord_flip() +
  ggplot2::labs(x = "Event sequence", y = "Number of patients",
               title = "Most common trajectories, by outcome") +
  ggplot2::theme_minimal()
ggplot2::ggsave(file.path(config$output_dir, "top_trajectories_by_outcome.svg"),
               p1, width = 10, height = 6)

## ---- trajectory frequency stratified by demographic/genetic subgroup ------
if (all(c("sex", "genetic_risk") %in% names(patient_sequences))) {
  sequence_counts_by_group <- patient_sequences %>%
    dplyr::count(sex, genetic_risk, sequence, outcome, name = "patients")

  p2 <- sequence_counts_by_group %>%
    dplyr::filter(outcome == "Outcome") %>%
    dplyr::group_by(sex, genetic_risk) %>%
    dplyr::slice_max(patients, n = 10) %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(sequence, patients), y = patients, fill = sex)) +
    ggplot2::geom_col(width = 0.7) + ggplot2::coord_flip() +
    ggplot2::facet_wrap(~genetic_risk) +
    ggplot2::labs(title = "Top trajectories leading to outcome, by subgroup",
                 x = NULL, y = "Number of patients") +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(file.path(config$output_dir, "trajectory_frequency_by_group.svg"),
                 p2, width = 9, height = 6.5)
}

## ---- sequence-state visualization (TraMineR) -------------------------------
events_list <- strsplit(gsub(" -> ", "-", patient_sequences$sequence), split = "-")
max_len <- max(lengths(events_list))
seq_matrix <- do.call(rbind, lapply(events_list, function(x) c(x, rep(NA, max_len - length(x)))))
colnames(seq_matrix) <- paste0("T", seq_len(max_len))

seq_obj <- TraMineR::seqdef(seq_matrix, right = "DEL", left = "DEL", gaps = "NA")

svg(file.path(config$output_dir, "sequence_state_distribution.svg"), width = 10, height = 6)
if ("genetic_risk" %in% names(patient_sequences)) {
  TraMineR::seqdplot(seq_obj, group = patient_sequences$genetic_risk)
} else {
  TraMineR::seqdplot(seq_obj)
}
dev.off()
