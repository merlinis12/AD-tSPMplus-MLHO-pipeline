## =============================================================================
## 02_iterative_feature_selection.R
##
## Step 2 of the pipeline: repeated, cross-validated feature selection over
## the candidate transitions from step 1, using MLHO. Running many
## iterations on random train/test splits and keeping features that are
## selected consistently gives a more stable feature set than a single run.
##
## Inputs: outputs of 01_temporal_sequence_mining.R
## Outputs:
##   - selected_features.csv       (features selected in > min_selection_freq
##                                   iterations, ranked by frequency/importance)
##   - feature_importance_plot.svg
## =============================================================================

source("R/00_setup.R")
source("config.R")
source("R/utils_sequence_mining.R")
source("R/utils_plotting.R")

dbmart_numeric <- readRDS(file.path(config$output_dir, "numeric_event_table.rds"))
patients       <- readRDS(file.path(config$output_dir, "patient_metadata_numeric_id.rds"))
events         <- readRDS(file.path(config$output_dir, "event_table_numeric_id.rds"))
candidate_sequences <- read.csv(file.path(config$output_dir, "candidate_sequences.csv"))

labels <- patients %>%
  dplyr::select(dplyr::all_of(c(config$id_col, config$outcome_col))) %>%
  dplyr::rename(label = dplyr::all_of(config$outcome_col))

## Lookup used to translate numeric transition codes back to readable labels
concept_lookup <- events[!duplicated(events$phenx), c("phenx", "DESCRIPTION")]

## Sequence table in the format MLHO expects: one row per candidate
## transition-occurrence, with `phenx` holding the transition code.
seq_table <- candidate_sequences %>% dplyr::rename(phenx = sequence)

## ---- train / test split -----------------------------------------------------
all_ids  <- as.character(unique(seq_table[[config$id_col]]))
test_ids <- sample(all_ids, round(config$test_fraction * length(all_ids)))

data.table::setDT(seq_table)
train_dat <- seq_table[!(get(config$id_col) %in% test_ids)]
test_dat  <- seq_table[get(config$id_col) %in% test_ids]

train_dat[, row := .I]; train_dat$value.var <- 1

## ---- repeated feature selection --------------------------------------------
all_features <- vector("list", config$n_iterations)

for (i in seq_len(config$n_iterations)) {
  message("MLHO iteration ", i, " / ", config$n_iterations)

  test_dat_sub <- subset(test_dat, phenx %in% colnames(train_dat))
  test_dat_sub[, row := .I]; test_dat_sub$value.var <- 1

  model <- run_mlho_iteration(
    train_dat = train_dat, test_dat = test_dat_sub, labels = labels,
    outcome_col = "label", top_n = config$top_n_features,
    classifier = config$classifier, cv_folds = config$cv_folds,
    note = paste0("iteration_", i)
  )

  model$features$feature_name <- decode_transitions(
    model$features$features, dbmart_numeric$phenxLookUp, concept_lookup
  )
  all_features[[i]] <- model$features
}

## ---- aggregate importance across iterations --------------------------------
features_dt <- data.table::rbindlist(all_features, idcol = "iteration")

feature_summary <- features_dt[, .(
  mean_importance = mean(Overall, na.rm = TRUE),
  code = unique(features),
  freq_selected = .N
), by = feature_name][order(-freq_selected, -mean_importance)]

selected_features <- feature_summary[freq_selected > config$min_selection_freq]

write.csv(selected_features,
         file.path(config$output_dir, "selected_features.csv"), row.names = FALSE)

p <- plot_feature_importance(selected_features, "feature_name", "mean_importance",
                             top_n = 15, title = "Top selected transitions")
ggplot2::ggsave(file.path(config$output_dir, "feature_importance_plot.svg"),
               plot = p, width = 6, height = 4)
