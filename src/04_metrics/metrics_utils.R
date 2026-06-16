compute_auprc <- function(truth, score) {
  truth <- as.integer(truth)
  score <- as.numeric(score)

  if (length(truth) == 0L || all(is.na(truth)) || length(unique(truth[!is.na(truth)])) < 2L) {
    return(NA_real_)
  }

  valid_idx <- which(!is.na(truth) & !is.na(score))
  truth <- truth[valid_idx]
  score <- score[valid_idx]

  if (!any(truth == 1L)) {
    return(NA_real_)
  }

  ordering <- order(score, decreasing = TRUE)
  truth <- truth[ordering]

  cumulative_positives <- cumsum(truth == 1L)
  precision <- cumulative_positives / seq_along(truth)
  recall <- cumulative_positives / sum(truth == 1L)
  recall_previous <- c(0, recall[-length(recall)])

  sum((recall - recall_previous) * precision)
}

compute_kendall_tau <- function(truth_rank, predicted_rank) {
  stats::cor(truth_rank, predicted_rank, method = "kendall", use = "complete.obs")
}

summarize_runtime <- function(runtime_sec, peak_ram_gb) {
  data.frame(
    runtime_sec = runtime_sec,
    peak_ram_gb = peak_ram_gb,
    stringsAsFactors = FALSE
  )
}

summarize_prediction_table <- function(prediction_table) {
  if (!all(c("truth", "score") %in% names(prediction_table))) {
    stop(
      "`prediction_table` must contain columns: truth, score.",
      call. = FALSE
    )
  }

  data.frame(
    auprc = compute_auprc(prediction_table$truth, prediction_table$score),
    kendall_tau = if (all(c("truth_rank", "predicted_rank") %in% names(prediction_table))) {
      compute_kendall_tau(prediction_table$truth_rank, prediction_table$predicted_rank)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

build_metric_summary <- function(results_manifest) {
  if (!all(c("method", "angle_degrees", "result_file") %in% names(results_manifest))) {
    stop(
      "`results_manifest` must contain columns: method, angle_degrees, result_file.",
      call. = FALSE
    )
  }

  summary_list <- lapply(seq_len(nrow(results_manifest)), function(i) {
    result_file <- results_manifest$result_file[[i]]
    if (!file.exists(result_file)) {
      return(data.frame(
        method = results_manifest$method[[i]],
        angle_degrees = results_manifest$angle_degrees[[i]],
        auprc = NA_real_,
        kendall_tau = NA_real_,
        runtime_sec = if ("runtime_sec" %in% names(results_manifest)) results_manifest$runtime_sec[[i]] else NA_real_,
        peak_ram_gb = if ("peak_ram_gb" %in% names(results_manifest)) results_manifest$peak_ram_gb[[i]] else NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    prediction_table <- read.csv(result_file, stringsAsFactors = FALSE)
    metrics <- summarize_prediction_table(prediction_table)

    data.frame(
      method = results_manifest$method[[i]],
      angle_degrees = results_manifest$angle_degrees[[i]],
      auprc = metrics$auprc,
      kendall_tau = metrics$kendall_tau,
      runtime_sec = if ("runtime_sec" %in% names(results_manifest)) results_manifest$runtime_sec[[i]] else NA_real_,
      peak_ram_gb = if ("peak_ram_gb" %in% names(results_manifest)) results_manifest$peak_ram_gb[[i]] else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, summary_list)
}
