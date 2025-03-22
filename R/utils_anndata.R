library(readxl)
library(stats)
library(ggplot2)
library(anndata)
library(KEGGREST)

#' Create AnnData object from Excel workbook
#'
#' @param excel_path Character, path to Excel file
#' @param data_sheet Character, name of the data matrix sheet (default: "data")
#' @param samples_sheet Character, name of the sample information sheet (default: "samples")
#' @param metabolites_sheet Character, name of the metabolite information sheet (default: "metabolites")
#' @param row_names_col Character, column name for row names in data sheet (default: NULL)
#' @return AnnData object
#' @import readxl
#' @import anndata
#' @export
create_anndata_from_excel <- function(excel_path, data_sheet = "data", samples_sheet = "samples", metabolites_sheet = "metabolites", row_names_col = NULL) {
  if (!file.exists(excel_path)) {
    stop("Excel file not found: ", excel_path)
  }

  sheets <- readxl::excel_sheets(excel_path)
  required_sheets <- c(data_sheet, samples_sheet, metabolites_sheet)
  missing_sheets <- setdiff(required_sheets, sheets)
  if (length(missing_sheets) > 0) {
    stop("Missing sheets: ", paste(missing_sheets, collapse = ", "))
  }

  tryCatch({
    data_df <- readxl::read_excel(excel_path, sheet = data_sheet)
    if (!is.null(row_names_col)) {
      row_names <- data_df[[row_names_col]]
      data_df <- data_df[, !colnames(data_df) %in% row_names_col]
    }

    samples_df <- readxl::read_excel(excel_path, sheet = samples_sheet)

    metabolites_df <- readxl::read_excel(excel_path, sheet = metabolites_sheet)

    X <- as.matrix(data_df)
    var_names <- colnames(data_df)
    obs_names <- if (!is.null(row_names_col)) row_names else rownames(data_df)

    adata <- AnnData(
      X = X,
      obs = metabolites_df,
      var = samples_df,
      obs_names = obs_names,
      var_names = var_names
    )

    return(adata)

  }, error = function(e) {
    stop("Error creating AnnData object: ", e$message)
  })
}

#' Validate Excel data format
#'
#' @param data_df Data frame of the main data matrix
#' @param samples_df Data frame of sample information
#' @param metabolites_df Data frame of metabolite information
#' @return List of validation messages
#' @keywords internal
.validate_excel_data <- function(data_df, samples_df, metabolites_df) {
  messages <- list()

  if (any(is.na(data_df))) {
    messages$missing_data <- "Missing values found in data matrix"
  }

  if (!all(sapply(data_df, is.numeric))) {
    messages$non_numeric <- "Non-numeric values found in data matrix"
  }

  if (any(duplicated(samples_df[[1]]))) {
    messages$duplicate_samples <- "Duplicate sample IDs found"
  }
  if (any(duplicated(metabolites_df[[1]]))) {
    messages$duplicate_metabolites <- "Duplicate metabolite IDs found"
  }

  return(messages)
}

#' Process Excel data before creating AnnData
#'
#' @param data_df Data frame of the main data matrix
#' @param log_transform Logical, whether to log transform the data
#' @param scale Logical, whether to scale the data
#' @return Processed data frame
#' @keywords internal
.process_data_matrix <- function(data_df, log_transform = FALSE, scale = FALSE) {
  data_matrix <- as.matrix(data_df)

  if (log_transform) {
    data_matrix <- log2(data_matrix + 1)
  }

  if (scale) {
    data_matrix <- scale(data_matrix)
  }

  return(data_matrix)
}

#' Clean sample metadata
#'
#' @param samples_df Data frame of sample information
#' @return Cleaned data frame
#' @keywords internal
.clean_sample_metadata <- function(samples_df) {
  samples_df[] <- lapply(samples_df, function(x) {
    if (is.character(x)) trimws(x)
    else x
  })

  categorical_cols <- sapply(samples_df, function(x) {
    length(unique(x)) < length(x) / 2
  })
  samples_df[categorical_cols] <- lapply(
    samples_df[categorical_cols],
    factor
  )

  return(samples_df)
}

#' Extract samples by condition from AnnData object
#'
#' @param adata AnnData object
#' @param column Character, column name to filter on
#' @param value Character or vector, value(s) to filter for
#' @param exact Logical, whether to use exact matching (default: TRUE)
#' @return AnnData object
#' @export
extract_samples_by_condition <- function(adata, column, value, exact = TRUE) {
  if (!.validate_anndata(adata)) {
    stop("Invalid AnnData object")
  }

  if (!column %in% colnames(adata$obs)) {
    stop(sprintf("Column '%s' not found", column))
  }

  if (exact) {
    idx <- which(adata$obs[[column]] %in% value)
  } else {
    pattern <- paste(value, collapse = "|")
    idx <- grep(pattern, adata$obs[[column]])
  }

  if (length(idx) == 0) {
    warning("No matching samples found")
    return(NULL)
  }

  adata_subset <- adata[idx, ]

  return(adata_subset)
}

#' @keywords internal
.validate_anndata <- function(adata) {
  if (!inherits(adata, "AnnData")) {
    return(FALSE)
  }

  if (is.null(adata$obs) || is.null(adata$X)) {
    return(FALSE)
  }

  return(TRUE)
}

#' Calculate RSD (Relative Standard Deviation) for QC samples
#'
#' @param adata AnnData object
#' @param sample_type_col Character, name of the sample type column
#' @param qc_value Character, value indicating QC samples
#' @param rsd_threshold Numeric, threshold value for RSD (default: 30)
#' @return data.frame with metabolite information and RSD values
#' @import stats
#' @export
calculate_qc_rsd <- function(adata,
                             sample_type_col,
                             qc_value,
                             rsd_threshold = 30) {

  adata_qc <- extract_samples_by_condition(
    adata = adata,
    column = sample_type_col,
    value = qc_value
  )

  if (is.null(adata_qc)) {
    stop("No QC samples found")
  }

  calculate_rsd <- function(x) {
    (sd(x) / mean(x)) * 100
  }

  rsd_values <- apply(adata_qc$X, 1, calculate_rsd)

  result <- data.frame(
    metabolite = rownames(adata_qc$X),
    rsd = rsd_values
  )

  if (!is.null(adata_qc$obs)) {
    result <- cbind(result, adata_qc$obs)
  }

  result <- result[order(result$rsd), ]

  attr(result, "summary") <- list(
    n_qc_samples = ncol(adata_qc$X),
    median_rsd = median(rsd_values),
    mean_rsd = mean(rsd_values),
    below = sum(rsd_values <= rsd_threshold) / length(rsd_values) * 100
  )

  return(result)
}

#' Print summary of QC RSD results
#'
#' @param rsd_result Result from calculate_qc_rsd
#' @export
print_rsd_summary <- function(rsd_result) {
  summary <- attr(rsd_result, "summary")
  cat("QC RSD Summary:\n")
  cat("Number of QC samples:", summary$n_qc_samples, "\n")
  cat("Median RSD:", round(summary$median_rsd, 2), "%\n")
  cat("Mean RSD:", round(summary$mean_rsd, 2), "%\n")
  cat("Percentage of metabolites with RSD <= 30%:",
      round(summary$rsd_threshold_20, 1), "%\n")
}

#' Plot RSD distribution
#'
#' @param rsd_result Result from calculate_qc_rsd
#' @param threshold Numeric, RSD threshold for highlighting (default: 20)
#' @import ggplot2
#' @export
plot_rsd_distribution <- function(rsd_result, threshold = 30) {
  ggplot(rsd_result, aes(x = rsd)) +
    geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
    geom_vline(xintercept = threshold, color = "red", linetype = "dashed") +
    labs(
      title = "Distribution of RSD values in QC samples",
      x = "RSD (%)",
      y = "Count"
    ) +
    theme_minimal()
}

#' Filter AnnData object based on RSD threshold
#'
#' @param adata AnnData object
#' @param rsd_result Result from calculate_qc_rsd
#' @param threshold Numeric, RSD threshold for filtering
#' @return Filtered AnnData object
#' @export
filter_by_rsd <- function(adata, rsd_result, threshold = 30) {
  keep_metabolites <- rsd_result$metabolite[rsd_result$rsd <= threshold]
  adata_filtered <- adata[rownames(adata$X) %in% keep_metabolites, ]
  return(adata_filtered)
}

#' Calculate pathway coverage from metabolites in AnnData
#'
#' @param adata AnnData object
#' @param kegg_id_col Character, column name containing KEGG compound IDs in obs
#' @param organism Character, KEGG organism code (default: "hsa" for human)
#' @return data.frame with pathway coverage information
#' @import KEGGREST
#' @export
calculate_pathway_coverage <- function(adata,
                                       kegg_id_col,
                                       organism = "hsa") {
  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    stop("Please install KEGGREST package: BiocManager::install('KEGGREST')")
  }

  if (!.validate_anndata(adata)) {
    stop("Invalid AnnData object")
  }
  if (!kegg_id_col %in% colnames(adata$obs)) {
    stop(sprintf("Column '%s' not found in metabolite table", kegg_id_col))
  }

  measured_compounds <- unique(adata$obs[[kegg_id_col]])
  measured_compounds <- measured_compounds[!is.na(measured_compounds)]

  tryCatch({
    pathways <- KEGGREST::keggList("pathway", organism)

    results <- data.frame(
      pathway_id = names(pathways),
      pathway_name = unname(pathways),
      total_compounds = NA,
      measured_compounds = NA,
      coverage = NA,
      measured_ids = "",
      stringsAsFactors = FALSE
    )

    for (i in seq_len(nrow(results))) {
      pathway_id <- results$pathway_id[i]
      pathway_compounds <- try({
        pathway_info <- KEGGREST::keggGet(pathway_id)
        compound_ids <- pathway_info[[1]]$COMPOUND
        if (is.null(compound_ids)) return(NULL)
        names(compound_ids)
      }, silent = TRUE)

    if (!inherits(pathway_compounds, "try-error") && !is.null(pathway_compounds)) {
      measured_in_pathway <- measured_compounds[measured_compounds %in% pathway_compounds]

        results$total_compounds[i] <- length(pathway_compounds)
        results$measured_compounds[i] <- length(measured_in_pathway)
        results$coverage[i] <- length(measured_in_pathway) / length(pathway_compounds) * 100
        results$measured_ids[i] <- paste(measured_in_pathway, collapse = ";")
      }
    }

    results <- results[!is.na(results$coverage), ]
    results <- results[order(-results$coverage), ]

    attr(results, "summary") <- list(
      total_pathways = nrow(results),
      mean_coverage = mean(results$coverage, na.rm = TRUE),
      median_coverage = median(results$coverage, na.rm = TRUE),
      measured_compounds = length(measured_compounds)
    )

    return(results)

  }, error = function(e) {
    stop("Error accessing KEGG database: ", e$message)
  })
}

#' Print summary of pathway coverage results
#'
#' @param coverage_result Result from calculate_pathway_coverage
#' @export
print_pathway_summary <- function(coverage_result) {
  summary <- attr(coverage_result, "summary")
  cat("Pathway Coverage Summary:\n")
  cat("Total pathways analyzed:", summary$total_pathways, "\n")
  cat("Mean coverage:", round(summary$mean_coverage, 2), "%\n")
  cat("Median coverage:", round(summary$median_coverage, 2), "%\n")
  cat("Total measured compounds:", summary$measured_compounds, "\n")
}

#' Plot pathway coverage distribution
#'
#' @param coverage_result Result from calculate_pathway_coverage
#' @param top_n Integer, number of top pathways to show (default: 20)
#' @import ggplot2
#' @export
plot_pathway_coverage <- function(coverage_result, top_n = 20) {
  top_pathways <- head(coverage_result, top_n)

  top_pathways$short_name <- substr(top_pathways$pathway_name, 1, 50)

  ggplot(top_pathways,
         aes(x = reorder(short_name, coverage), y = coverage)) +
    geom_bar(stat = "identity", fill = "lightblue") +
    coord_flip() +
    labs(
      title = sprintf("Top %d Pathways by Coverage", top_n),
      x = "Pathway",
      y = "Coverage (%)"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8)
    )
}
