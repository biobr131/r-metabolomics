source("R/utils_anndata.R")

adata <- create_anndata_from_excel(
  excel_path = "metabolomics_data.xlsx",
  data_sheet = "data",
  samples_sheet = "samples",
  metabolites_sheet = "metabolites",
  row_names_col = "Compound ID"
)

rsd_results <- calculate_qc_rsd(
  adata,
  sample_type_col = "sample_type",
  qc_value = "pooled_qc"
)

head(rsd_results)
print_rsd_summary(rsd_results)

plot_rsd_distribution(rsd_results)

poor_quality_metabolites <- subset(rsd_results, rsd > 30)
print(poor_quality_metabolites)

write.csv(rsd_results, "qc_rsd_results.csv", row.names = FALSE)
