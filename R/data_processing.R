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

# ライブラリの読み込み
library(KEGGREST)
source("R/pathway_coverage.R")

# パスウェイカバレッジの計算
coverage_results <- calculate_pathway_coverage(
    adata,
    kegg_id_col = "KEGG_ID",
    organism = "hsa"
)

# 結果の表示
head(coverage_results)
print_pathway_summary(coverage_results)

# カバレッジ分布のプロット
plot_pathway_coverage(coverage_results, top_n = 20)

# 特定のカバレッジ以上のパスウェイを抽出
high_coverage_pathways <- subset(coverage_results, coverage >= 50)
print(high_coverage_pathways)

# 結果の保存
write.csv(coverage_results, "pathway_coverage_results.csv", row.names = FALSE)
