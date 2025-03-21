install.packages(c(
  "tidyverse",
  "ellipse",
  "rjson",
  "loadings",
  "qvalue",
  "anndata",
  "Seurat"
), dependencies = TRUE)

BiocManager::install(c(
  "impute",
  "pcaMethods",
  "globaltest",
  "GlobalAncova",
  "Rgraphviz",
  "preprocessCore",
  "genefilter",
  "sva",
  "limma",
  "KEGGgraph",
  "siggenes",
  "BiocParallel",
  "MSnbase",
  "multtest",
  "RBGL",
  "edgeR",
  "fgsea",
  "devtools",
  "crmn",
  "httr",
  "qs"
), dependencies = TRUE)
