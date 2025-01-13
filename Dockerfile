FROM bioconductor/bioconductor_docker:latest
RUN apt-get update && \
	apt-get install -y libcairo2-dev libnetcdf-dev libxml2 libxt-dev libssl-dev && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/*
RUN R -e 'install.packages(c("tidyverse", "ellipse", "rjson", "loadings"), dependencies = TRUE)' && \
	R -e 'install.packages(c("qvalue"), dependencies = TRUE)' && \
	R -e 'BiocManager::install(c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "sva", "limma", "KEGGgraph", "siggenes", "BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR", "fgsea", "devtools", "crmn", "httr", "qs"), dependencies = TRUE)' && \
    R -e 'library(devtools)' && \
	R -e 'devtools::install_github("xia-lab/OptiLCMS", build = TRUE, build_vignettes = FALSE)' &&\
    R -e 'devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = FALSE)'
