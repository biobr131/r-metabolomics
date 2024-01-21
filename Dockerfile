FROM bioconductor/bioconductor_docker:latest
RUN apt-get update && \
	apt-get install -y --no-install-recommends apt-utils && \
	apt-get install -y --no-install-recommends \
	texlive \
	texlive-latex-extra \
	texlive-fonts-extra \
	texlive-bibtex-extra \
	texlive-science \
	texi2html \
	texinfo && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/*
RUN R -e 'BiocManager::install(c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR", "fgsea", "devtools", "crmn", "hmdbQuery", "metabolomicsWorkbenchR", "qvalue"))' && \
    R -e 'install.packages(c("devtools", "tidyverse"))' && \
    R -e 'library(devtools)' && \
    R -e 'devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = TRUE, build_manual =T)'
