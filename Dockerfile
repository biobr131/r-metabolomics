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
RUN R -e 'BiocManager::install(c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR", "fgsea", "devtools", "crmn", "hmdbQuery", "metabolomicsWorkbenchR", "qvalue"), dependencies=T)' && \
    R -e 'install.packages(c("devtools", "tidyverse"), dependencies=T)' && \
    R -e 'library(devtools)' && \
    R -e 'devtools::install_github("xia-lab/MetaboAnalystR", build=T, build_vignettes=T, build_manual=T)' && \
	R -e 'devtools::install_github("xia-lab/OptiLCMS", build=T, build_vignettes=F, build_manual=T)'
