FROM bioconductor/bioconductor_docker:latest

RUN apt-get update && \
    apt-get install -y \
        libcairo2-dev \
        libnetcdf-dev \
        libxml2 \
        libxt-dev \
        libssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY install_packages.R /tmp/
RUN Rscript /tmp/install_packages.R && \
    rm /tmp/install_packages.R

RUN R -e 'devtools::install_github(c( \
        "xia-lab/OptiLCMS", \
        "xia-lab/MetaboAnalystR" \
    ), \
    build = TRUE, \
    build_vignettes = FALSE)'
