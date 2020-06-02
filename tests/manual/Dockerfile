# install latest base from Rocker
FROM rocker/r-base:latest

# Install Ubuntu packages
RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev/unstable \
    libxt-dev \
    libssl-dev \
	libxml2-dev
	

# install Imports
RUN R -e "install.packages(c('abind', 'dplyr','methods','grDevices', 'ggplot2', 'parallel', 'Rcpp', 'stats', 'utils'))"

# install Suggests
RUN R -e "install.packages(c('boot', 'broom', 'crayon', 'devtools', 'DT', 'fmsb','ggrepel','gridExtra','kableExtra', 'knitr','MASS', 'mvtnorm', 'openxlsx','purrr', 'readxl','rmarkdown', 'rfishbase', 'shiny', 'testthat', 'tidyr'))"


CMD build DLMtool

CMD check --as-cran DLMtool_5.4.4
			

# docker build --tag dlmtool_test .







