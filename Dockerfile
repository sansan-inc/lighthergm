# If you would like to utilize parallelized computation, the R version should be less than 4.0.0.
# There seems some problem with OpenMP for R >= 4.0.0 when you use macOS.
# https://github.com/RcppCore/RcppArmadillo/issues/290
FROM rocker/verse:3.6.3

# Install Python for installing Python's infomap library.
RUN apt-get update && apt-get install -y python3-dev python3-pip

# Install infomap by pip3.
RUN pip3 install -U infomap

# Install the latest version of ergm.
# When you include the version tag, it seems you cannot install the latest version (at least >= 3.11.0) of ergm by default, which cause a dependency issue with lighthergm.
# To avoid that problem, ergm is downloaded from http://cran.rstudio.com/.
RUN Rscript -e "install.packages(c('RcppArmadillo', 'ergm'), repos = 'http://cran.rstudio.com/')"

# Install lighthergm.
COPY ./lighthergm /lighthergm
RUN Rscript -e "devtools::install('lighthergm')"
