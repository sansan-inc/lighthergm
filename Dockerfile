# If you would like to utilize parallelized computation, the R version should be less than 4.0.0.
# There seems some problem with OpenMP for R >= 4.0.0 when you use macOS.
# https://github.com/RcppCore/RcppArmadillo/issues/290
FROM rocker/verse:3.6.3

# Install Python for installing Python's infomap library.
RUN apt-get update && apt-get install -y python3-dev python3-pip

# Install infomap by pip3.
RUN pip3 install -U infomap

# Make it possible to install latest packages
RUN echo "options(repos = c(REPO_NAME = 'https://packagemanager.rstudio.com/all/__linux__/centos7/latest'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site

# Install lighthergm.
COPY ./lighthergm /lighthergm
RUN Rscript -e "devtools::install('lighthergm')"
