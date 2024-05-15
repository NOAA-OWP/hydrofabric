FROM rocker/verse:4.4.0

LABEL org.opencontainers.image.licenses="GPL-2.0-or-later" \
      org.opencontainers.image.source="https://github.com/rocker-org/rocker-versioned2" \
      org.opencontainers.image.vendor="Rocker Project" \
      org.opencontainers.image.authors="Carl Boettiger <cboettig@ropensci.org>"
      
ENV MY_ENV_VAR="hydrofabric"
ENV PROJ_VERSION=9.3.1
ENV GDAL_VERSION=3.8.4
ENV GEOS_VERSION=3.12.1
ENV S6_VERSION=v2.1.0.2
ENV RSTUDIO_VERSION=2024.04.1+748
ENV DEFAULT_USER=hydrofabric_user
ENV PANDOC_VERSION=default
ENV SHINY_SERVER_VERSION=latest
ENV QUARTO_VERSION=latest

COPY scripts/experimental/install_geospatial_unstable.sh /rocker_scripts/experimental/install_geospatial_unstable.sh

RUN /rocker_scripts/experimental/install_geospatial_unstable.sh
RUN /rocker_scripts/install_rstudio.sh
RUN /rocker_scripts/install_shiny_server.sh
RUN /rocker_scripts/install_pandoc.sh
RUN /rocker_scripts/install_quarto.sh

# Update and install dependencies
RUN apt-get update && apt-get install -y \
    r-cran-devtools \
    cargo \
    binutils \
    libproj-dev \
    gdal-bin \
    && rm -rf /var/lib/apt/lists/*


RUN add-apt-repository ppa:saiarcot895/chromium-beta
RUN apt-get install -y chromium-browser


RUN R -e "install.packages('gifski', repos='http://cran.rstudio.com/')"
RUN R -q -e 'install.packages("devtools")'
RUN R -q -e 'install.packages("remotes")'
RUN R -q -e 'remotes::install_github("NOAA-OWP/hydrofabric")'
RUN R -q -e 'remotes::install_github("mikejohnson51/ngen.hydrofab")'
RUN R -q -e 'install.packages("mikejohnson51/climateR")'
RUN R -q -e 'install.packages("NOAA-OWP/hydrofabric")'
RUN R -q -e 'devtools::install_github("dblodgett-usgs/nhdplusTools")'
RUN R -q -e 'remotes::install_github("rstudio/webshot2")'

EXPOSE 8787

CMD ["/init"]