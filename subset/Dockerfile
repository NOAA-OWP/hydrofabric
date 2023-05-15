FROM rocker/geospatial:4.3.0

COPY --from=public.ecr.aws/lambda/provided:al2.2023.05.13.00 /lambda-entrypoint.sh /lambda-entrypoint.sh
COPY --from=public.ecr.aws/lambda/provided:al2.2023.05.13.00 /usr/local/bin/aws-lambda-rie /usr/local/bin/aws-lambda-rie
ENV LAMBDA_TASK_ROOT=/var/task
ENV LAMBDA_RUNTIME_DIR=/var/runtime
RUN mkdir /var/runtime
RUN mkdir /var/task
WORKDIR /var/task
ENTRYPOINT ["/lambda-entrypoint.sh"]

# Tnstall Apache Arrow/build tools
ENV HF_BUILD_PKGS="build-essential git"
ENV HF_ARROW_PKGS="libarrow-dev libarrow-glib-dev libarrow-dataset-dev libarrow-dataset-glib-dev libarrow-flight-dev libarrow-flight-glib-dev libparquet-dev libparquet-glib-dev"
RUN apt update \
    && apt install -y -V ca-certificates lsb-release wget \
    && wget https://apache.jfrog.io/artifactory/arrow/$(lsb_release --id --short | tr 'A-Z' 'a-z')/apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb \
    && apt install -y -V ./apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb \
    && apt update \
    && apt install -y -V ${HF_BUILD_PKGS} ${HF_ARROW_PKGS}


# Setup hydrofabric location
RUN mkdir -p /hydrofabric

# Install CRAN Packages
ENV HF_CRAN_PKGS="cli arrow crayon dplyr DBI RSQLite sf terra lambdr glue rstudioapi purrr magrittr"
RUN cd /hydrofabric \
    && . /etc/lsb-release \
    && echo "options(Ncpus = $(nproc --all), repos=c(CRAN = 'https://packagemanager.rstudio.com/cran/__linux__/${DISTRIB_CODENAME}/latest'))" >> .Rprofile \
    && install2.r devtools \
    && install2.r -n 6 -s ${HF_CRAN_PKGS}

# Install GH Packages
ENV HF_GH_PKGS="mikejohnson51/hydrofab mikejohnson51/ngen.hydrofab mikejohnson51/zonal mikejohnson51/climateR DOI-USGS/nhdplusTools"
RUN cd /hydrofabric \
    && installGithub.r ${HF_GH_PKGS}

COPY . /hydrofabric

# Setup Lambdr
RUN cd /hydrofabric \
    && Rscript -e "devtools::install_local()" \
    && chmod 755 subset/runtime.R \
    && printf "#!/bin/sh\ncd /hydrofabric/subset\nRscript runtime.R" > /var/runtime/bootstrap \
    && chmod +x /var/runtime/bootstrap

CMD ["subset"]