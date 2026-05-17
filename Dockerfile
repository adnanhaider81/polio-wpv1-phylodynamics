FROM rocker/r-ver:4.5.2

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       build-essential \
       git \
       libcurl4-openssl-dev \
       libfontconfig1-dev \
       libfreetype6-dev \
       libfribidi-dev \
       libharfbuzz-dev \
       libjpeg-dev \
       libpng-dev \
       libssl-dev \
       libtiff5-dev \
       libxml2-dev \
       make \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace
COPY . /workspace

RUN Rscript scripts/install_packages.R

CMD ["Rscript", "run_pipeline.R", "--root", ".", "--mode", "distance", "--nsim", "3", "--seed", "1", "--skip_figures"]
