# The analysis image contains R and its pinned libraries. The project and Codex
# state are persistent host mounts, independent of the image's lifecycle.
FROM rocker/r-ver:4.3.2

ENV TZ=Etc/UTC \
    LC_ALL=C.UTF-8 \
    R_LIBS_USER=/opt/tumorcast-library \
    OPENBLAS_NUM_THREADS=1 \
    OMP_NUM_THREADS=1 \
    LIBARROW_BINARY=true

RUN apt-get update && apt-get install -y --no-install-recommends \
    git ripgrep curl ca-certificates cmake make g++ gfortran \
    libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
    libhdf5-dev libglpk-dev libpng-dev libjpeg-dev libtiff-dev \
    libfontconfig1-dev libfreetype6-dev libcairo2-dev librsvg2-dev \
    ghostscript poppler-utils fonts-urw-base35 \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /opt/tumorcast-library /opt/tumorcast
COPY env/renv.lock /opt/tumorcast/renv.lock
COPY env/install_packages.R /opt/tumorcast/install_packages.R
RUN Rscript --vanilla /opt/tumorcast/install_packages.R

WORKDIR /work/TUMORCAST
CMD ["sleep", "infinity"]
