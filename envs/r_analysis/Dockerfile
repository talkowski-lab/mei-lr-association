## Environment for rendering ad-hoc Rmarkdown reports.
## rocker/verse ships R, tidyverse, pandoc, and TinyTeX (for PDF output) already.
FROM rocker/verse:4.4.1
RUN R -e "install.packages(c('optparse', 'jsonlite', 'AnVIL', 'AnVILGCP', 'data.table', 'ggplot2', 'ggthemes', 'glue'), repos = 'https://cloud.r-project.org')"

## Build samtools/bcftools from source with --enable-libcurl so htslib's remote
## file support can read gs:// (and s3://) URLs directly, not just local paths.
ARG HTS_VERSION=1.21
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev libncurses-dev \
    && rm -rf /var/lib/apt/lists/* \
    && cd /tmp \
    && for pkg in htslib samtools bcftools; do \
        curl -sSL "https://github.com/samtools/${pkg}/releases/download/${HTS_VERSION}/${pkg}-${HTS_VERSION}.tar.bz2" | tar xj; \
    done \
    && cd /tmp/htslib-${HTS_VERSION} && ./configure --enable-libcurl && make -j"$(nproc)" && make install \
    && cd /tmp/samtools-${HTS_VERSION} && ./configure --with-htslib=system && make -j"$(nproc)" && make install \
    && cd /tmp/bcftools-${HTS_VERSION} && ./configure --with-htslib=system && make -j"$(nproc)" && make install \
    && ldconfig \
    && rm -rf /tmp/htslib-${HTS_VERSION} /tmp/samtools-${HTS_VERSION} /tmp/bcftools-${HTS_VERSION}

COPY scripts/ ./scripts/
CMD ["bash"]
