FROM r-base:3.6.3 as base

FROM base as builder

RUN apt-get update && apt-get install -y \
 git \
 && apt-get clean && rm -rf /var/lib/apt/lists/*
WORKDIR /install

# setup app
RUN git clone --depth 1 https://bitbucket.org/bromberglab/avadx.git && \
 mkdir -p /app/R && mv avadx/R /app/R/avadx && \
 rm -rf avadx

FROM base

COPY --from=builder /install /usr/local
COPY --from=builder /app /app

RUN apt-get update && apt-get install -y \
 file \
 curl \
 libcurl4-openssl-dev \
 openssl \
 libssl-dev \
 libxml2-dev \
 openjdk-11-jdk-headless \
 && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN ln -s /usr/lib/jvm/java-11-openjdk-amd64 /usr/lib/jvm/default-java

RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")'
RUN R -e 'BiocManager::install(c("SNPRelate"))'
RUN R -e "install.packages(c('optparse', 'data.table', 'tidyverse', 'seqinr', 'stringr', 'EthSEQ', 'e1071', 'caret', 'ggfortify', 'R.utils', 'PRROC', 'xlsx', 'ranger'), repos = 'http://cran.us.r-project.org', lib = '/usr/local/lib/R/site-library')"

# setup bio-node
LABEL bio-node=v1.0

# set environment variables
WORKDIR /app
ENV R_LIBS_USER=/app/R

# set app ENTRYPOINT
ENTRYPOINT ["Rscript"]

# set app CMD
CMD ["--version"]
