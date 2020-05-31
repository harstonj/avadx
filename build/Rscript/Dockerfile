FROM r-base:3.6.3 as base

RUN apt-get update && apt-get install -y \
 curl \
 libcurl4-openssl-dev \
 openssl \
 libssl-dev \
 libxml2-dev \
 && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")'
RUN R -e 'BiocManager::install(c("SNPRelate"))'
RUN R -e "install.packages(c('optparse', 'data.table', 'tidyverse', 'seqinr', 'stringr', 'EthSEQ', 'e1071', 'caret', 'ggfortify'), repos = 'http://cran.us.r-project.org', lib = '/usr/local/lib/R/site-library')"

COPY ./R /app/R/avadx

# setup bio-node
LABEL bio-node=v1.0

# set environment variables
WORKDIR /app
ENV R_LIBS_USER=/app/R

# set app ENTRYPOINT
ENTRYPOINT ["Rscript"]

# set app CMD
CMD ["--version"]
