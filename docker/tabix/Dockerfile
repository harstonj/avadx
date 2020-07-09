FROM biocontainers/tabix:v1.9-11-deb_cv1 as base

FROM base as builder

# setup system
USER root
RUN apt-get update && apt-get install -y \
 git \
 && apt-get clean && rm -rf /var/lib/apt/lists/*
WORKDIR /install

# setup app
RUN git clone --depth 1 https://bitbucket.org/bromberglab/avadx-meta.git && \
 mkdir -p /app/python && mv avadx-meta/python /app/python/avadx && \
 mkdir -p /app/bash && mv avadx-meta/bash /app/bash/avadx && \
 rm -rf avadx-meta

FROM base

COPY --from=builder /install /usr/local

USER root
RUN apt-get update && apt-get install -y \
 python3 \
 && apt-get clean && rm -rf /var/lib/apt/lists/*

USER biodocker

COPY --from=builder /app /app

# setup bio-node
LABEL bio-node=v1.0

# set environment variables
WORKDIR /app
ENV PYTHONPATH=/app/python

# set app ENTRYPOINT
ENTRYPOINT ["python3", "-m"]

# set app CMD
CMD ["avadx.filterVCF_by_gnomAD.py"]
