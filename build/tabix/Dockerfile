FROM biocontainers/tabix:v1.9-11-deb_cv1 as base

USER root
RUN apt-get update && apt-get install -y \
 python3 \
 && apt-get clean && rm -rf /var/lib/apt/lists/*
USER biodocker

COPY ./python /app/python/avadx
COPY ./bash /app/bash/avadx

# setup bio-node
LABEL bio-node=v1.0

# set environment variables
WORKDIR /app
ENV PYTHONPATH=/app/python

# set app ENTRYPOINT
ENTRYPOINT ["python3", "-m"]

# set app CMD
CMD ["avadx.filterVCF_by_gnomAD.py"]
