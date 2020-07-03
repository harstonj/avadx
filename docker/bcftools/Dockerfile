FROM biocontainers/bcftools:v1.9-1-deb_cv1 as base

# setup bio-node
LABEL bio-node=v1.0

# set environment variables
WORKDIR /app

# set app ENTRYPOINT
ENTRYPOINT ["bcftools"]

# set app CMD
CMD ["--version"]
