FROM python:3.8-alpine as base

LABEL maintainer="mmiller@bromberglab.org" \
      description="avadx docker image (https://services.bromberglab.org/avadx)"

FROM base as builder

# setup system
WORKDIR /install
RUN apk --no-cache add build-base gcc make gfortran linux-headers openblas-dev lapack-dev jpeg-dev zlib-dev wget git git-lfs

# setup app
RUN git clone --depth 1 https://bitbucket.org/bromberglab/avadx.git && \
 mkdir -p /app/python && mv avadx/python /app/python/avadx && rm -rf avadx && \
 pip install --upgrade pip && pip install --no-warn-script-location --prefix=/install https://bitbucket.org/bromberglab/avadx/get/master.zip scipy  && \
 git ls-remote https://bitbucket.org/bromberglab/avadx HEAD > $(python -c "import avadx as _; print(_.__path__[0])")/.build && \
 python -c 'import compileall; compileall.compile_dir("/app/python/avadx", force=1)'

FROM base

COPY --from=builder /install /usr/local
COPY --from=builder /app /app

RUN apk --no-cache add bash zip p7zip outils-md5 openblas-dev lapack-dev libjpeg-turbo

# setup bio-node
LABEL bio-node=v1.0

# set environment variables
WORKDIR /app
ENV PYTHONPATH=/app/python

# set app ENTRYPOINT
ENTRYPOINT ["avadx"]

# set app CMD
CMD ["--help"]
