FROM python:3.8-alpine as base

LABEL maintainer="mmiller@bromberglab.org" \
      description="avadx docker image (https://services.bromberglab.org/avadx)"

FROM base as builder

# setup system
WORKDIR /install

# setup app
RUN apk --no-cache add build-base gcc make wget git
RUN git clone --depth 1 https://bitbucket.org/bromberglab/avadx.git && \
 mkdir -p /app/python && mv avadx/python /app/python/avadx && \
 rm -rf avadx

RUN pip install --upgrade pip && pip install --no-warn-script-location --prefix=/install https://bitbucket.org/bromberglab/avadx/get/master.zip

FROM base

COPY --from=builder /install /usr/local
COPY --from=builder /app /app

RUN apk --no-cache add bash zip p7zip outils-md5

# setup bio-node
LABEL bio-node=v1.0

# set environment variables
WORKDIR /app
ENV PYTHONPATH=/app/python

# set app ENTRYPOINT
ENTRYPOINT ["avadx"]

# set app CMD
CMD ["--help"]
