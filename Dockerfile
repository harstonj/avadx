FROM python:3.8-alpine as base

LABEL maintainer="mmiller@bromberglab.org" \
      description="avadx-meta docker image (https://services.bromberglab.org/avadx-meta)"

FROM base as builder

# setup system
RUN mkdir /install
WORKDIR /install
# RUN apk update && apk add git && rm -rf /var/cache/apk/*

# setup app
COPY ./build/app /app/python/app
COPY ./python /app/python/avadx
# RUN pip install --upgrade pip && pip install --prefix=/install -r /app/requirements.txt
# RUN (chmod a+x /app/setup.sh; /app/setup.sh)

FROM base

COPY --from=builder /install /usr/local
COPY --from=builder /app /app

RUN apk --no-cache add bash

# setup bio-node
LABEL bio-node=v1.0

# set environment variables
WORKDIR /app
ENV PYTHONPATH=/app/python

# set app ENTRYPOINT
ENTRYPOINT ["python", "-m"]

# set app CMD
CMD ["app.pipeline", "--help"]
