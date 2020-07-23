FROM python:3.8-alpine as base

LABEL maintainer="mmiller@bromberglab.org" \
      description="avadx docker image (https://services.bromberglab.org/avadx)"

FROM base as builder

# setup system
WORKDIR /install

# setup app
COPY ./avadx /app/python/app
COPY ./python /app/python/avadx
COPY ./requirements.txt /tmp/requirements.txt
RUN pip install --upgrade pip && pip install --no-warn-script-location --prefix=/install -r /tmp/requirements.txt

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
ENTRYPOINT ["python", "-m"]

# set app CMD
CMD ["avadx.pipeline", "--help"]
