# hash:sha256:a28f27f165a10b984e0ef86a1fe6bad14fdeb53867a037e9c9c65d144fabff8f
FROM registry.codeocean.com/codeocean/ubuntu:20.04.2

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        git=1:2.25.1-1ubuntu3.6 \
        libcairo2-dev=1.16.0-4ubuntu1 \
        wget=1.20.3-1ubuntu2 \
    && rm -rf /var/lib/apt/lists/*

COPY postInstall /
RUN /postInstall
