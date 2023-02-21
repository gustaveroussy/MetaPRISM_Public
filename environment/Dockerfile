FROM registry.codeocean.com/codeocean/ubuntu:20.04.2

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        git \
        libcairo2-dev \
        wget \
    && rm -rf /var/lib/apt/lists/*

COPY postInstall /
RUN /postInstall
