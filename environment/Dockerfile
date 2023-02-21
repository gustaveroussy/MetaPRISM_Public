# hash:sha256:4a385a7611d077163db456211aa3b2ffd00e4e38eb98ca082a39651e1b3cf959
FROM registry.codeocean.com/codeocean/ubuntu:20.04.2

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        git=1:2.25.1-1ubuntu3.10 \
        libcairo2-dev=1.16.0-4ubuntu1 \
        wget=1.20.3-1ubuntu2 \
    && rm -rf /var/lib/apt/lists/*

COPY postInstall /
RUN /postInstall

