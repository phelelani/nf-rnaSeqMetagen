FROM ubuntu:18.04
MAINTAINER Phelelani.Mpangase@wits.ac.za

RUN apt-get update && apt-get install -y build-essential \
        git \
        software-properties-common \
        unzip \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN cd /opt \
    && git clone https://github.com/alexdobin/STAR.git \
    && cd STAR/source \
    && make STAR

ENV PATH=/opt/STAR/source:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
