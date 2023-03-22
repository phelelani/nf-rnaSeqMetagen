FROM ubuntu:18.04
MAINTAINER Phelelani.Mpangase@wits.ac.za

RUN apt-get update && apt-get install -y build-essential \
        curl \
        git \
        software-properties-common \
    && add-apt-repository -y ppa:deadsnakes/ppa \
    && apt-get install -y python3.7 \
               python3.7-dev \
               python3.7-distutils \
    && rm -rf /var/lib/apt/lists/*

RUN ln -sf /usr/bin/python3.7 /usr/bin/python \
    && ln -sf /usr/bin/python3.7 /usr/bin/python3

RUN cd /opt \
    && curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py \
    && python3.7 get-pip.py

RUN cd /opt \
    && git clone https://github.com/ewels/MultiQC.git \
    && cd MultiQC && python3.7 setup.py install \
    && cd /opt && rm -rf MultiQC

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
