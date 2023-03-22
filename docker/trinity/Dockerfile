FROM ubuntu:18.04
MAINTAINER Phelelani.Mpangase@wits.ac.za

RUN apt-get update \
    && apt-get install -y build-essential \
               apt-utils \
               autoconf \
               automake \
               cmake \
               curl \
               debconf-utils \
               git \
               libbz2-dev \
               libcurl4-openssl-dev \
               libgsl-dev \
               liblzma-dev \
               libncurses5-dev \
               libncursesw5-dev \
               libssl-dev \
               libtbb-dev \
               libz-dev \
               rsync \
               software-properties-common \
               unzip \
               wget \
               zlib1g-dev \
    && apt-add-repository -y ppa:openjdk-r/ppa \
    && apt-get install -y openjdk-8-jre \
    && add-apt-repository -y ppa:deadsnakes/ppa \
    && apt-get install -y python3.6 \
               python3.6-dev \
               python3.6-distutils \
               python3.6-apt \
    && rm -rf /var/lib/apt/lists/*

RUN ln -sf /usr/bin/python3.6 /usr/bin/python \
    && ln -sf /usr/bin/python3.6 /usr/bin/python3

RUN cd /opt \
    && curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py \
    && python3.6 get-pip.py \
    && pip3 install numpy

RUN cd /opt \
    && wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.9/jellyfish-2.2.9.tar.gz \
    && tar -vxf jellyfish-2.2.9.tar.gz \
    && cd jellyfish-2.2.9 \
    && ./configure --prefix=/opt/jellyfish-2.2.9 \
    && make \
    && make install \
    && ldconfig \
    && rm /opt/jellyfish-2.2.9.tar.gz

RUN cd /opt \
    && git clone https://github.com/samtools/htslib.git \
    && cd /opt/htslib/ \
    && git submodule update --init --recursive \
    && autoheader \
    && autoconf \
    && ./configure \
    && make \
    && make install

RUN cd /opt \
    && git clone https://github.com/samtools/samtools.git \
    && cd /opt/samtools/ \
    && autoheader \
    && autoconf \
    && ./configure \
    && make \
    && make install

RUN cd /opt \
    && wget --no-check-certificate https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-source.zip \
    && unzip bowtie2-2.4.2-source.zip \
    && cd bowtie2-2.4.2 \
    && make \
    && rm /opt/bowtie2-2.4.2-source.zip

RUN cd /opt \
    && wget -c https://github.com/COMBINE-lab/salmon/releases/download/v1.4.0/salmon-1.4.0_linux_x86_64.tar.gz \
    && tar -vxf salmon-1.4.0_linux_x86_64.tar.gz \
    && rm /opt/salmon-1.4.0_linux_x86_64.tar.gz

RUN cd /opt \
    && wget -c https://github.com/trinityrnaseq/trinityrnaseq/releases/download/v2.12.0/trinityrnaseq-v2.12.0.FULL.tar.gz \
    && tar -vxf trinityrnaseq-v2.12.0.FULL.tar.gz \
    && cd trinityrnaseq-v2.12.0 \
    && make \
    && make plugins \
    && make install \
    && rm /opt/trinityrnaseq-v2.12.0.FULL.tar.gz

ENV PATH=/opt/jellyfish-2.2.9/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/jellyfish-2.2.9/lib
ENV PATH=/opt/bowtie2-2.4.2:$PATH
ENV PATH=/opt/salmon-latest_linux_x86_64/bin:$PATH
ENV PATH=/usr/local/bin/trinityrnaseq-v2.12.0:$PATH
ENV JAVA_HOME=/usr/lib/jvm/java-8-oracle
ENV PYTHONPATH=/usr/local/lib/python3.6/dist-packages
