FROM ubuntu:18.04
MAINTAINER Phelelani.Mpangase@wits.ac.za

## UPDATE AND INSTALL REQUEIRED LIBRARIES ON THE BASE IMAGE
RUN apt-get update \
    && apt-get install --no-install-recommends --no-install-suggests -y build-essential \
               wget \
               git \
               curl \
               libz-dev \
               zlib1g-dev \
               ca-certificates \
               lftp \
               rsync \
    && rm -rf /var/lib/apt/lists/*

RUN cd /opt \
    && wget http://www.cbcb.umd.edu/software/jellyfish/jellyfish-1.1.11.tar.gz \
    && tar -vxf jellyfish-1.1.11.tar.gz \
    && cd jellyfish-1.1.11 \
    && ./configure \
    && make \
    && make install \
    && rm /opt/jellyfish-1.1.11.tar.gz

RUN cd /opt \
    && git clone https://github.com/DerrickWood/kraken2.git \
    && cd kraken2 \
    && ./install_kraken2.sh bin

RUN cd /opt \
    && wget https://github.com/marbl/Krona/releases/download/v2.8/KronaTools-2.8.tar \
    && tar -vxf KronaTools-2.8.tar \
    && cd KronaTools-2.8 \
    && ./install.pl \
    && rm /opt/KronaTools-2.8.tar

RUN cd /opt \
    && wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-x64-linux.tar.gz \
    && tar -vxf ncbi-blast-2.11.0+-x64-linux.tar.gz \
    && rm /opt/ncbi-blast-2.11.0+-x64-linux.tar.gz

## COMMENT OUT LINES 46 - 48 in the file: rsync_from_ncbi.pl
    # if (! ($full_path =~ s#^ftp://${qm_server}${qm_server_path}/##)) {
    #     die "$PROG: unexpected FTP path whatttt the duc (new server?) for $ftp_path\n";
    # }
## ADD NEW LINE 41 (ABOVE: next if $ftp_path eq "na";  # Skip if no provided path) in the file: rsync_from_ncbi.pl
    # $ftp_path = substr($ftp_path,37);
## CHANGE WGET TO LFTP IN LINE 33in the file: download_genomic_library.sh
    #  lftp -e 'pget -n10 ${FTP_SERVER}${file} . ; exit'
COPY download_genomic_library.sh /opt/kraken2/scripts/
COPY rsync_from_ncbi.pl /opt/kraken2/scripts/
COPY download_genomic_library.sh /opt/kraken2/bin/
COPY rsync_from_ncbi.pl /opt/kraken2/bin/

ENV PATH=/opt/kraken2/bin:$PATH
ENV LD_LIBRARY_PATH=/usr/local/lib
ENV PATH=/opt/ncbi-blast-2.11.0+/bin:$PATH
