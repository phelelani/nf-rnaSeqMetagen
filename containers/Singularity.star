Bootstrap:shub
From:singularityhub/ubuntu

%labels
Maintainer Phelelani.Mpangase@wits.ac.za
    
%post	       
#### Updates and essentials
apt-get update
apt-get install -y build-essential
apt-get install -y wget git zlib1g-dev unzip

## Install STAR Aligner
## From Source: https://github.com/alexdobin/STAR
cd /opt \
    && wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz \
    && tar -vxf 2.5.3a.tar.gz \
    && make STAR -C STAR-2.5.3a/source \
    && rm /opt/2.5.3a.tar.gz

%environment
### Add paths to environment
export PATH=/opt/STAR-2.5.3a/source:$PATH
