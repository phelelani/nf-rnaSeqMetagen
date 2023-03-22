FROM ubuntu:18.04
MAINTAINER Phelelani.Mpangase@wits.ac.za

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		       ed \
		       less \
		       locales \
		       vim-tiny \
		       wget \
		       ca-certificates \
		       fonts-texgyre \
	&& rm -rf /var/lib/apt/lists/*

RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8

RUN echo "deb http://http.debian.net/debian sid main" > /etc/apt/sources.list.d/debian-unstable.list \
    && echo 'APT::Default-Release "testing";' > /etc/apt/apt.conf.d/default \
    && echo 'APT::Install-Recommends "false";' > /etc/apt/apt.conf.d/90local-no-recommends

ENV R_BASE_VERSION=4.0.5

RUN echo "deb http://deb.debian.org/debian experimental main" > /etc/apt/sources.list.d/experimental.list \
    && echo "deb [trusted=yes] https://eddelbuettel.github.io/ppaR400 ./" > /etc/apt/sources.list.d/edd-r4.list

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
               gcc-9-base \
               libopenblas0-pthread \
		       littler \
               r-cran-littler \
		       r-base=${R_BASE_VERSION}-* \
		       r-base-dev=${R_BASE_VERSION}-* \
               r-base-core=${R_BASE_VERSION}-* \
		       r-recommended=${R_BASE_VERSION}-* \
	&& ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/installBioc.r /usr/local/bin/installBioc.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/installDeps.r /usr/local/bin/installDeps.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
	&& install.r docopt \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
	&& rm -rf /var/lib/apt/lists/*

CMD ["R"]

RUN apt-get update \
    && apt-get install -y build-essential \
               git \
               curl \
               apt-transport-https

RUN R -e 'install.packages("stringr", repos="http://cloud.r-project.org/", dependencies=TRUE)'

RUN cd /opt \
    && git clone https://github.com/VCG/upset.git
