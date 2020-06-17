FROM ubuntu:20.04 as builder

USER root

# ALL tool versions used by opt-build.sh
# ensure updated in setup.sh too
ENV VER_BEDTOOLS="2.28.0"
ENV VER_VCFTOOLS="0.1.16"
ENV VER_BIODBHTS="3.01"
ENV VER_HTSLIB="1.10.2"
ENV VER_SAMTOOLS="1.10"

RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
locales \
libperlio-gzip-perl \
curl \
ca-certificates \
zlib1g-dev \
libbz2-dev \
liblzma-dev \
libcurl4-gnutls-dev \
libncurses5-dev \
g++ \
make \
gcc \
pkg-config \
libgd-dev \
libdb-dev

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV OPT /opt/wtsi-cgp
ENV PATH $OPT/bin:$PATH
ENV PERL5LIB $OPT/lib/perl5
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

# don't work in the default location, it can cause problems
WORKDIR /tmp/builder

# build tools from other repos
ADD build/opt-build.sh build/
RUN bash build/opt-build.sh $OPT

# build the tools in this repo, separate to reduce build time on errors
COPY . .
RUN bash build/opt-build-local.sh $OPT

FROM ubuntu:20.04

LABEL maintainer="cgphelp@sanger.ac.uk" \
      uk.ac.sanger.cgp="Cancer, Ageing and Somatic Mutation, Wellcome Trust Sanger Institute" \
      version="3.6.1" \
      description="VAGrENT genome annotation docker"

RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
locales \
libperlio-gzip-perl \
curl \
ca-certificates \
bzip2 \
time \
zlib1g \
liblzma5 \
libncurses5 \
unattended-upgrades && \
unattended-upgrade -d -v && \
apt-get remove -yq unattended-upgrades && \
apt-get autoremove -yq

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV OPT /opt/wtsi-cgp
ENV PATH $OPT/bin:$PATH
ENV PERL5LIB $OPT/lib/perl5
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

RUN mkdir -p $OPT
COPY --from=builder $OPT $OPT

## USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER    ubuntu
WORKDIR /home/ubuntu

CMD ["/bin/bash"]
