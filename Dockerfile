FROM ubuntu:14.04
RUN apt-get update && apt-get install -y \
	autoconf \
	automake \
	make \
	g++ \
	gcc \
	build-essential \ 
	zlib1g-dev \
	libgsl0-dev \
	perl \
	curl \
	git \
	wget \
	unzip \
	tabix \
	libncurses5-dev

RUN apt-get install -y cpanminus
RUN apt-get install -y libmysqlclient-dev
RUN cpanm CPAN::Meta \
	Archive::Zip \
	DBI \
	DBD::mysql \ 
	JSON \
	DBD::SQLite \
	Set::IntervalTree \
	LWP \
	LWP::Simple \
	Archive::Extract \
	Archive::Tar \
	Archive::Zip \
	CGI \
	Time::HiRes \
	Encode \
	File::Copy::Recursive \
	Perl::OSType \
	Module::Metadata version \
	Bio::Root::Version \
	TAP::Harness \
	Module::Build

WORKDIR /opt
RUN wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2
RUN tar jxf samtools-1.3.tar.bz2
WORKDIR /opt/samtools-1.3
RUN make
RUN make install

WORKDIR /opt
RUN rm samtools-1.3.tar.bz2

RUN wget https://github.com/Ensembl/ensembl-tools/archive/release/85.zip
RUN mkdir variant_effect_predictor_85
RUN mkdir variant_effect_predictor_85/cache
RUN unzip 85.zip -d variant_effect_predictor_85
RUN rm 85.zip 
WORKDIR /opt/variant_effect_predictor_85/ensembl-tools-release-85/scripts/variant_effect_predictor/
RUN perl INSTALL.pl --AUTO ap --PLUGINS LoF --CACHEDIR /opt/variant_effect_predictor_85/cache
WORKDIR /opt/variant_effect_predictor_85/cache/Plugins
RUN wget https://raw.githubusercontent.com/konradjk/loftee/master/splice_module.pl

WORKDIR /opt
ADD . /opt/vcf2maf 

COPY Dockerfile /opt/
MAINTAINER Michele Mattioni, Seven Bridges, <michele.mattioni@sbgenomics.com>
