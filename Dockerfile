FROM ubuntu:bionic

RUN apt update && apt install -y build-essential git cpanminus curl wget unzip automake samtools tabix
RUN apt install -y libmysqlclient-dev libncurses5-dev zlib1g-dev libgsl0-dev libexpat1-dev libgd-dev

RUN cpanm --notest LWP::Simple DBI DBD::mysql Archive::Zip Archive::Extract HTTP::Tiny Test::Simple
RUN cpanm --notest File::Copy::Recursive Perl::OSType Module::Metadata version TAP::Harness CGI Encode
RUN cpanm --notest CPAN::Meta JSON DBD::SQLite Set::IntervalTree Archive::Tar Time::HiRes Module::Build
RUN cpanm --notest Bio::Root::Version

WORKDIR /opt
RUN mkdir -p variant_effect_predictor_95/cache
RUN wget -q https://github.com/Ensembl/ensembl-vep/archive/release/95.3.tar.gz
RUN tar -zxf 95.3.tar.gz -C variant_effect_predictor_95 && rm 95.3.tar.gz

WORKDIR /opt/variant_effect_predictor_95
RUN perl ensembl-vep-release-95.3/INSTALL.pl --NO_TEST --NO_UPDATE --AUTO ap --PLUGINS LoF --DESTDIR ensembl-vep-release-95.3 --CACHEDIR cache

WORKDIR /opt/variant_effect_predictor_95/cache/Plugins
RUN wget -q https://raw.githubusercontent.com/konradjk/loftee/v0.3-beta/splice_module.pl

WORKDIR /opt
ADD . /opt/vcf2maf 
COPY Dockerfile /opt/

MAINTAINER Michele Mattioni, Seven Bridges, <michele.mattioni@sbgenomics.com>
