vcf2maf
=======

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.14107.svg)](http://dx.doi.org/10.5281/zenodo.14107)

To convert a [VCF](http://samtools.github.io/hts-specs/) into a [MAF](https://wiki.nci.nih.gov/x/eJaPAQ), each variant must be mapped to only one of all possible gene transcripts/isoforms that it might affect. This selection of a single effect per variant, is often subjective. So this project is an attempt to make the selection criteria smarter, reproducible, and more configurable. And the default criteria must lean towards best practices. Per the current defaults, a single affected transcript is selected per variant, as follows:
 1. Sort effects first by transcript [biotype priority](https://github.com/ckandoth/vcf2maf/blob/master/vcf2maf.pl#L81), then by [effect severity](https://github.com/ckandoth/vcf2maf/blob/master/vcf2maf.pl#L24), and finally by decreasing transcript length
 2. Pick the gene on the top of the list (worst-affected), and choose it's [canonical](http://www.ensembl.org/Help/Glossary?id=346) transcript (VEP picks the longest [CCDS](http://www.ncbi.nlm.nih.gov/CCDS/) isoform)
 3. If the gene has no canonical transcript tagged (if you used snpEff), choose its longest transcript instead

Quick start
-----------

Download the latest release of vcf2maf, and view the detailed usage manual:

    curl -LO https://github.com/ckandoth/vcf2maf/archive/master.zip; unzip master.zip; cd vcf2maf-master
    perl vcf2maf.pl --man

To download properly versioned releases, [click here](https://github.com/ckandoth/vcf2maf/releases) for a list.

If you don't have [VEP](http://useast.ensembl.org/info/docs/tools/vep/index.html) or [snpEff](http://snpeff.sourceforge.net/) installed, see the sections below. VEP is preferred for it's CLIA-compliant [HGVS formats](http://www.hgvs.org/mutnomen/recs.html), and is used by default. So after installing VEP, you can test the script like so:

    perl vcf2maf.pl --input-vcf data/test.vcf --output-maf data/test.maf

If you'd rather use snpEff, there's an option for that:

    perl vcf2maf.pl --input-vcf data/test.vcf --output-maf data/test.snpeff.maf --use-snpeff

If you already have a VCF annotated with either VEP or snpEff, you can use those directly. You should have ran VEP with at least these options: `--everything --check_existing --total_length --allele_number --xref_refseq`. And for snpEff use these options: `-hgvs -sequenceOntology`. In older versions of snpEff, `-sequenceOntology` was incorrectly spelled `-sequenceOntolgy`. Feed your VEP/snpEff annotated VCFs into vcf2maf as follows:

    perl vcf2maf.pl --input-vep data/test.vep.vcf --output-maf data/test.maf
    perl vcf2maf.pl --input-snpeff data/test.snpeff.vcf --output-maf data/test.maf

To fill columns 16 and 17 of the output MAF with tumor/normal sample IDs, and to parse out genotypes and allele counts from matched genotype columns in the VCF, use options `--tumor-id` and `--normal-id`. Skip option `--normal-id` if you didn't have a matched normal:

    perl vcf2maf.pl --input-vcf data/test.vcf --output-maf data/test.maf --tumor-id WD1309 --normal-id NB1308

VCFs from variant callers like [VarScan](http://varscan.sourceforge.net/somatic-calling.html#somatic-output) use hardcoded sample IDs TUMOR/NORMAL in the genotype columns of the VCF. To have this script correctly parse the correct genotype columns, while still printing the proper IDs in the output MAF:

    perl vcf2maf.pl --input-vcf data/test_varscan.vcf --output-maf data/test_varscan.maf --tumor-id WD1309 --normal-id NB1308 --vcf-tumor-id TUMOR --vcf-normal-id NORMAL

If you have VEP in a different folder like `/opt/vep`, and cached in `/srv/vep`, there are options available to point the script there. Similar options available for snpEff too:

    perl vcf2maf.pl --input-vcf data/test.vcf --output-maf data/test.maf --vep-path /opt/vep --vep-data /srv/vep
    perl vcf2maf.pl --input-vcf data/test.vcf --output-maf data/test.maf --snpeff-path /opt/snpEff --snpeff-data /opt/snpEff/data --use-snpeff

Install VEP
-----------

Ensembl's VEP ([Variant Effect Predictor](http://useast.ensembl.org/info/docs/tools/vep/index.html)) is popular for how it selects a single "canonical transcript" per gene as [detailed here](http://useast.ensembl.org/Help/Glossary?id=346), its CLIA-compliant [HGVS variant format](http://www.hgvs.org/mutnomen/recs.html), and [Sequence Ontology nomenclature](http://useast.ensembl.org/info/genome/variation/predicted_data.html#consequences) for variant effects. It's download-able as a Perl script, so make sure you have [Perl installed](http://www.perl.org/get.html).

On Ubuntu/Debian, sudoers can install VEP's pre-requisite packages like this:

    sudo apt-get install -y curl rsync samtools tabix libarchive-extract-perl libarchive-zip-perl libwww-perl libcgi-pm-perl libdbi-perl libdbd-mysql-perl

On CentOS/Redhat/Fedora, here is how to install the equivalent packages:

    sudo yum -y install curl rsync samtools tabix perl-Archive-Extract perl-Archive-Zip perl-libwww-perl perl-CGI perl-DBI perl-DBD-mysql perl-Time-HiRes

Create temporary shell variables pointing to where we'll store VEP and its cache data (non default paths can be used, but specify `--vep-path` and `--vep-data` when running vcf2maf):

    export VEP_PATH=~/vep
    export VEP_DATA=~/.vep

Download the v78 release of VEP:

    mkdir $VEP_PATH; cd $VEP_PATH
    curl -LO https://github.com/Ensembl/ensembl-tools/archive/release/78.tar.gz
    tar -zxf 78.tar.gz --starting-file variant_effect_predictor --transform='s|.*/|./|g'

Download and unpack VEP's offline cache for GRCh37 and GRCh38:

    rsync -zvh rsync://ftp.ensembl.org/ensembl/pub/release-78/variation/VEP/homo_sapiens_vep_78_GRCh{37,38}.tar.gz $VEP_DATA
    cat $VEP_DATA/*.tar.gz | tar -izxf - -C $VEP_DATA

Install the Ensembl v78 API and download the reference FASTAs for GRCh37 and GRCh38:

    cd $VEP_PATH
    perl INSTALL.pl --AUTO af --SPECIES homo_sapiens --ASSEMBLY GRCh37 --DESTDIR $VEP_PATH --CACHEDIR $VEP_DATA
    perl INSTALL.pl --AUTO af --SPECIES homo_sapiens --ASSEMBLY GRCh38 --DESTDIR $VEP_PATH --CACHEDIR $VEP_DATA

Convert the offline cache for use with tabix, that significantly speeds up the lookup of known variants:

    perl convert_cache.pl --species homo_sapiens --version 78_GRCh37 --dir $VEP_DATA
    perl convert_cache.pl --species homo_sapiens --version 78_GRCh38 --dir $VEP_DATA

Test running VEP in offline mode, on the provided sample GRCh37 and GRCh38 VCFs:

    perl variant_effect_predictor.pl --offline --gencode_basic --everything --total_length --allele_number --no_escape --check_existing --xref_refseq --dir $VEP_DATA --fasta $VEP_DATA/homo_sapiens/78_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --assembly GRCh37 --input_file example_GRCh37.vcf --output_file example_GRCh37.vep.txt
    perl variant_effect_predictor.pl --offline --gencode_basic --everything --total_length --allele_number --no_escape --check_existing --xref_refseq --dir $VEP_DATA --fasta $VEP_DATA/homo_sapiens/78_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --assembly GRCh38 --input_file example_GRCh38.vcf --output_file example_GRCh38.vep.txt

Install snpEff
--------------

snpEff ([snpeff.sourceforge.net](http://snpeff.sourceforge.net/)) is popular because of its portability and speed at mapping effects on all possible transcripts in a database like [Ensembl](http://useast.ensembl.org/Homo_sapiens/Info/Annotation) or [Refseq](http://www.ncbi.nlm.nih.gov/refseq/). It's download-able as a java archive, so make sure you have [Java installed](https://www.java.com/en/download/help/download_options.xml).

Create temporary shell variables pointing to where we'll store snpEff and its cache data (non default paths can be used, but specify `--snpeff-path` and `--snpeff-data` when running vcf2maf):

    export SNPEFF_PATH=~/snpEff
    export SNPEFF_DATA=~/snpEff/data

Download the latest release of snpEff into your home directory:

    mkdir $SNPEFF_PATH; cd $SNPEFF_PATH/..
    curl -LO http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
    unzip snpEff_latest_core.zip

Import the Ensembl v75 (Gencode v19) database for GRCh37, and Ensembl v78 (Gencode v21) for GRCh38 (writes to `snpEff/data` by default):

    cd $SNPEFF_PATH
    java -Xmx2g -jar snpEff.jar download GRCh37.75
    java -Xmx2g -jar snpEff.jar download GRCh38.78

Test running snpEff on any available GRCh37 and GRCh38 VCFs:

    java -Xmx4g -jar snpEff.jar eff -sequenceOntology -hgvs GRCh37.75 ~/vep/example_GRCh37.vcf > example_GRCh37.snpeff.vcf
    java -Xmx4g -jar snpEff.jar eff -sequenceOntology -hgvs GRCh38.78 ~/vep/example_GRCh38.vcf > example_GRCh38.snpeff.vcf

Authors
-------

    Cyriac Kandoth (ckandoth@gmail.com)

License
-------

    Apache-2.0 | Apache License, Version 2.0 | https://www.apache.org/licenses/LICENSE-2.0
