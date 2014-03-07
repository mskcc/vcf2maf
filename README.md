vcf2maf
=======

To convert a VCF into a MAF, each variant must be annotated to only one of all possible gene transcripts/isoforms that it might affect. This selection of a single affected transcript/isoform per variant, is often subjective. For now, we try to follow best-practices, but over time, the selection process will be made smarter and more configurable.

Quick start
-----------

Download the script, and view the detailed usage manual:

    curl -LO https://github.com/ckandoth/vcf2maf/archive/master.zip; unzip master.zip
    perl vcf2maf-master/vcf2maf.pl --man

If you don't already have snpEff, see the next section. And then test the script:

    perl vcf2maf-master/vcf2maf.pl --input-vcf vcf2maf-master/test.vcf --output-maf vcf2maf-master/test.maf

Install snpEff
--------------

This script needs snpEff ([snpeff.sourceforge.net](http://snpeff.sourceforge.net/)), a variant annotator that can quickly map each variant to all possible transcripts in a database. It also includes a downloader/importer for thousands of popular transcript databases like from Ensembl and UCSC. snpEff is downloadable as a java archive, so make sure you also have Java installed.

Get to your home directory, and download snpEff's jar files and supporting scripts:

    cd ~/
    curl -LO http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
    unzip snpEff_latest_core.zip
    cd snpEff

There are over 2500 pre-indexed transcript databases available through snpEff (different species/references/versions). They can be listed out as follows:

    java -jar snpEff.jar databases

Grep out the databases for humans:

    java -jar snpEff.jar databases | grep -i homo_sapiens

Download and import the Ensembl v74 database for humans, and also the UCSC Known Genes database:

    java -jar snpEff.jar download GRCh37.74
    java -jar snpEff.jar download hg19kg

To test out snpEff, download and unzip VCFs from the NHLBI [Exome Sequencing Project](http://evs.gs.washington.edu/EVS/):

    curl -LO http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.dbSNP138-rsIDs.snps_indels.vcf.tar.gz
    mkdir nhlbi_esp
    tar -zxvf ESP6500SI-V2-SSA137.dbSNP138-rsIDs.snps_indels.vcf.tar.gz -C nhlbi_esp

Benchmark snpEff on the chr17 NHLBI ESP VCF with ~120K variants (Give it 2GB JVM heap size, use the Ensembl v74 database, and use HGVS output format for codon changes):

    time java -Xmx2g -jar snpEff.jar eff -hgvs GRCh37.74 nhlbi_esp/ESP6500SI-V2-SSA137.updatedRsIds.chr17.snps_indels.vcf > nhlbi_esp/ESP6500SI-V2-SSA137.updatedRsIds.chr17.snps_indels.anno.vcf

It may report some warnings, but that's OK. The Ensembl transcript database is very comprehensive and lists many putative transcripts that need more curation. But our script will appropriately prioritize the "best" transcript to annotate each variant to.

Authors
-------

    Cyriac Kandoth (ckandoth@gmail.com)
    William Lee (leew1@cbio.mskcc.org)

License
-------

    LGPLv3, Memorial Sloan-Kettering Cancer Center, New York, NY 10065, USA
