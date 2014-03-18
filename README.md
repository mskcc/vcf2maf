vcf2maf
=======

To convert a [VCF](http://samtools.github.io/hts-specs/) into a [MAF](https://wiki.nci.nih.gov/x/eJaPAQ), each variant must be annotated to only one of all possible gene transcripts/isoforms that it might affect. This selection of a single effect per variant, is often subjective. So this project is an attempt to make the selection criteria smarter, reproducible, and more configurable. And the default criteria must lean toward best practices.

Quick start
-----------

Download the latest release, and view the detailed usage manual:

    curl -LO https://github.com/ckandoth/vcf2maf/archive/master.zip; unzip master.zip; cd vcf2maf-master
    perl vcf2maf.pl --man

To download properly versioned releases, [click here](https://github.com/ckandoth/vcf2maf/releases) for a list.

If you don't have [VEP](http://useast.ensembl.org/info/docs/tools/vep/index.html) or [snpEff](http://snpeff.sourceforge.net/) installed, follow the sections below. VEP is preferred for it's CLIA-compliant HGVS format, and is used by default. So after installing VEP following the section below, you can test the script on the provided `test.vcf`:

    perl vcf2maf.pl --input-vcf test.vcf --output-maf test.maf

If you had VEP installed in a different folder like `/opt/vep`, and cached in `/srv/vep`, then point the script there:

    perl vcf2maf.pl --input-vcf test.vcf --output-maf test.maf --vep-path /opt/vep --vep-data /srv/vep

If you'd rather run snpEff, which runs much faster than VEP, there's an option for that:

    perl vcf2maf.pl --input-vcf test.vcf --output-maf test.snpeff.maf --use-snpeff

If you already annotated a VCF with either VEP or snpEff, you can use those directly:

    perl vcf2maf.pl --input-vep test.vep.vcf --output-maf test.maf
    perl vcf2maf.pl --input-snpeff test.snpeff.vcf --output-maf test.maf

Install VEP
-----------

Ensembl's VEP ([Variant Effect Predictor](http://useast.ensembl.org/info/docs/tools/vep/index.html)) is popular for how it selects a single "canonical transcript" per gene as [detailed here](http://useast.ensembl.org/Help/Glossary?id=346), its CLIA-compliant [HGVS variant format](http://www.hgvs.org/mutnomen/recs.html), and [Sequence Ontology nomenclature](http://useast.ensembl.org/info/genome/variation/predicted_data.html#consequences) for variant effects. It's download-able as a Perl script, so make sure you have [Perl installed](http://www.perl.org/get.html).

Download the latest release of VEP into your home directory:

    cd ~/
    curl -LO https://github.com/Ensembl/ensembl-tools/archive/release/75.tar.gz
    tar -zxf 75.tar.gz --starting-file variant_effect_predictor --transform='s|.*/|vep/|g'
    cd vep

Import the Ensembl v75 (Gencode v19) database for humans:

    perl INSTALL.pl --AUTO acf --SPECIES homo_sapiens_vep --VERSION 75

Test running VEP in offline mode with 4 parallel threads, on the provided example VCF:

    perl variant_effect_predictor.pl --offline --no_stats --everything --xref_refseq --check_existing --total_length --allele_number --no_escape --fork 4 --fasta ~/.vep --input_file example.vcf --output_file example.vep.txt

Install snpEff
--------------

snpEff ([snpeff.sourceforge.net](http://snpeff.sourceforge.net/)) is popular because of its portability and speed at mapping effects on all possible transcripts in a database like [Ensembl](http://useast.ensembl.org/Homo_sapiens/Info/Annotation) or [Refseq](http://www.ncbi.nlm.nih.gov/refseq/). It's download-able as a java archive, so make sure you have [Java installed](https://www.java.com/en/download/help/download_options.xml).

Download the latest release of snpEff into your home directory:

    cd ~/
    curl -LO http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
    unzip snpEff_latest_core.zip
    cd snpEff

Import the Ensembl v75 (Gencode v19) database for humans:

    java -Xmx4g -jar snpEff.jar download GRCh37.75

Test running snpEff on a sample VCF:

    java -Xmx4g -jar snpEff.jar eff -noStats -sequenceOntolgy -hgvs GRCh37.75 example.vcf > example.snpeff.vcf

Authors
-------

    Cyriac Kandoth (ckandoth@gmail.com)
    William Lee, Senior Research Scientist, Memorial Sloan Kettering Cancer Center

License
-------

    LGPLv3, Memorial Sloan Kettering Cancer Center, New York, NY 10065, USA
