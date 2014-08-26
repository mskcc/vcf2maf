vcf2maf
=======

To convert a [VCF](http://samtools.github.io/hts-specs/) into a [MAF](https://wiki.nci.nih.gov/x/eJaPAQ), each variant must be mapped to only one of all possible gene transcripts/isoforms that it might affect. This selection of a single effect per variant, is often subjective. So this project is an attempt to make the selection criteria smarter, reproducible, and more configurable. And the default criteria must lean towards best practices. Per the current defaults, a single affected transcript is selected per variant, as follows:
 1. Sort effects first by transcript [biotype priority](https://github.com/ckandoth/vcf2maf/blob/master/vcf2maf.pl#L75), then by [effect severity](https://github.com/ckandoth/vcf2maf/blob/master/vcf2maf.pl#L17), and finally by transcript length
 2. Pick the gene on the top of the list (worst-affected), and choose it's [canonical](http://www.ensembl.org/Help/Glossary?id=346) transcript (if you used VEP, that relies on [CCDS](http://www.ncbi.nlm.nih.gov/CCDS/))
 3. If the gene has no canonical transcript tagged (if you used snpEff), choose its longest transcript instead

Quick start
-----------

Download the latest release of vcf2maf, and view the detailed usage manual:

    curl -LO https://github.com/ckandoth/vcf2maf/archive/master.zip; unzip master.zip; cd vcf2maf-master
    perl vcf2maf.pl --man

To download properly versioned releases, [click here](https://github.com/ckandoth/vcf2maf/releases) for a list.

If you don't have [VEP](http://useast.ensembl.org/info/docs/tools/vep/index.html) or [snpEff](http://snpeff.sourceforge.net/) installed, see the sections below. VEP is preferred for it's CLIA-compliant [HGVS formats](http://www.hgvs.org/mutnomen/recs.html), and is used by default. So after installing VEP, you can test the script like so:

    perl vcf2maf.pl --input-vcf data/test.vcf --output-maf data/test.maf

If you'd rather use snpEff, there's an option for that, but make sure you ran snpEff with the `-sequenceOntology` option. In older versions of snpEff, this option was incorrectly spelled as `-sequenceOntolgy`:

    perl vcf2maf.pl --input-vcf data/test.vcf --output-maf data/test.snpeff.maf --use-snpeff

If you already have a VCF annotated with either VEP or snpEff, you can use those directly:

    perl vcf2maf.pl --input-vep data/test.vep.vcf --output-maf data/test.maf
    perl vcf2maf.pl --input-snpeff data/test.snpeff.vcf --output-maf data/test.maf

To fill columns 16 and 17 of the output MAF with tumor/normal sample IDs, and to parse out genotypes and allele counts from the corresponding tumor/normal genotype columns in the VCF:

    perl vcf2maf.pl --input-vcf data/test.vcf --output-maf data/test.maf --tumor-id WD1309 --normal-id NB1308

VCFs from variant callers like [VarScan](http://varscan.sourceforge.net/somatic-calling.html#somatic-output) use hardcoded sample IDs TUMOR/NORMAL in the genotype columns of the VCF. To have this script correctly parse the correct genotype columns, while still printing the proper IDs in the output MAF:

    perl vcf2maf.pl --input-vcf data/test.vcf --output-maf data/test.maf --tumor-id WD1309 --normal-id NB1308 --vcf-tumor-id TUMOR --vcf-normal-id NORMAL

If you have VEP in a different folder like `/opt/vep`, and cached in `/srv/vep`, there are options available to point the script there. Similar options available for snpEff too:

    perl vcf2maf.pl --input-vcf data/test.vcf --output-maf data/test.maf --vep-path /opt/vep --vep-data /srv/vep
    perl vcf2maf.pl --input-vcf data/test.vcf --output-maf data/test.maf --snpeff-path /opt/snpEff --snpeff-data /srv/snpEff/data --use-snpeff

Install VEP
-----------

Ensembl's VEP ([Variant Effect Predictor](http://useast.ensembl.org/info/docs/tools/vep/index.html)) is popular for how it selects a single "canonical transcript" per gene as [detailed here](http://useast.ensembl.org/Help/Glossary?id=346), its CLIA-compliant [HGVS variant format](http://www.hgvs.org/mutnomen/recs.html), and [Sequence Ontology nomenclature](http://useast.ensembl.org/info/genome/variation/predicted_data.html#consequences) for variant effects. It's download-able as a Perl script, so make sure you have [Perl installed](http://www.perl.org/get.html).

Download the v76 release of VEP into your home directory:

    cd ~/
    curl -LO https://github.com/Ensembl/ensembl-tools/archive/release/76.tar.gz
    tar -zxf 76.tar.gz --starting-file variant_effect_predictor --transform='s|.*/|vep/|g'
    cd vep

To significantly speed up the step where VEP looks up known variants in COSMIC/dbSNP/etc., it's strongly recommended to have [tabix](https://github.com/samtools/tabix) installed, or available in your `$PATH`. Install it with `apt-get install tabix` on Ubuntu/Debian or `yum install tabix` on CentOS/Redhat/Fedora.

Import the Ensembl v76 (Gencode v20) cache for GRCh37 and GRCh38 (writes to `~/.vep` by default). it is packaged for  If you were not able to set up tabix, then skip argument `--CONVERT` below:

    perl INSTALL.pl --AUTO acf --SPECIES homo_sapiens --ASSEMBLY GRCh37 --CONVERT
    perl INSTALL.pl --AUTO acf --SPECIES homo_sapiens --ASSEMBLY GRCh38 --CONVERT

Test running VEP in offline mode, on the provided sample GRCh38 VCF:

    perl variant_effect_predictor.pl --offline --gencode_basic --everything --total_length --allele_number --no_escape --check_existing --xref_refseq --fasta ~/.vep --assembly GRCh38 --input_file example_GRCh38.vcf --output_file example_GRCh38.vep.txt

Install snpEff
--------------

snpEff ([snpeff.sourceforge.net](http://snpeff.sourceforge.net/)) is popular because of its portability and speed at mapping effects on all possible transcripts in a database like [Ensembl](http://useast.ensembl.org/Homo_sapiens/Info/Annotation) or [Refseq](http://www.ncbi.nlm.nih.gov/refseq/). It's download-able as a java archive, so make sure you have [Java installed](https://www.java.com/en/download/help/download_options.xml).

Download the latest release of snpEff into your home directory:

    cd ~/
    curl -LO http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
    unzip snpEff_latest_core.zip
    cd snpEff

Import the Ensembl v75 (Gencode v19) database for humans (writes to `snpEff/data` by default):

    java -Xmx4g -jar snpEff.jar download GRCh37.75

Test running snpEff on any available sample GRCh37 VCF:

    java -Xmx4g -jar snpEff.jar eff -sequenceOntology -hgvs GRCh37.75 ~/vep/example_GRCh37.vcf > example_GRCh37.snpeff.vcf

Authors
-------

    Cyriac Kandoth (ckandoth@gmail.com)

License
-------

    LGPLv3, Memorial Sloan Kettering Cancer Center, New York, NY 10065, USA
