vcf<img src="https://i.giphy.com/R6X7GehJWQYms.gif" width="28">maf
=======

To convert a [VCF](https://samtools.github.io/hts-specs//) into a [MAF](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format), each variant must be mapped to only one of all possible gene transcripts/isoforms that it might affect. But even within a single isoform, a `Missense_Mutation` close enough to a `Splice_Site`, can be labeled as either in MAF format, but not as both. **This selection of a single effect per variant, is often subjective. And that's what this project attempts to standardize.** The `vcf2maf` and `maf2maf` scripts leave most of that responsibility to [Ensembl's VEP](http://ensembl.org/info/docs/tools/vep/index.html), but allows you to override their "canonical" isoforms, or use a custom ExAC VCF for annotation. Though the most useful feature is the **extensive support in parsing a wide range of crappy MAF-like or VCF-like formats** we've seen out in the wild.

Quick start
-----------

Find the [latest release](https://github.com/mskcc/vcf2maf/releases), download it, and view the detailed usage manuals for `vcf2maf` and `maf2maf`:

    export VCF2MAF_URL=`curl -sL https://api.github.com/repos/mskcc/vcf2maf/releases | grep -m1 tarball_url | cut -d\" -f4`
    curl -L -o mskcc-vcf2maf.tar.gz $VCF2MAF_URL; tar -zxf mskcc-vcf2maf.tar.gz; cd mskcc-vcf2maf-*
    perl vcf2maf.pl --man
    perl maf2maf.pl --man

If you don't have VEP installed, then [follow this gist](https://gist.github.com/ckandoth/4bccadcacd58aad055ed369a78bf2e7c). Of the many annotators out there, VEP is preferred for its large team of active coders, and its CLIA-compliant [HGVS formats](http://www.hgvs.org/mutnomen/recs.html). After installing VEP, test out `vcf2maf` like this:

    perl vcf2maf.pl --input-vcf tests/test.vcf --output-maf tests/test.vep.maf

To fill columns 16 and 17 of the output MAF with tumor/normal sample IDs, and to parse out genotypes and allele counts from matched genotype columns in the VCF, use options `--tumor-id` and `--normal-id`. Skip option `--normal-id` if you didn't have a matched normal:

    perl vcf2maf.pl --input-vcf tests/test.vcf --output-maf tests/test.vep.maf --tumor-id WD1309 --normal-id NB1308

VCFs from variant callers like [VarScan](http://varscan.sourceforge.net/somatic-calling.html#somatic-output) use hardcoded sample IDs TUMOR/NORMAL to name genotype columns. To have `vcf2maf` correctly locate the columns to parse genotypes, while still printing proper sample IDs in the output MAF:

    perl vcf2maf.pl --input-vcf tests/test_varscan.vcf --output-maf tests/test_varscan.vep.maf --tumor-id WD1309 --normal-id NB1308 --vcf-tumor-id TUMOR --vcf-normal-id NORMAL

If VEP is installed under `/opt/vep` and the VEP cache is under `/srv/vep`, there are options available to tell `vcf2maf` where to find them:

    perl vcf2maf.pl --input-vcf tests/test.vcf --output-maf tests/test.vep.maf --vep-path /opt/vep --vep-data /srv/vep

If you want to skip running VEP and need a minimalist MAF-like file listing data from the input VCF only, then use the `--inhibit-vep` option. If your input VCF contains VEP annotation, then `vcf2maf` will try to extract it. But be warned that the accuracy of your resulting MAF depends on how VEP was operated upstream. In standard operation, `vcf2maf` runs VEP with very specific parameters to make sure everyone produces comparable MAFs. So, it is strongly recommended to avoid `--inhibit-vep` unless you know what you're doing.

maf2maf
-------

If you have a MAF or a MAF-like file that you want to reannotate, then use `maf2maf`, which simply runs `maf2vcf` followed by `vcf2maf`:

    perl maf2maf.pl --input-maf tests/test.maf --output-maf tests/test.vep.maf

After tests on variant lists from many sources, `maf2vcf` and `maf2maf` are quite good at dealing with formatting errors or "MAF-like" files. It even supports VCF-style alleles, as long as `Start_Position == POS`. But it's OK if the input format is imperfect. Any variants with a reference allele mismatch are kept aside in a separate file for debugging. The bare minimum columns that `maf2maf` expects as input are:

    Chromosome	Start_Position	Reference_Allele	Tumor_Seq_Allele2	Tumor_Sample_Barcode
    1	3599659	C	T	TCGA-A1-A0SF-01
    1	6676836	A	AGC	TCGA-A1-A0SF-01
    1	7886690	G	A	TCGA-A1-A0SI-01

See `data/minimalist_test_maf.tsv` for a sampler. Addition of `Tumor_Seq_Allele1` will be used to determine zygosity. Otherwise, it will try to determine zygosity from variant allele fractions, assuming that arguments `--tum-vad-col` and `--tum-depth-col` are set correctly to the names of columns containing those read counts. Specifying the `Matched_Norm_Sample_Barcode` with its respective columns containing read-counts, is also strongly recommended. Columns containing normal allele read counts can be specified using argument `--nrm-vad-col` and `--nrm-depth-col`.

Docker
------

Assuming you have a recent version of docker, clone the main branch and build an image as follows:

    git clone git@github.com:mskcc/vcf2maf.git
    cd vcf2maf
    docker build -t vcf2maf:main .
    docker builder prune -f

Now you run the scripts in docker as follows:

    docker run --rm vcf2maf:main perl vcf2maf.pl --help
    docker run --rm vcf2maf:main perl maf2maf.pl --help

Testing
-------

A small standalone test dataset was created by restricting VEP v112 cache/fasta to chr21 in GRCh38 and hosting that on a private server for download by CI services. We can manually fetch those as follows:

    wget -P tests https://data.cyri.ac/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
    gzip -d tests/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
    wget -P tests https://data.cyri.ac/homo_sapiens_vep_112_GRCh38_chr21.tar.gz
    tar -zxf tests/homo_sapiens_vep_112_GRCh38_chr21.tar.gz -C tests

And the following scripts test the docker image on predefined inputs and compare outputs against expected outputs:

    perl tests/vcf2maf.t
    perl tests/vcf2vcf.t
    perl tests/maf2vcf.t

License
-------

    Apache-2.0 | Apache License, Version 2.0 | https://www.apache.org/licenses/LICENSE-2.0

Citation
--------

    Cyriac Kandoth. mskcc/vcf2maf: vcf2maf v1.6. (2020). doi:10.5281/zenodo.593251
