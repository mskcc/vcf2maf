#!/usr/bin/env perl

use strict;
use warnings;

# Find and chdir into the parent folder containing the scripts to test
use File::Basename 'dirname';
use Cwd 'abs_path';
my $test_dir = dirname( abs_path( __FILE__ ));
my $script_dir = dirname( $test_dir );
chdir $script_dir;

# Set the number of tests we'll run, and run them
use Test::Simple tests => 6;
ok( system( "docker run --rm vcf2maf:main perl vcf2maf.pl --help > /dev/null" ) == 0 );
ok( system( "docker run --rm vcf2maf:main perl vcf2maf.pl --man > /dev/null" ) == 0 );

# Test standard operation, diff, and cleanup
ok( system( "docker run --rm -v $test_dir:/opt/tests vcf2maf:main perl vcf2maf.pl --vep-path /usr/local/bin --vep-data tests --ncbi-build GRCh38 --input-vcf tests/test_b38.vcf --output-maf tests/test_b38_output.new.maf --ref-fasta tests/Homo_sapiens.GRCh38.dna.chromosome.21.fa" ) == 0 );
# Skip column 76 (DOMAINS) in the diff, since it is a randomly ordered comma-delimited list
ok( system( "bash -c 'diff <(cut -f1-75,77- tests/test_b38_output.maf) <(cut -f1-75,77- tests/test_b38_output.new.maf)'" ) == 0 );
system( "rm -f tests/test_b38_output.new.maf" );

# Test some more options, diff, and cleanup
ok( system( "docker run --rm -v $test_dir:/opt/tests vcf2maf:main perl vcf2maf.pl --vep-path /usr/local/bin --vep-data tests --ncbi-build GRCh38 --input-vcf tests/test_b38.vcf --output-maf tests/test_b38_output.more.new.maf --ref-fasta tests/Homo_sapiens.GRCh38.dna.chromosome.21.fa --vep-overwrite --vcf-tumor-id TUMOR --vcf-normal-id NORMAL --tumor-id MSK_T001 --normal-id MSK_N001 --maf-center mskcc.org --vep-forks 1 --buffer-size 50 --vep-custom tests/test_b38.gnomad.exomes.r2.1.1.sites.vcf.gz,gnomAD,vcf,exact,,AC --retain-ann gnomAD_AC --retain-fmt GT" ) == 0 );
ok( system( "bash -c 'diff <(cut -f1-75,77- tests/test_b38_output.more.maf) <(cut -f1-75,77- tests/test_b38_output.more.new.maf)'" ) == 0 );
system( "rm -f tests/test_b38_output.more.new.maf tests/test_b38.vep.vcf" );
