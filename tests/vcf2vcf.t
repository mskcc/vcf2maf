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
use Test::Simple tests => 4;
ok( system( "docker run --rm vcf2maf:main perl vcf2vcf.pl --help > /dev/null" ) == 0 );
ok( system( "docker run --rm vcf2maf:main perl vcf2vcf.pl --man > /dev/null" ) == 0 );
ok( system( "docker run --rm -v $test_dir:/opt/tests vcf2maf:main perl vcf2vcf.pl --input-vcf tests/test_b38.vcf --output-vcf tests/test_b37.new.vcf --remap-chain data/GRCh38_to_GRCh37.chain --ref-fasta tests/Homo_sapiens.GRCh38.dna.chromosome.21.fa" ) == 0 );
ok( system( "bash -c 'diff <(grep -v ^##fileDate tests/test_b37.vcf) <(grep -v ^##fileDate tests/test_b37.new.vcf)'" ) == 0 );
system( "rm -f tests/test_b37.new.vcf" );
