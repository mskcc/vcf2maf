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
ok( system( "docker run --rm vcf2maf:main perl maf2vcf.pl --help > /dev/null" ) == 0 );
ok( system( "docker run --rm vcf2maf:main perl maf2vcf.pl --man > /dev/null" ) == 0 );

# Test standard operation, diff, and cleanup
ok( system( "docker run --rm -v $test_dir:/opt/tests vcf2maf:main perl maf2vcf.pl --input-maf tests/test_b38_output.maf --output-dir tests --output-vcf tests/test_b38.new.vcf --ref-fasta tests/Homo_sapiens.GRCh38.dna.chromosome.21.fa" ) == 0 );
ok( system( "bash -c 'diff <(cat tests/test_b38.vcf) <(grep -v ^##reference tests/test_b38.new.vcf)'" ) == 0 );
system( "rm -f tests/test_b38.new.vcf tests/test_b38_output.pairs.tsv" );
