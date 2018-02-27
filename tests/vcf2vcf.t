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
ok( system( "perl vcf2vcf.pl --help > /dev/null" ) == 0 );
ok( system( "perl vcf2vcf.pl --man > /dev/null" ) == 0 );
ok( system( "perl vcf2vcf.pl --input-vcf tests/test.vcf --output-vcf tests/test_grch38.new.vcf --remap-chain data/GRCh37_to_GRCh38.chain" ) == 0 );
ok( system( "bash -c 'diff <(grep -v ^##fileDate tests/test_grch38.vcf) <(grep -v ^##fileDate tests/test_grch38.new.vcf)'" ) == 0 );

# Cleanup
system( "rm -f tests/test_grch38.new.vcf" );
