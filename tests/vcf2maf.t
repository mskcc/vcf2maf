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
ok( system( "perl vcf2maf.pl --help > /dev/null" ) == 0 );
ok( system( "perl vcf2maf.pl --man > /dev/null" ) == 0 );

# Test standard operation, diff, and cleanup
ok( system( "perl vcf2maf.pl --online --ncbi-build GRCh38 --input-vcf tests/test_b38.vcf --output-maf tests/test_b38_output.online.new.maf --ref-fasta data/Homo_sapiens.GRCh38.dna.chromosome.13.fa" ) == 0 );
# Skip column 76 (DOMAINS) in the diff, since it is a randomly ordered comma-delimited list
ok( system( "bash -c 'diff <(cut -f1-75,77- tests/test_b38_output.online.maf) <(cut -f1-75,77- tests/test_b38_output.online.new.maf)'" ) == 0 );
system( "rm -f tests/test_b38_output.online.new.maf" );

# Test some more options, diff, and cleanup
ok( system( "perl vcf2maf.pl --online --vcf-tumor-id TUMOR --vcf-normal-id NORMAL --tumor-id MSK_T001 --normal-id MSK_N001 --maf-center mskcc.org --vep-forks 1 --buffer-size 50 --ncbi-build GRCh38 --input-vcf tests/test_b38.vcf --output-maf tests/test_b38_output.more.new.maf --ref-fasta data/Homo_sapiens.GRCh38.dna.chromosome.13.fa --custom-enst data/isoform_overrides_uniprot --retain-fmt GT" ) == 0 );
ok( system( "bash -c 'diff <(cut -f1-75,77- tests/test_b38_output.more.maf) <(cut -f1-75,77- tests/test_b38_output.more.new.maf)'" ) == 0 );
system( "rm -f tests/test_b38_output.more.new.maf" );
