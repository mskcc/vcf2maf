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
use Test::Simple tests => 3;
ok( system( "perl maf2vcf.pl --help > /dev/null" ) == 0 );
ok( system( "perl maf2vcf.pl --man > /dev/null" ) == 0 );
ok( system( "perl maf2vcf.pl --input-maf tests/test.maf --output-dir tests/vcfs --per-tn-vcfs" ) == 0 );
