#!/usr/bin/env perl

# vcf2vcf - Create a minimalist VCF from a given VCF, fixing various problems in the process

use strict;
use warnings;
use IO::File;
use Getopt::Long qw( GetOptions );
use Pod::Usage qw( pod2usage );
use POSIX qw( strftime );

# Set any default paths and constants
my ( $vcf_tumor_id, $vcf_normal_id ) = ( "TUMOR", "NORMAL" );

# Define the minimal INFO and FORMAT fields that we want to retain
my @retain_info = qw( SOMATIC SS );
my @retain_format = qw( GT DP AD );

# Define FILTER descriptors that we'll add if user specified the --add-filters option
my ( $min_tum_depth, $min_nrm_depth ) = ( 14, 8 );
my ( $min_tum_support, $max_nrm_support ) = ( 3, 2 );
my %filter_tags = (
    "LowTotalDepth" => "Less than $min_tum_depth total reads in tumor, or less than $min_nrm_depth reads in normal",
    "LowTumorSupport" => "Less than $min_tum_support allele supporting read(s) in tumor",
    "HighNormalSupport" => "More than $max_nrm_support allele supporting read(s) in the normal"
);

# Check for missing or crappy arguments
unless( @ARGV and $ARGV[0] =~ m/^-/ ) {
    pod2usage( -verbose => 0, -message => "ERROR: Missing or invalid arguments!\n", -exitval => 2 );
}

# Parse options and print usage if there is a syntax error, or if usage was explicitly requested
my ( $man, $help, $add_filters ) = ( 0, 0, 0 );
my ( $input_vcf, $output_vcf, $new_tumor_id, $new_normal_id );
GetOptions(
    'help!' => \$help,
    'man!' => \$man,
    'input-vcf=s' => \$input_vcf,
    'output-vcf=s' => \$output_vcf,
    'vcf-tumor-id=s' => \$vcf_tumor_id,
    'vcf-normal-id=s' => \$vcf_normal_id,
    'new-tumor-id=s' => \$new_tumor_id,
    'new-normal-id=s' => \$new_normal_id,
    'add-filters!' => \$add_filters
) or pod2usage( -verbose => 1, -input => \*DATA, -exitval => 2 );
pod2usage( -verbose => 1, -input => \*DATA, -exitval => 0 ) if( $help );
pod2usage( -verbose => 2, -input => \*DATA, -exitval => 0 ) if( $man );

# Make sure we were given an inout VCF file
die "ERROR: Please specify an input file with --input-vcf. STDIN is not supported.\n" unless( $input_vcf );

# Unless specified, assume that the new VCF will use the same sample IDs as the input VCF
$new_tumor_id = $vcf_tumor_id unless( $new_tumor_id );
$new_normal_id = $vcf_normal_id unless( $new_normal_id );

# Initialize the output VCF file with the mandatory first line of a VCF header
my $vcf_out_fh = *STDOUT; # Use STDOUT if an output VCF file was not defined
$vcf_out_fh = IO::File->new( $output_vcf, ">" ) or die "ERROR: Couldn't open output VCF file: $output_vcf!\n" if( $output_vcf );
$vcf_out_fh->print( "##fileformat=VCFv4.2\n" );

# Also add today's date just cuz we can
my $file_date = strftime( "%Y%m%d", localtime );
$vcf_out_fh->print( "##fileDate=$file_date\n" );

# And add our custom FILTER tag descriptors, if --add-filters was specified
if( $add_filters ) {
    $vcf_out_fh->print( "##FILTER=<ID=$_,Description=\"" . $filter_tags{$_} . "\">\n" ) foreach ( sort keys %filter_tags );
}

# Parse through each variant in the annotated VCF, and fix 'em
my $vcf_in_fh = IO::File->new( $input_vcf ) or die "ERROR: Couldn't open input VCF file: $input_vcf!\n";
my ( $vcf_tumor_idx, $vcf_normal_idx, %vcf_header_lines );
while( my $line = $vcf_in_fh->getline ) {

    # Skip all but a few INFO/FORMAT/FILTER descriptors in the header
    if( $line =~ m/^##/ ) {
        # Retain only the minimal INFO amd FORMAT descriptor lines in the header
        if( $line =~ m/^##(INFO|FORMAT)=<ID=([^,]+)/ ) {
            my ( $type, $tag ) = ( $1, $2 );
            if(( $type eq "INFO" and grep( $tag, @retain_info )) or ( $type eq "FORMAT" and grep( $tag, @retain_format ))) {
                $vcf_out_fh->print( $line );
            }
        }
        elsif( $line =~ m/^##FILTER=<ID=([^,]+)/ and !defined $filter_tags{$1}) {
            $vcf_out_fh->print( $line );
        }
        next;
    }

    chomp( $line );
    my ( $chrom, $pos, $ids, $ref, $alt, $qual, $filter, $info_line, $format_line, @rest ) = split( /\t/, $line );

    # If FORMATted genotype fields are available, find the sample with the variant, and matched normal
    if( $line =~ m/^#CHROM/ ) {
        if( $format_line and scalar( @rest ) > 0 ) {
            for( my $i = 0; $i <= $#rest; ++$i ) {
                $vcf_tumor_idx = $i if( $rest[$i] eq $vcf_tumor_id );
                $vcf_normal_idx = $i if( $rest[$i] eq $vcf_normal_id );
            }
            ( defined $vcf_tumor_idx ) or warn "WARNING: No genotype column for $vcf_tumor_id in VCF!\n";
            ( defined $vcf_normal_idx ) or warn "WARNING: No genotype column for $vcf_normal_id in VCF!\n";
        }
        $vcf_out_fh->print( "#", join( "\t", qw( CHROM POS ID REF ALT QUAL FILTER INFO FORMAT ), $new_tumor_id, $new_normal_id ), "\n" );
        next;
    }

    # Set QUAL and FILTER to "." unless defined and non-empty
    $qual = "." unless( defined $qual and $qual ne "" );
    $filter = "." unless( defined $filter and $filter ne "" );

    # Parse out the data in the info column, and retain only the minimal fields
    my %info = map {( m/=/ ? ( split( /=/, $_, 2 )) : ( $_, "" ))} split( /\;/, $info_line );
    $info_line = join( ";", map{( $info{$_} eq "" ? "$_" : "$_=$info{$_}" )} grep{defined $info{$_}} @retain_info );

    # By default, the variant allele is the first (usually the only) allele listed under ALT
    # If there are multiple ALT alleles, choose the allele specified in the tumor GT field
    # If tumor GT is undefined or ambiguous, choose the one with the most supporting read depth
    my @alleles = ( $ref, split( /,/, $alt ));
    my $var_allele_idx = 1;

    # Parse through info from the tumor genotype field
    my %tum_info;
    if( defined $vcf_tumor_idx ) {
        my @format_keys = split( /\:/, $format_line );
        my $idx = 0;
        %tum_info = map {( $format_keys[$idx++], $_ )} split( /\:/, $rest[$vcf_tumor_idx] );

        # If possible, parse the genotype to identify the variant allele in this tumor
        if( defined $tum_info{GT} and $tum_info{GT} ne "." and $tum_info{GT} ne "./." ) {
            my @genotype = split( /[\/|]/, $tum_info{GT} );
            # In case of polyploid calls, choose the first non-REF allele, if any
            ( $var_allele_idx ) = grep {$_ ne "0"} @genotype;
            $var_allele_idx = 1 unless( defined $var_allele_idx and $var_allele_idx =~ m/^\d+$/ );
        }

        # Standardize AD and DP based on data in the genotype fields
        FixAlleleDepths( \@alleles, $var_allele_idx, \%tum_info );
        $tum_info{GT} = "./." unless( defined $tum_info{GT} and $tum_info{GT} ne '.' );
    }

    # Parse through info from the normal genotype field
    my %nrm_info;
    if( defined $vcf_normal_idx ) {
        my @format_keys = split( /\:/, $format_line );
        my $idx = 0;
        %nrm_info = map {( $format_keys[$idx++], $_ )} split( /\:/, $rest[$vcf_normal_idx] );

        # Standardize AD and DP based on data in the genotype fields
        FixAlleleDepths( \@alleles, $var_allele_idx, \%nrm_info );
        $nrm_info{GT} = "./." unless( defined $nrm_info{GT} and $nrm_info{GT} ne '.' );
    }

    # Add more filter tags to the FILTER field, if --add-filters was specified
    if( $add_filters ) {
        my %tags = ();
        map{ $tags{$_} = 1 unless( $_ eq "PASS" or $_ eq "." )} split( ",", $filter );
        $tags{LowTotalDepth} = 1 if(( $tum_info{DP} ne "." and $tum_info{DP} < $min_tum_depth ) or ( $nrm_info{DP} ne "." and $nrm_info{DP} < $min_nrm_depth ));
        my @tum_depths = split( /,/, $tum_info{AD} );
        my @nrm_depths = split( /,/, $nrm_info{AD} );
        my $tum_alt_depth = $tum_depths[$var_allele_idx];
        my $nrm_alt_depth = $nrm_depths[$var_allele_idx];
        $tags{LowTumorSupport} = 1 if( $tum_alt_depth ne "." and $tum_alt_depth < $min_tum_support );
        $tags{HighNormalSupport} = 1 if( $nrm_alt_depth ne "." and $nrm_alt_depth > $max_nrm_support );
        my $tags_to_add = join( ",", sort keys %tags );
        $filter = ( $tags_to_add ? $tags_to_add : $filter );
    }

    # Print a VCF file with minimal genotype fields for tumor and normal
    $format_line = join( ":", @retain_format );
    my $tum_fmt = join( ":", map{( $tum_info{$_} ? $tum_info{$_} : "." )} @retain_format );
    my $nrm_fmt = join( ":", map{( $nrm_info{$_} ? $nrm_info{$_} : "." )} @retain_format );
    $vcf_out_fh->print( join( "\t", $chrom, $pos, $ids, $ref, $alt, $qual, $filter, $info_line, $format_line, $tum_fmt, $nrm_fmt ), "\n" );
}
$vcf_out_fh->close if( $output_vcf );
$vcf_in_fh->close;

# Fix the AD and DP fields, given data from a FORMATted genotype string
sub FixAlleleDepths {
    my ( $alleles_ref, $var_allele_idx, $fmt_info_ref ) = @_;
    my %fmt_info = %{$fmt_info_ref};
    my @alleles = @{$alleles_ref};
    my @depths = ();

    # If AD is defined, then parse out all REF/ALT allele depths, or whatever is in it
    if( defined $fmt_info{AD} and $fmt_info{AD} ne "." ) {
        @depths = map{( m/^\d+$/ ? $_ : "" )}split( /,/, $fmt_info{AD} );
    }

    # Handle VarScan VCF lines where AD contains only 1 depth, and REF allele depth is in RD
    if( scalar( @depths ) == 1 and defined $fmt_info{RD} ) {
        @depths = map{""} @alleles;
        $depths[0] = $fmt_info{RD};
        $depths[$var_allele_idx] = $fmt_info{AD};
    }
    # Handle SomaticSniper VCF lines, where allele depths must be extracted from BCOUNT
    elsif( !defined $fmt_info{AD} and defined $fmt_info{BCOUNT} ) {
        my %b_idx = ( A=>0, C=>1, G=>2, T=>3 );
        my @bcount = split( /,/, $fmt_info{BCOUNT} );
        @depths = map{(( defined $b_idx{$_} and defined $bcount[$b_idx{$_}] ) ? $bcount[$b_idx{$_}] : "" )} @alleles;
    }
    # Handle VCF SNV lines by Strelka, where allele depths are in AU:CU:GU:TU
    elsif( !defined $fmt_info{AD} and scalar( grep{defined $fmt_info{$_}} qw/AU CU GU TU/ ) == 4 ) {
        # Strelka allele depths come in tiers 1,2. We'll use tier1 cuz it's stricter, and DP already is
        map{( $fmt_info{$_.'U'} ) = split( ",", $fmt_info{$_.'U'} )} qw( A C G T );

        # If the only ALT allele is N, then set it to the allele with the highest non-ref readcount
        if( scalar( @alleles ) == 2 and $alleles[1] eq "N" ) {
            $var_allele_idx = 1;
            my %acgt_depths = map{( defined $fmt_info{$_.'U'} ? ( $_, $fmt_info{$_.'U'} ) : ( $_, "" ))} qw( A C G T );
            my @deepest = sort {$acgt_depths{$b} <=> $acgt_depths{$a}} keys %acgt_depths;
            ( $alleles[$var_allele_idx] ) = ( $deepest[0] ne $alleles[0] ? $deepest[0] : $deepest[1] );
        }
        @depths = map{( defined $fmt_info{$_.'U'} ? $fmt_info{$_.'U'} : "" )} @alleles;
    }
    # Handle VCF Indel lines by Strelka, where variant allele depth is in TIR
    elsif( !defined $fmt_info{AD} and $fmt_info{TIR} ) {
        # Reference allele depth is not provided by Strelka for indels, so we have to skip it
        @depths = ( "", ( split /,/, $fmt_info{TIR} )[0] );
    }
    # Handle VCF lines from the Ion Torrent Suite where ALT depths are in AO and REF depths are in RO
    elsif( !defined $fmt_info{AD} and defined $fmt_info{AO} and defined $fmt_info{RO} ) {
        @depths = ( $fmt_info{RO}, map{( m/^\d+$/ ? $_ : "" )}split( /,/, $fmt_info{AO} ));
    }
    # Handle VCF lines with ALT allele fraction in FA, which needs to be multiplied by DP to get AD
    elsif( !defined $fmt_info{AD} and defined $fmt_info{FA} and defined $fmt_info{DP} and $fmt_info{DP} ne '.' ) {
        # Reference allele depth and depths for any other ALT alleles must be left undefined
        @depths = map{""} @alleles;
        $depths[$var_allele_idx] = sprintf( "%.0f", $fmt_info{FA} * $fmt_info{DP} );
    }
    # Handle VCF lines from mpileup/bcftools where DV contains the ALT allele depth
    elsif( !defined $fmt_info{AD} and defined $fmt_info{DV} and defined $fmt_info{DP} ) {
        # Reference allele depth and depths for any other ALT alleles must be left undefined
        @depths = map{""} @alleles;
        $depths[$var_allele_idx] = $fmt_info{DV};
    }
    # Handle VCF lines where AD contains only 1 value, that we can assume is the variant allele
    elsif( defined $fmt_info{AD} and @depths and scalar( @depths ) == 1 ) {
        # Reference allele depth and depths for any other ALT alleles must be left undefined
        @depths = map{""} @alleles;
        $depths[$var_allele_idx] = $fmt_info{AD};
    }
    # For all other lines where #depths is not equal to #alleles, blank out the depths
    elsif( @depths and scalar( @depths ) ne scalar( @alleles )) {
        @depths = map{""} @alleles;
    }

    # Sanity check that REF/ALT allele depths are lower than the total depth
    if( defined $fmt_info{DP} and $fmt_info{DP} ne '.' and (( $depths[0] and $depths[0] > $fmt_info{DP} ) or
        ( $depths[$var_allele_idx] and $depths[$var_allele_idx] > $fmt_info{DP} ) or
        ( $depths[0] and $depths[$var_allele_idx] and $depths[0] + $depths[$var_allele_idx] > $fmt_info{DP} ))) {
        $fmt_info{DP} = 0;
        map{$fmt_info{DP} += $_ if($_ and $_ ne '.')} @depths;
    }

    # If we have REF/ALT allele depths, but no DP, then set DP equal to the sum of all ADs
    if(( defined $depths[0] and defined $depths[$var_allele_idx] ) and ( !defined $fmt_info{DP} or $fmt_info{DP} eq '.' )) {
        $fmt_info{DP} = 0;
        map{$fmt_info{DP} += $_ if($_ and $_ ne '.')} @depths;
    }

    # Put all our changes back into the hash/array references that were passed over
    $fmt_info{AD} = join( ",", map{( $_ ne "" ? $_ : "." )} @depths );
    %{$fmt_info_ref} = %fmt_info;
    @{$alleles_ref} = @alleles;

    return 1;
}

__DATA__

=head1 NAME

 vcf2vcf.pl - Create a minimalist VCF from a given VCF, fixing various problems in the process

=head1 SYNOPSIS

 perl vcf2vcf.pl --help
 perl vcf2vcf.pl --input-vcf weird.vcf --output-vcf not_weird.vcf --vcf-tumor-id WD4086 --vcf-normal-id NB4086

=head1 OPTIONS

 --input-vcf      Path to input file in VCF format
 --output-vcf     Path to output VCF file [Default: STDOUT]
 --vcf-tumor-id   Tumor sample ID used in VCF's genotype column [Default: TUMOR]
 --vcf-normal-id  Matched normal ID used in VCF's genotype column [Default: NORMAL]
 --new-tumor-id   Tumor sample ID to use in the new VCF [Default: --vcf-tumor-id]
 --new-normal-id  Matched normal ID to use in the new VCF [Default: --vcf-normal-id]
 --add-filters    Use this to add some extra tags under FILTER [Default: 0]
 --help           Print a brief help message and quit
 --man            Print the detailed manual

=head1 DESCRIPTION

The VCF files that variant callers generate are rarely compliant with VCF specifications. This script fixes the most serious grievances, and creates a VCF with only important fields in INFO and FORMAT, like GT:DP:AD. Specify the --add-filters option, which will use the allele depths and fractions to add more tags under FILTER.

=head2 Relevant links:

 Homepage: https://github.com/ckandoth/vcf2maf
 VCF format: http://samtools.github.io/hts-specs/

=head1 AUTHORS

 Cyriac Kandoth (ckandoth@gmail.com)

=head1 LICENSE

 Apache-2.0 | Apache License, Version 2.0 | https://www.apache.org/licenses/LICENSE-2.0

=cut
