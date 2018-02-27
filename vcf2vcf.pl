#!/usr/bin/env perl

# vcf2vcf - Create a standardized VCF from a given VCF, fixing various problems in the process

use strict;
use warnings;
use IO::File;
use Getopt::Long qw( GetOptions );
use Pod::Usage qw( pod2usage );
use POSIX qw( strftime );
use File::Temp qw( tempdir );

# Set any default paths and constants
my ( $vcf_tumor_id, $vcf_normal_id ) = ( "TUMOR", "NORMAL" );
my ( $add_header, $add_info ) = ( "", "" );
my ( $retain_info, $retain_format ) = ( "SOMATIC,SS,I16,MQSB", "GT,AD,DP" );

# Find out if samtools and bcftools are properly installed, and warn the user if it's not
my ( $samtools ) = map{chomp; $_}`which samtools`;
( $samtools and -e $samtools ) or die "ERROR: Please install samtools, and make sure it's in your PATH\n";
my ( $bcftools ) = map{chomp; $_}`which bcftools`;
( $bcftools and -e $bcftools ) or die "ERROR: Please install bcftools, and make sure it's in your PATH\n";

# E.g. Get DP:AD from tumor/normal BAMs for a TGG-insertion at GRCh38 loci 21:46302071-46302072:
# ::NOTE:: This returns multiple lines+alleles, so the matching allele needs to be parsed out
# samtools mpileup --region chr21:46302071-46302071 --count-orphans --no-BAQ --min-MQ 1 --min-BQ 5 --ignore-RG --excl-flags UNMAP,SECONDARY,QCFAIL,DUP --VCF --uncompressed --output-tags DP,AD,ADF,ADR --ext-prob 20 --gap-frac 0.002 --tandem-qual 80 --min-ireads 1 --open-prob 30 --fasta-ref /ifs/depot/assemblies/H.sapiens/GRCh38_GDC/GRCh38.d1.vd1.fa /ifs/tcga/gdc/files/95a08b09-c46a-458c-bd21-deb43a309b00/69f9c49d8f6376a7092cff2a3bd2922b_gdc_realn.bam /ifs/tcga/gdc/files/47ac9742-74bc-4d76-a2ac-46c708e9cbbd/e30ade9704fbc29ccd9e6b69c91db237_gdc_realn.bam

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
my ( $input_vcf, $output_vcf, $new_tumor_id, $new_normal_id, $remap_chain );
my ( $tumor_bam, $normal_bam, $ref_fasta ) = ( "", "", "$ENV{HOME}/.vep/homo_sapiens/91_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz" );
GetOptions(
    'help!' => \$help,
    'man!' => \$man,
    'input-vcf=s' => \$input_vcf,
    'output-vcf=s' => \$output_vcf,
    'vcf-tumor-id=s' => \$vcf_tumor_id,
    'vcf-normal-id=s' => \$vcf_normal_id,
    'new-tumor-id=s' => \$new_tumor_id,
    'new-normal-id=s' => \$new_normal_id,
    'tumor-bam=s' => \$tumor_bam,
    'normal-bam=s' => \$normal_bam,
    'ref-fasta=s' => \$ref_fasta,
    'add-header=s' => \$add_header,
    'add-info=s' => \$add_info,
    'retain-info=s' => \$retain_info,
    'retain-format=s' => \$retain_format,
    'remap-chain=s' => \$remap_chain,
    'add-filters!' => \$add_filters
) or pod2usage( -verbose => 1, -input => \*DATA, -exitval => 2 );
pod2usage( -verbose => 1, -input => \*DATA, -exitval => 0 ) if( $help );
pod2usage( -verbose => 2, -input => \*DATA, -exitval => 0 ) if( $man );

# Make sure we are given paths to the input and output VCF files, and a matched reference FASTA
die "ERROR: Please specify an input file with --input-vcf. STDIN is not supported.\n" unless( $input_vcf );
die "ERROR: Please specify an output file with --output-vcf. STDOUT is not supported.\n" unless( $output_vcf );
die "ERROR: Provided --ref-fasta is missing or empty: $ref_fasta\n" unless( -s $ref_fasta );

# Unless specified, assume that the new VCF will use the same sample IDs as the input VCF
$new_tumor_id = $vcf_tumor_id unless( $new_tumor_id );
$new_normal_id = $vcf_normal_id unless( $new_normal_id );

# To grep them easily, make vectors containing the INFO/FORMAT tags
my @retain_info_tags = split( ",", $retain_info );
my @retain_format_tags = split( ",", $retain_format );

# If user defined additional INFO tags, then wipe clean any whitespace, and make it ";" delimited
if( $add_info ) {
    $add_info =~ s/\s//g;
    $add_info =~ s/,/;/g;
}

# If needed, create a temporary folder that will get deleted after a clean exit
my $tmp_dir = tempdir( CLEANUP => 1 ) if( $remap_chain or $tumor_bam or $normal_bam );

# If a liftOver chain was provided, remap and switch the input VCF before proceeding
my %remap;
if( $remap_chain ) {
    # Find out if liftOver is properly installed, and warn the user if it's not
    my $liftover = `which liftOver`;
    chomp( $liftover );
    ( $liftover and -e $liftover ) or die "ERROR: Please install liftOver, and make sure it's in your PATH\n";

    # Make a BED file from the VCF, run liftOver on it, and create a hash mapping old to new loci
    `grep -v ^# $input_vcf | cut -f1,2 | awk '{OFS="\\t"; print \$1,\$2-1,\$2,\$1":"\$2}' > $tmp_dir/input.bed`;
    %remap = map{chomp; my @c=split("\t"); ($c[3], "$c[0]:$c[2]")}`$liftover $tmp_dir/input.bed $remap_chain /dev/stdout /dev/null 2> /dev/null`;
    unlink( "$tmp_dir/input.bed" );

    # Create a new VCF in the temp folder, with remapped loci
    my $vcf_fh = IO::File->new( $input_vcf ) or die "ERROR: Couldn't open --input-vcf: $input_vcf!\n";
    my $remap_vcf_fh = IO::File->new( "$tmp_dir/input.remap.vcf", "w" ) or die "ERROR: Couldn't open VCF: $tmp_dir/input.remap.vcf!\n";
    while( my $line = $vcf_fh->getline ) {
        # If the file uses Mac OS 9 newlines, quit with an error
        ( $line !~ m/\r$/ ) or die "ERROR: Your VCF uses CR line breaks, which we can't support. Please use LF or CRLF.\n";

        if( $line =~ m/^#/ ) {
            $remap_vcf_fh->print( $line ); # Write header lines unchanged
        }
        else {
            chomp( $line );
            my @cols = split( "\t", $line );
            my $locus = $cols[0] . ":" . $cols[1];
            if( defined $remap{$locus} ) {
                @cols[0,1] = split( ":", $remap{$locus} );
                $remap_vcf_fh->print( join( "\t", @cols ), "\n" );
            }
            else {
                warn "WARNING: Skipping variant at $locus; Unable to liftOver using $remap_chain\n";
            }
        }
    }
    $remap_vcf_fh->close;
    $vcf_fh->close;
    $input_vcf = "$tmp_dir/input.remap.vcf";
}

# Parse through each variant in the input VCF, and fix 'em up
my $vcf_in_fh = IO::File->new( $input_vcf ) or die "ERROR: Couldn't open input VCF file: $input_vcf!\n";
my ( $vcf_tumor_idx, $vcf_normal_idx, %vcf_header, @vcf_lines );
while( my $line = $vcf_in_fh->getline ) {

    # Skip all but a few INFO/FORMAT/FILTER descriptors in the header
    if( $line =~ m/^##/ ) {
        # Retain only the minimal INFO and FORMAT descriptor lines in the header
        if( $line =~ m/^##(INFO|FORMAT|FILTER)=<ID=([^,]+)/ ) {
            my ( $type, $tag ) = ( $1, $2 );
            if(( $type eq "INFO" and grep( /^$tag$/, @retain_info_tags )) or ( $type eq "FORMAT" and grep( /^$tag$/, @retain_format_tags )) or $type eq "FILTER" ) {
                $vcf_header{$type}{$tag} = $line;
            }
        }
        elsif( $line =~ m/^##(source|reference|assembly|phasing)=/ ) {
            my $tag = $1;
            $vcf_header{MAIN}{$tag} = $line;
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
        push( @vcf_lines, "#" . join( "\t", qw( CHROM POS ID REF ALT QUAL FILTER INFO FORMAT ), $new_tumor_id, $new_normal_id ) . "\n" );
        next;
    }

    # Set QUAL and FILTER to "." unless defined and non-empty
    $qual = "." unless( defined $qual and $qual ne "" );
    $filter = "." unless( defined $filter and $filter ne "" );

    # Parse out the data in the info column
    my %info = map {( m/=/ ? ( split( /=/, $_, 2 )) : ( $_, "" ))} split( /\;/, $info_line );

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

    # If we have either/both tumor/normal BAMs, then let's generate mpileup in VCF format
    if( $tumor_bam or $normal_bam ) {

        # Let's wipe clean values for DP:AD:ADF:ADR since we're about to override them
        $tum_info{DP} = $tum_info{AD} = $tum_info{ADF} = $tum_info{ADR} = ".";
        $nrm_info{DP} = $nrm_info{AD} = $nrm_info{ADF} = $nrm_info{ADR} = ".";

        # Generate mpileup and parse out only DP,AD,ADF,ADR for tumor/normal samples
        my @p_lines = `samtools mpileup --region $chrom:$pos-$pos --count-orphans --no-BAQ --min-MQ 1 --min-BQ 5 --ignore-RG --excl-flags UNMAP,SECONDARY,QCFAIL,DUP --VCF --uncompressed --output-tags DP,AD,ADF,ADR --ext-prob 20 --gap-frac 0.002 --tandem-qual 80 --min-ireads 1 --open-prob 30 --fasta-ref $ref_fasta $tumor_bam $normal_bam 2> /dev/null`;

        my ( $p_vcf_tumor_idx, $p_vcf_normal_idx ) = ( 0, 1 );
        foreach my $p_line ( @p_lines ) {

            # Retain only the minimal INFO and FORMAT descriptor lines in the header
            if( $p_line =~ m/^##/ ) {
                if( $p_line =~ m/^##(INFO|FORMAT)=<ID=([^,]+)/ ) {
                    my ( $type, $tag ) = ( $1, $2 );
                    if(( $type eq "INFO" and grep( /^$tag$/, @retain_info_tags )) or ( $type eq "FORMAT" and grep( /^$tag$/, @retain_format_tags ))) {
                        $vcf_header{$type}{$tag} = $p_line;
                    }
                }
                next;
            }

            # Parse the column headers to locate the tumor/normal genotype fields
            chomp( $p_line );
            my ( $p_chrom, $p_pos, undef, $p_ref, $p_alt, undef, undef, $p_info_line, $p_format_line, @p_rest ) = split( "\t", $p_line );
            if( $p_line =~ m/^#CHROM/ and $p_format_line and scalar( @p_rest ) > 0 ) {
                for( my $i = 0; $i <= $#p_rest; ++$i ) {
                    $p_vcf_tumor_idx = $i if( $p_rest[$i] eq $tumor_bam );
                    $p_vcf_normal_idx = $i if( $p_rest[$i] eq $normal_bam );
                }
                next;
            }

            # Parse out the data in the info column, and add it to the INFO fields from input
            %info = ( %info, map {( m/=/ ? ( split( /=/, $_, 2 )) : ( $_, "" ))} split( /\;/, $p_info_line ));

            # Parse out the mpileup REF/ALT alleles and match them to the input REF/ALT alleles
            my @p_alts = split( /,/, $p_alt );
            my %p_allele_idx = ( $alleles[0], 0 ); # The reference allele is always the first one
            for( my $i = 1; $i <= $#alleles; ++$i ) {
                foreach my $p_alt ( @p_alts ) {
                    # De-pad suffixed bps that are identical between ref/var alleles
                    my $p_ref_norm = $p_ref;
                    while( $p_ref_norm and $p_alt and substr( $p_ref_norm, -1, 1 ) eq substr( $p_alt, -1, 1 ) and $p_ref_norm ne $p_alt ) {
                        # Just in case the input also has one or more suffixed bps
                        last if( $p_ref_norm eq $alleles[0] and $p_alt eq $alleles[$i] );
                        ( $p_ref_norm, $p_alt ) = map{substr( $_, 0, -1 )} ( $p_ref_norm, $p_alt );
                    }
                    # Make sure that both REF and ALT alleles match, because deletions can vary
                    if( $p_ref_norm eq $alleles[0] and $p_alt eq $alleles[$i] ) {
                        $p_allele_idx{$p_alt} = $i;
                        next;
                    }
                }
            }

            # If genotype fields were generated for this variant, then parse out relevant data
            my ( $idx, %p_tum_info, %p_nrm_info );
            my @format_keys = split( /\:/, $p_format_line );
            if( defined $p_rest[$p_vcf_tumor_idx] ) {
                $idx = 0;
                %p_tum_info = map{( $format_keys[$idx++], $_ )} split( /\:/, $p_rest[$p_vcf_tumor_idx] );
                $tum_info{DP} = $p_tum_info{DP};
                my @depths = split( ",", $p_tum_info{AD} );
                $tum_info{AD} = join( ",", map{(defined $p_allele_idx{$_} ? $depths[$p_allele_idx{$_}] : 0)} @alleles );
                @depths = split( ",", $p_tum_info{ADF} );
                $tum_info{ADF} = join( ",", map{(defined $p_allele_idx{$_} ? $depths[$p_allele_idx{$_}] : 0)} @alleles );
                @depths = split( ",", $p_tum_info{ADR} );
                $tum_info{ADR} = join( ",", map{(defined $p_allele_idx{$_} ? $depths[$p_allele_idx{$_}] : 0)} @alleles );
            }
            if( defined $p_rest[$p_vcf_normal_idx] ) {
                $idx = 0;
                %p_nrm_info = map{( $format_keys[$idx++], $_ )} split( /\:/, $p_rest[$p_vcf_normal_idx] );
                $nrm_info{DP} = $p_nrm_info{DP};
                my @depths = split( ",", $p_nrm_info{AD} );
                $nrm_info{AD} = join( ",", map{(defined $p_allele_idx{$_} ? $depths[$p_allele_idx{$_}] : 0)} @alleles );
                @depths = split( ",", $p_nrm_info{ADF} );
                $nrm_info{ADF} = join( ",", map{(defined $p_allele_idx{$_} ? $depths[$p_allele_idx{$_}] : 0)} @alleles );
                @depths = split( ",", $p_nrm_info{ADR} );
                $nrm_info{ADR} = join( ",", map{(defined $p_allele_idx{$_} ? $depths[$p_allele_idx{$_}] : 0)} @alleles );
            }
        }
    }

    # Add more filter tags to the FILTER field, if --add-filters was specified
    if( $add_filters ) {
        my %tags = ();
        map{ $tags{$_} = 1 unless( $_ eq "PASS" or $_ eq "." )} split( /,|;/, $filter );
        $tags{LowTotalDepth} = 1 if(( $tum_info{DP} ne "." and $tum_info{DP} < $min_tum_depth ) or ( $nrm_info{DP} ne "." and $nrm_info{DP} < $min_nrm_depth ));
        my @tum_depths = split( /,/, $tum_info{AD} );
        my @nrm_depths = split( /,/, $nrm_info{AD} );
        my $tum_alt_depth = $tum_depths[$var_allele_idx];
        my $nrm_alt_depth = $nrm_depths[$var_allele_idx];
        $tags{LowTumorSupport} = 1 if( $tum_alt_depth ne "." and $tum_alt_depth < $min_tum_support );
        $tags{HighNormalSupport} = 1 if( $nrm_alt_depth ne "." and $nrm_alt_depth > $max_nrm_support );
        my $tags_to_add = join( ";", sort keys %tags );
        $filter = ( $tags_to_add ? $tags_to_add : $filter );
    }

    # Retain only the default or user-defined INFO fields, and add any additional fields if defined
    $info_line = join( ";", map{( $info{$_} eq "" ? "$_" : "$_=$info{$_}" )} grep{defined $info{$_}} @retain_info_tags );
    $info_line = ( $info_line ? ( $add_info ? join( ";", $add_info, $info_line ) : $info_line ) : $add_info );

    # Print a VCF file with default or user-defined genotype fields for tumor and normal
    $format_line = join( ":", @retain_format_tags );
    my $tum_fmt = join( ":", map{( defined $tum_info{$_} ? $tum_info{$_} : "." )} @retain_format_tags );
    my $nrm_fmt = join( ":", map{( defined $nrm_info{$_} ? $nrm_info{$_} : "." )} @retain_format_tags );
    push( @vcf_lines, join( "\t", $chrom, $pos, $ids, $ref, $alt, $qual, $filter, $info_line, $format_line, $tum_fmt, $nrm_fmt ) . "\n" );
}
$vcf_in_fh->close;

# Initialize the output VCF file with the mandatory first line of a VCF header
my $vcf_out_fh = IO::File->new( $output_vcf, ">" ) or die "ERROR: Couldn't open output VCF file: $output_vcf!\n";
$vcf_out_fh->print( "##fileformat=VCFv4.2\n" );

# Also add today's date just cuz we can
my $file_date = strftime( "%Y%m%d", localtime );
$vcf_out_fh->print( "##fileDate=$file_date\n" );

# Append all shortlisted MAIN/INFO/FORMAT/FILTER descriptors collected earlier
$vcf_out_fh->print( map{$vcf_header{MAIN}{$_}} sort keys %{$vcf_header{MAIN}} );
$vcf_out_fh->print( map{$vcf_header{INFO}{$_}} sort keys %{$vcf_header{INFO}} );
$vcf_out_fh->print( map{$vcf_header{FORMAT}{$_}} sort keys %{$vcf_header{FORMAT}} );
$vcf_out_fh->print( map{$vcf_header{FILTER}{$_}} sort keys %{$vcf_header{FILTER}} );

# Append any user-defined header lines
if( $add_header ) {
    $add_header =~ s/\s*$/\n/;
    $vcf_out_fh->print( $add_header );
}

# Append our custom FILTER tag descriptors, if --add-filters was specified
if( $add_filters ) {
    $vcf_out_fh->print( "##FILTER=<ID=$_,Description=\"" . $filter_tags{$_} . "\">\n" ) foreach ( sort keys %filter_tags );
}

# Append all the column header and variant lines accumulated earlier
$vcf_out_fh->print( @vcf_lines );
$vcf_out_fh->close;

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
    # Handle VCF lines by CaVEMan, where allele depths are in FAZ:FCZ:FGZ:FTZ:RAZ:RCZ:RGZ:RTZ
    elsif( !defined $fmt_info{AD} and scalar( grep{defined $fmt_info{$_}} qw/FAZ FCZ FGZ FTZ RAZ RCZ RGZ RTZ/ ) == 8 ) {
        # Create tags for forward+reverse strand reads, and use those to determine REF/ALT depths
        map{ $fmt_info{$_} = $fmt_info{'F'.$_} + $fmt_info{'R'.$_} } qw( AZ CZ GZ TZ );
        @depths = map{( defined $fmt_info{$_.'Z'} ? $fmt_info{$_.'Z'} : "" )} @alleles;
    }
    # Handle VCF lines from the Ion Torrent Suite where ALT depths are in AO and REF depths are in RO
    elsif( !defined $fmt_info{AD} and defined $fmt_info{AO} and defined $fmt_info{RO} ) {
        @depths = ( $fmt_info{RO}, map{( m/^\d+$/ ? $_ : "" )}split( /,/, $fmt_info{AO} ));
    }
    # Handle VCF lines from Delly where REF/ALT SV junction read counts are in RR/RV respectively
    elsif( !defined $fmt_info{AD} and defined $fmt_info{RR} and defined $fmt_info{RV} ) {
        # Reference allele depth and depths for any other ALT alleles must be left undefined
        @depths = map{""} @alleles;
        $depths[0] = $fmt_info{RR};
        $depths[$var_allele_idx] = $fmt_info{RV};
    }
    # Handle VCF lines from cgpPindel, where ALT depth and total depth are in PP:NP:PR:NR
    elsif( !defined $fmt_info{AD} and scalar( grep{defined $fmt_info{$_}} qw/PP NP PR NR/ ) == 4 ) {
        # Reference allele depth and depths for any other ALT alleles must be left undefined
        @depths = map{""} @alleles;
        $depths[$var_allele_idx] = $fmt_info{PP} + $fmt_info{NP};
        $fmt_info{DP} = $fmt_info{PR} + $fmt_info{NR};
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

 vcf2vcf.pl - Create a standardized VCF from a given VCF, fixing common problems in the process

=head1 SYNOPSIS

 perl vcf2vcf.pl --help
 perl vcf2vcf.pl --input-vcf weird.vcf --output-vcf not_weird.vcf --vcf-tumor-id WD4086 --vcf-normal-id NB4086

=head1 OPTIONS

 --input-vcf      Path to input file, a VCF or VCF-like format
 --output-vcf     Path to output file, a standardized VCF format
 --vcf-tumor-id   Tumor sample ID used in VCF's genotype column [TUMOR]
 --vcf-normal-id  Matched normal ID used in VCF's genotype column [NORMAL]
 --new-tumor-id   Tumor sample ID to use in the new VCF [--vcf-tumor-id]
 --new-normal-id  Matched normal ID to use in the new VCF [--vcf-normal-id]
 --tumor-bam      Path to tumor BAM, if provided, will add or override DP:AD:ADF:ADR in output VCF
 --normal-bam     Path to normal BAM, if provided, will add or override DP:AD:ADF:ADR in output VCF
 --ref-fasta      Reference FASTA file [~/.vep/homo_sapiens/91_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz]
 --add-header     VCF-style header lines to add to the output VCF; Use "\n" to separate lines []
 --add-info       Comma-delimited tag=value pairs to add as INFO fields in the output VCF []
 --retain-info    Comma-delimited names of INFO fields to retain in output VCF [SOMATIC,SS,I16,MQSB]
 --retain-format  Comma-delimited names of FORMATted genotype fields to retain in output VCF [GT,AD,DP]
 --remap-chain    Chain file to remap variants to a different assembly before standardizing VCF
 --add-filters    Use this to add some extra tags under FILTER [Default: 0]
 --help           Print a brief help message and quit
 --man            Print the detailed manual

=head1 DESCRIPTION

The VCF files that variant callers generate are rarely compliant with VCF specifications. This script fixes the most serious grievances, and creates a VCF with only important fields in INFO and FORMAT, like GT:AD:DP. You may optionally specify --add-filters, to use these allele depths and fractions to add more tags under FILTER.

=head2 Relevant links:

 Homepage: https://github.com/ckandoth/vcf2maf
 VCF format: http://samtools.github.io/hts-specs/

=head1 AUTHORS

 Cyriac Kandoth (ckandoth@gmail.com)

=head1 LICENSE

 Apache-2.0 | Apache License, Version 2.0 | https://www.apache.org/licenses/LICENSE-2.0

=cut
