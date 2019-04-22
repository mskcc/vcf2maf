#!/usr/bin/env perl

# maf2maf - Reannotate the effects of variants in a MAF by running maf2vcf followed by vcf2maf

use strict;
use warnings;
use IO::File;
use Getopt::Long qw( GetOptions );
use Pod::Usage qw( pod2usage );
use File::Temp qw( tempdir );
use File::Copy qw( move );
use File::Path qw( mkpath );
use Config;

# Set any default paths and constants
my ( $tum_depth_col, $tum_rad_col, $tum_vad_col ) = qw( t_depth t_ref_count t_alt_count );
my ( $nrm_depth_col, $nrm_rad_col, $nrm_vad_col ) = qw( n_depth n_ref_count n_alt_count );
my ( $vep_path, $vep_data, $vep_forks, $buffer_size, $any_allele ) = ( "$ENV{HOME}/vep", "$ENV{HOME}/.vep", 4, 5000, 0 );
my ( $ref_fasta, $filter_vcf ) = ( "$ENV{HOME}/.vep/homo_sapiens/91_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz", "$ENV{HOME}/.vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz" );
my ( $species, $ncbi_build, $cache_version, $maf_center, $max_filter_ac ) = ( "homo_sapiens", "GRCh37", "", ".", 10 );
my $perl_bin = $Config{perlpath};

# Columns that can be safely borrowed from the input MAF
my $retain_cols = "Center,Verification_Status,Validation_Status,Mutation_Status,Sequencing_Phase" .
    ",Sequence_Source,Validation_Method,Score,BAM_file,Sequencer,Tumor_Sample_UUID" .
    ",Matched_Norm_Sample_UUID";

# Columns that should never be overridden since they are results of re-annotation
my %force_new_cols = map{ my $c = lc; ( $c, 1 )} qw( Hugo_Symbol Entrez_Gene_Id NCBI_Build
    Chromosome Start_Position End_Position Strand Variant_Classification Variant_Type
    Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 dbSNP_RS Tumor_Sample_Barcode
    Matched_Norm_Sample_Barcode Match_Norm_Seq_Allele1 Match_Norm_Seq_Allele2
    Tumor_Validation_Allele1 Tumor_Validation_Allele2 Match_Norm_Validation_Allele1
    Match_Norm_Validation_Allele2 HGVSc HGVSp HGVSp_Short Transcript_ID Exon_Number t_depth
    t_ref_count t_alt_count n_depth n_ref_count n_alt_count all_effects Allele Gene Feature
    Feature_type Consequence cDNA_position CDS_position Protein_position Amino_acids Codons
    Existing_variation ALLELE_NUM DISTANCE STRAND_VEP SYMBOL SYMBOL_SOURCE HGNC_ID BIOTYPE CANONICAL
    CCDS ENSP SWISSPROT TREMBL UNIPARC RefSeq SIFT PolyPhen EXON INTRON DOMAINS AF AFR_AF
    AMR_AF ASN_AF EAS_AF EUR_AF SAS_AF AA_AF EA_AF CLIN_SIG SOMATIC PUBMED MOTIF_NAME
    MOTIF_POS HIGH_INF_POS MOTIF_SCORE_CHANGE IMPACT PICK VARIANT_CLASS TSL HGVS_OFFSET PHENO
    MINIMISED ExAC_AF ExAC_AF_AFR ExAC_AF_AMR ExAC_AF_EAS ExAC_AF_FIN ExAC_AF_NFE ExAC_AF_OTH
    ExAC_AF_SAS GENE_PHENO FILTER flanking_bps variant_id variant_qual ExAC_AF_Adj ExAC_AC_AN_Adj
    ExAC_AC_AN ExAC_AC_AN_AFR ExAC_AC_AN_AMR ExAC_AC_AN_EAS ExAC_AC_AN_FIN ExAC_AC_AN_NFE
    ExAC_AC_AN_OTH ExAC_AC_AN_SAS ExAC_FILTER gnomAD_AF gnomAD_AFR_AF gnomAD_AMR_AF gnomAD_ASJ_AF
    gnomAD_EAS_AF gnomAD_FIN_AF gnomAD_NFE_AF gnomAD_OTH_AF gnomAD_SAS_AF );

# Check for missing or crappy arguments
unless( @ARGV and $ARGV[0]=~m/^-/ ) {
    pod2usage( -verbose => 0, -message => "$0: Missing or invalid arguments!\n", -exitval => 2 );
}

# Parse options and print usage syntax on a syntax error, or if help was explicitly requested
my ( $man, $help ) = ( 0, 0 );
my ( $input_maf, $output_maf, $tmp_dir, $custom_enst_file );
GetOptions(
    'help!' => \$help,
    'man!' => \$man,
    'input-maf=s' => \$input_maf,
    'output-maf=s' => \$output_maf,
    'tmp-dir=s' => \$tmp_dir,
    'tum-depth-col=s' => \$tum_depth_col,
    'tum-rad-col=s' => \$tum_rad_col,
    'tum-vad-col=s' => \$tum_vad_col,
    'nrm-depth-col=s' => \$nrm_depth_col,
    'nrm-rad-col=s' => \$nrm_rad_col,
    'nrm-vad-col=s' => \$nrm_vad_col,
    'retain-cols=s' => \$retain_cols,
    'custom-enst=s' => \$custom_enst_file,
    'vep-path=s' => \$vep_path,
    'vep-data=s' => \$vep_data,
    'vep-forks=s' => \$vep_forks,
    'buffer-size=i' => \$buffer_size,
    'any-allele!' => \$any_allele,
    'species=s' => \$species,
    'ncbi-build=s' => \$ncbi_build,
    'cache-version=s' => \$cache_version,
    'ref-fasta=s' => \$ref_fasta,
    'filter-vcf=s' => \$filter_vcf,
    'max-filter-ac=i' => \$max_filter_ac
) or pod2usage( -verbose => 1, -input => \*DATA, -exitval => 2 );
pod2usage( -verbose => 1, -input => \*DATA, -exitval => 0 ) if( $help );
pod2usage( -verbose => 2, -input => \*DATA, -exitval => 0 ) if( $man );

# Locate the maf2vcf and vcf2maf scripts that should be next to this script
my ( $script_dir ) = $0 =~ m/^(.*)\/maf2maf/;
$script_dir = "." unless( $script_dir );
my ( $maf2vcf_path, $vcf2maf_path ) = ( "$script_dir/maf2vcf.pl", "$script_dir/vcf2maf.pl" );
( -s $maf2vcf_path ) or die "ERROR: Couldn't locate maf2vcf.pl! Must be beside maf2maf.pl\n";
( -s $vcf2maf_path ) or die "ERROR: Couldn't locate vcf2maf.pl! Must be beside maf2maf.pl\n";

# Check if required arguments are missing or problematic
( defined $input_maf ) or die "ERROR: --input-maf must be defined!\n";
( -s $input_maf ) or die "ERROR: Provided --input-maf is missing or empty: $input_maf\n";
( -s $ref_fasta ) or die "ERROR: Provided --ref-fasta is missing or empty: $ref_fasta\n";
( $input_maf !~ m/\.(gz|bz2|bcf)$/ ) or die "ERROR: Unfortunately, --input-maf cannot be in a compressed format\n";

# Create a temporary directory for our intermediate files, unless the user wants to use their own
if( $tmp_dir ) {
    mkpath( $tmp_dir );
}
else {
    $tmp_dir = tempdir( CLEANUP => 1 );
}

# Construct a maf2vcf command and run it
my $maf2vcf_cmd = "$perl_bin $maf2vcf_path --input-maf $input_maf --output-dir $tmp_dir " .
    "--ref-fasta $ref_fasta --tum-depth-col $tum_depth_col --tum-rad-col $tum_rad_col " .
    "--tum-vad-col $tum_vad_col --nrm-depth-col $nrm_depth_col --nrm-rad-col $nrm_rad_col ".
    "--nrm-vad-col $nrm_vad_col --per-tn-vcfs";
system( $maf2vcf_cmd ) == 0 or die "\nERROR: Failed to run maf2vcf! Command: $maf2vcf_cmd\n";

my $vcf_file = "$tmp_dir/" . substr( $input_maf, rindex( $input_maf, '/' ) + 1 );
$vcf_file =~ s/(\.)?(maf|tsv|txt)?$/.vcf/;
my $vep_anno = $vcf_file;
$vep_anno =~ s/\.vcf$/.vep.vcf/;

# Skip running VEP if a VEP-annotated VCF already exists
if( -s $vep_anno ) {
    warn "WARNING: Annotated VCF already exists ($vep_anno). Skipping re-annotation.\n";
}
else {
    warn "STATUS: Running VEP and writing to: $vep_anno\n";
    # Make sure we can find the VEP script
    my $vep_script = ( -s "$vep_path/vep" ? "$vep_path/vep" : "$vep_path/variant_effect_predictor.pl" );
    ( -s $vep_script ) or die "ERROR: Cannot find VEP script in path: $vep_path\n";

    # Contruct VEP command using some default options and run it
    my $vep_cmd = "$perl_bin $vep_script --species $species --assembly $ncbi_build --offline --no_progress --no_stats --buffer_size $buffer_size --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --dir $vep_data --fasta $ref_fasta --format vcf --input_file $vcf_file --output_file $vep_anno";
    # VEP barks if --fork is set to 1. So don't use this argument unless it's >1
    $vep_cmd .= " --fork $vep_forks" if( $vep_forks > 1 );
    # Require allele match for co-located variants unless user-rejected or we're using a newer VEP
    $vep_cmd .= " --check_allele" unless( $any_allele or $vep_script =~ m/vep$/ );
    # Add --cache-version only if the user specifically asked for a version
    $vep_cmd .= " --cache_version $cache_version" if( $cache_version );
    # Add options that only work on human variants
    if( $species eq "homo_sapiens" ) {
        # Slight change in these arguments if using the newer VEP
        $vep_cmd .= " --polyphen b " . ( $vep_script =~ m/vep$/ ? "--af --af_1kg --af_esp --af_gnomad" : "--gmaf --maf_1kg --maf_esp" );
    }
    # Add options that work for most species, except a few
    $vep_cmd .= " --regulatory" unless( $species eq "canis_familiaris" );

    # Make sure it ran without error codes
    system( $vep_cmd ) == 0 or die "\nERROR: Failed to run the VEP annotator! Command: $vep_cmd\n";
    ( -s $vep_anno ) or warn "WARNING: VEP-annotated VCF file is missing or empty: $vep_anno\n";
}

# Load the tumor-normal pairs from the TSV created by maf2vcf
my $tsv_file = $vcf_file;
$tsv_file =~ s/(\.vcf)?$/.pairs.tsv/;

# Store the VEP annotated VCF header so we can duplicate it for per-TN VCFs
my $vep_vcf_header = `grep ^## $vep_anno`;

# Split the multi-sample VEP annotated VCF into per-TN VCFs
my ( %tn_pair, %t_col_idx, %n_col_idx, %tn_vep );
my $vep_fh = IO::File->new( $vep_anno ) or die "ERROR: Couldn't open file: $vep_anno\n";
while( my $line = $vep_fh->getline ) {

    # Skip comment lines, but parse everything else including the column headers
    next if( $line =~ m/^##/ );
    my @cols = map{s/^\s+|\s+$|\r|\n//g; $_} split( /\t/, $line );

    # Parse the header line to map column names to their indexes
    if( $line =~ m/^#CHROM/ ) {

        # Initialize VCF header and fill up %tn_pair for each tumor-normal pair
        foreach ( `grep -Ev ^# $tsv_file` ){
            chomp;
            my @ids = split( "\t", $_ );
            $t_col_idx{ $ids[ 0 ] } = 1;
            $n_col_idx{ $ids[ 1 ] } = 1;
            # If the same tumor is paired with different normals, treat them as separate TN-pairs
            $tn_pair{ $ids[ 0 ] }{ $ids[ 1 ] } = 1;
            my $tn_vcf_file = "$tmp_dir/$ids[0]\_vs_$ids[1].vep.vcf";
            $tn_vep{ $tn_vcf_file } = $vep_vcf_header . "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$ids[0]\t$ids[1]\n";
        }
        # Save VCF column indexes of tumors and normals in %t_col_idx and %n_col_idx, respectively.
        foreach my $idx ( 9..$#cols ){
            my $id = $cols[ $idx ];
            $t_col_idx{ $id } = $idx if( exists $t_col_idx{ $id } );
            $n_col_idx{ $id } = $idx if( exists $n_col_idx{ $id } );
        }
    }
    # For all other lines containing variants, cache it for the appropriate per-TN VCF
    else {
        my $GT_idx;
        my @format_keys = split( /\:/, $cols[8] );
        map{ $GT_idx = $_ if( $format_keys[ $_ ] eq "GT" ) } ( 0..$#format_keys );
        # Look for non-null normal genotypes
        my @n_cols;
        foreach my $n_id ( keys %n_col_idx ){
            my $n_idx = $n_col_idx{ $n_id };
            next if( $n_idx < 9 );
            my @n_info = split( /\:/, $cols[ $n_idx ] );
            ( $n_info[ $GT_idx ] eq './.' ) or push @n_cols, $n_id;
        }
        
        foreach my $t_id ( keys %t_col_idx ){
            my $t_idx = $t_col_idx{ $t_id };
            next if( $t_idx < 9 );
            my @t_info = split( /\:/, $cols[ $t_idx ] );
            # Skip variants for TN-pairs where the tumor genotype is null
            next if ( $t_info[ $GT_idx ] eq './.' );

            # Otherwise write it to the appropriate per-TN VCF file
            foreach my $n_id ( @n_cols ){
                my $n_idx = $n_col_idx{ $n_id };
                my $tn_vcf_file = "$tmp_dir/$t_id\_vs_$n_id.vep.vcf";
                $tn_vep{ $tn_vcf_file } .= join( "\t", @cols[0..8,$t_idx,$n_idx] ) . "\n" if( exists $tn_pair{ $t_id }{ $n_id } );
            }
        }
    }
}

# Write the cached contents of per-TN annotated VCFs into files
foreach my $tn_vcf_file ( keys %tn_vep ) {
    my $tn_vep_fh = IO::File->new( $tn_vcf_file, ">" ) or die "ERROR: Couldn't open file: $tn_vcf_file\n";
    $tn_vep_fh->print( $tn_vep{$tn_vcf_file} );
    $tn_vep_fh->close;
}

# For each VCF generated by maf2vcf above, contruct a vcf2maf command and run it
my @vcfs = grep{ !m/.vep.vcf$/ and !m/$vcf_file/ } glob( "$tmp_dir/*.vcf" ); # Avoid reannotating annotated VCFs
foreach my $tn_vcf ( @vcfs ) {
    my ( $tumor_id, $normal_id ) = $tn_vcf=~m/^.*\/(.*)_vs_(.*)\.vcf/;
    my $tn_maf = $tn_vcf;
    $tn_maf =~ s/.vcf$/.vep.maf/;
    my $vcf2maf_cmd = "$perl_bin $vcf2maf_path --input-vcf $tn_vcf --output-maf $tn_maf" .
        " --tumor-id $tumor_id --normal-id $normal_id --vep-path $vep_path --vep-data $vep_data" .
        " --vep-forks $vep_forks --ref-fasta $ref_fasta --ncbi-build $ncbi_build --species $species" .
        " --filter-vcf $filter_vcf --max-filter-ac $max_filter_ac";
    $vcf2maf_cmd .= " --custom-enst $custom_enst_file" if( $custom_enst_file );
    system( $vcf2maf_cmd ) == 0 or die "\nERROR: Failed to run vcf2maf! Command: $vcf2maf_cmd\n";
}

# Fetch the column header from one of the resulting MAFs
my @mafs = glob( "$tmp_dir/*.vep.maf" );
my $maf_header = `grep ^Hugo_Symbol $mafs[0] | head -n1`;
chomp( $maf_header );

# If user wants to retain some columns from the input MAF, fetch those and override
my %input_maf_data = ();
if( $retain_cols ) {

    # Parse the input MAF and fetch the data for columns that we need to retain/override
    my $input_maf_fh = IO::File->new( $input_maf ) or die "ERROR: Couldn't open --input-maf\n";
    my %input_maf_col_idx = (); # Hash to map column names to column indexes
    while( my $line = $input_maf_fh->getline ) {

        next if( $line =~ m/^#/ ); # Skip comments

        # Do a thorough removal of carriage returns, line feeds, prefixed/suffixed whitespace
        my @cols = map{s/^\s+|\s+$|\r|\n//g; $_} split( /\t/, $line );

        # Parse the header line to map column names to their indexes
        if( $line =~ m/^(Hugo_Symbol|Chromosome|Tumor_Sample_Barcode)/i ) {
            my $idx = 0;
            map{ my $c = lc; $input_maf_col_idx{$c} = $idx; ++$idx } @cols;

            # Check for columns to retain that are not in old MAF, or columns we shouldn't override in new MAF
            foreach my $c ( split( ",", $retain_cols )) {
                my $c_lc = lc( $c );
                if( !defined $input_maf_col_idx{$c_lc} ) {
                    warn "WARNING: Column '$c' not found in old MAF.\n";
                }
                elsif( $force_new_cols{$c_lc} ) {
                    warn "WARNING: Column '$c' cannot be overridden in new MAF.\n";
                }
            }
        }
        else {
            # Extract minimal variant info, and figure out which of the tumor alleles is non-REF
            my ( $chr, $pos, $ref, $al1, $al2, $sample_id ) = map{ my $c = lc; ( defined $input_maf_col_idx{$c} ? $cols[$input_maf_col_idx{$c}] : "" ) } qw( Chromosome Start_Position Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 Tumor_Sample_Barcode );
            # Prefer Tumor_Seq_Allele2 over Tumor_Seq_Allele1 if both are non-REF
            my $var = (( defined $al2 and $al2 and $al2 ne $ref ) ? $al2 : $al1 );

            # Remove any prefixed reference bps from alleles, using "-" for simple indels
            while( $ref and $var and substr( $ref, 0, 1 ) eq substr( $var, 0, 1 ) and $ref ne $var ) {
                ( $ref, $var ) = map{$_ = substr( $_, 1 ); ( $_ ? $_ : "-" )} ( $ref, $var );
                ++$pos unless( $ref eq "-" );
            }

            # Create a key for this variant using Chromosome:Start_Position:Tumor_Sample_Barcode:Reference_Allele:Variant_Allele
            my $key = join( ":", $chr, $pos, $sample_id, $ref, $var );
            %{$input_maf_data{$key}} = ();

            # Store values for this variant into a hash, adding column names to the key
            foreach my $c ( map{lc} split( ",", $retain_cols )) {
                if( defined $input_maf_col_idx{$c} and defined $cols[$input_maf_col_idx{$c}] ) {
                    $input_maf_data{$key}{$c} = $cols[$input_maf_col_idx{$c}];
                }
            }
        }
    }
    $input_maf_fh->close;

    # Add additional column headers for the output MAF, if any
    my %maf_cols = map{ my $c = lc; ( $c, 1 )} split( /\t/, $maf_header );
    my @addl_maf_cols = grep{ my $c = lc; !$maf_cols{$c} } split( ",", $retain_cols );
    map{ $maf_header .= "\t$_" } @addl_maf_cols;

    # Retain/override data in each of the per-TN-pair MAFs
    foreach my $tn_maf ( @mafs ) {
        my $tn_maf_fh = IO::File->new( $tn_maf ) or die "ERROR: Couldn't open file: $tn_maf\n";
        my %output_maf_col_idx = (); # Hash to map column names to column indexes
        my $tmp_tn_maf_fh = IO::File->new( "$tn_maf.tmp", ">" ) or die "ERROR: Couldn't open file: $tn_maf.tmp\n";
        while( my $line = $tn_maf_fh->getline ) {

            # Do a thorough removal of carriage returns, line feeds, prefixed/suffixed whitespace
            my @cols = map{ s/^\s+|\s+$|\r|\n//g; $_ } split( /\t/, $line );

            # Copy comment lines to the new MAF unchanged
            if( $line =~ m/^#/ ) {
                $tmp_tn_maf_fh->print( $line );
            }
            # Print the MAF header prepared earlier, but also create a hash with column indexes
            elsif( $line =~ m/^Hugo_Symbol/ ) {
                my $idx = 0;
                map{ my $c = lc; $output_maf_col_idx{$c} = $idx; ++$idx } ( @cols, @addl_maf_cols );
                $tmp_tn_maf_fh->print( "$maf_header\n" );
            }
            # For all other lines, insert the data collected from the original input MAF
            else {
                my $key = join( ":", map{ my $c = lc; $cols[$output_maf_col_idx{$c}] } qw( Chromosome Start_Position Tumor_Sample_Barcode Reference_Allele Tumor_Seq_Allele2 ));
                foreach my $c ( map{lc} split( /\t/, $maf_header )){
                    if( !$force_new_cols{$c} ) {
                        $cols[$output_maf_col_idx{$c}] = $input_maf_data{$key}{$c} if( defined $input_maf_data{$key}{$c} );
                    }
                }
                $tmp_tn_maf_fh->print( join( "\t", @cols ) . "\n" );
            }
        }
        $tmp_tn_maf_fh->close;
        $tn_maf_fh->close;

        # Overwrite the old MAF with the new one containing data from the original input MAF
        move( "$tn_maf.tmp", $tn_maf );
    }
}

# Concatenate the per-TN-pair MAFs into the user-specified final MAF
# Default to printing to screen if an output MAF was not defined
my $maf_fh = *STDOUT;
if( $output_maf ) {
    $maf_fh = IO::File->new( $output_maf, ">" ) or die "ERROR: Couldn't open --output-maf\n";
}
$maf_fh->print( "#version 2.4\n$maf_header\n" );
foreach my $tn_maf ( @mafs ) {
    my @maf_lines = `grep -Ev "^#|^Hugo_Symbol" $tn_maf`;
    $maf_fh->print( @maf_lines );
}
$maf_fh->close;

__DATA__

=head1 NAME

 maf2maf.pl - Reannotate the effects of variants in a MAF by running maf2vcf followed by vcf2maf

=head1 SYNOPSIS

 perl maf2maf.pl --help
 perl maf2maf.pl --input-maf test.maf --output-maf test.vep.maf

=head1 OPTIONS

 --input-maf      Path to input file in MAF format
 --output-maf     Path to output MAF file [Default: STDOUT]
 --tmp-dir        Folder to retain intermediate VCFs/MAFs after runtime [Default: usually /tmp]
 --tum-depth-col  Name of MAF column for read depth in tumor BAM [t_depth]
 --tum-rad-col    Name of MAF column for reference allele depth in tumor BAM [t_ref_count]
 --tum-vad-col    Name of MAF column for variant allele depth in tumor BAM [t_alt_count]
 --nrm-depth-col  Name of MAF column for read depth in normal BAM [n_depth]
 --nrm-rad-col    Name of MAF column for reference allele depth in normal BAM [n_ref_count]
 --nrm-vad-col    Name of MAF column for variant allele depth in normal BAM [n_alt_count]
 --retain-cols    Comma-delimited list of columns to retain from the input MAF [Center,Verification_Status,Validation_Status,Mutation_Status,Sequencing_Phase,Sequence_Source,Validation_Method,Score,BAM_file,Sequencer,Tumor_Sample_UUID,Matched_Norm_Sample_UUID]
 --custom-enst    List of custom ENST IDs that override canonical selection
 --vep-path       Folder containing the vep script [~/vep]
 --vep-data       VEP's base cache/plugin directory [~/.vep]
 --vep-forks      Number of forked processes to use when running VEP [4]
 --buffer-size    Number of variants VEP loads at a time; Reduce this for low memory systems [5000]
 --any-allele     When reporting co-located variants, allow mismatched variant alleles too
 --filter-vcf     A VCF for FILTER tag common_variant. Set to 0 to disable [~/.vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz]
 --max-filter-ac  Use tag common_variant if the filter-vcf reports a subpopulation AC higher than this [10]
 --species        Ensembl-friendly name of species (e.g. mus_musculus for mouse) [homo_sapiens]
 --ncbi-build     NCBI reference assembly of variants in MAF (e.g. GRCm38 for mouse) [GRCh37]
 --cache-version  Version of offline cache to use with VEP (e.g. 75, 84, 91) [Default: Installed version]
 --ref-fasta      Reference FASTA file [~/.vep/homo_sapiens/91_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz]
 --help           Print a brief help message and quit
 --man            Print the detailed manual

=head1 DESCRIPTION

This script runs a given MAF through maf2vcf to generate per-TN-pair VCFs in a temporary folder, and then runs vcf2maf on each VCF to reannotate variant effects and create a new combined MAF

=head2 Relevant links:

 Homepage: https://github.com/ckandoth/vcf2maf
 VCF format: http://samtools.github.io/hts-specs/
 MAF format: https://wiki.nci.nih.gov/x/eJaPAQ
 VEP: http://ensembl.org/info/docs/tools/vep/index.html
 VEP annotated VCF format: http://ensembl.org/info/docs/tools/vep/vep_formats.html#vcfout

=head1 AUTHORS

 Cyriac Kandoth (ckandoth@gmail.com)
 Qingguo Wang (josephw10000@gmail.com)

=head1 LICENSE

 Apache-2.0 | Apache License, Version 2.0 | https://www.apache.org/licenses/LICENSE-2.0

=cut
