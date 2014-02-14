#!/usr/bin/env perl

use strict;
use warnings;
use IO::File;
use File::Temp qw( tempdir );
use Getopt::Long qw( GetOptions );
use Pod::Usage qw( pod2usage );

# Set any default paths and constants
my $snpeff_cmd = "java -Xmx2g -jar ~/snpEff/snpEff.jar eff -config ~/snpEff/snpEff.config -noStats -hgvs GRCh37.74";

# Check for missing or crappy arguments
unless( @ARGV and $ARGV[0]=~m/^-/ ) {
    pod2usage( -verbose => 0, -message => "$0: Missing or invalid arguments!\n", -exitval => 2 );
}

# Parse options and print usage if there is a syntax error, or if usage was explicitly requested
my ( $man, $help ) = ( 0, 0 );
my ( $input_vcf, $output_maf, $tumor_id, $normal_id );
GetOptions(
    'help!' => \$help,
    'man!' => \$man,
    'input-vcf=s' => \$input_vcf,
    'output-maf=s' => \$output_maf,
    'tumor-id=s' => \$tumor_id,
    'normal-id=s' => \$normal_id,
    'snpeff-cmd=s' => \$snpeff_cmd
) or pod2usage( -verbose => 1, -input => \*DATA, -exitval => 2 );

pod2usage( -verbose => 1, -input => \*DATA, -exitval => 0 ) if( $help );
pod2usage( -verbose => 2, -input => \*DATA, -exitval => 0 ) if( $man );

# Make sure all our input files exist and are non-zero
( -s $input_vcf ) or die "Input VCF file missing or empty!\nPath: $input_vcf\n";

# Run snpEff to annotate the VCF variants to all possible isoforms, unless it was done earlier
my $snpeff_vcf = $input_vcf;
$snpeff_vcf =~ s/.vcf$//;
$snpeff_vcf .= ".anno.vcf";
if( -s $snpeff_vcf ) {
    print "Reusing this snpEff annotated VCF: $snpeff_vcf\n";
}
else {
    print "Running snpEff to annotate each variant to all possible isoforms...\n";
    print `$snpeff_cmd $input_vcf > $snpeff_vcf`;
}

# Use STDIN/STDOUT if input/output filenames are not defined
my ( $vcf_fh, $maf_fh ) = ( *STDIN, *STDOUT );
$vcf_fh = IO::File->new( $snpeff_vcf ) or die "Couldn't open file $snpeff_vcf! $!" if( $snpeff_vcf );
$maf_fh = IO::File->new( $output_maf, ">" ) or die "Couldn't open file $output_maf! $!" if( $output_maf );

# Predefine the columns of the MAF Header, and add 4 columns to the end for protein changes
# Ref: https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification
my @mafHeader = qw(
    Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome Start_Position End_Position Strand
    Variant_Classification Variant_Type Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2
    dbSNP_RS dbSNP_Val_Status Tumor_Sample_Barcode Matched_Norm_Sample_Barcode
    Match_Norm_Seq_Allele1 Match_Norm_Seq_Allele2 Tumor_Validation_Allele1 Tumor_Validation_Allele2
    Match_Norm_Validation_Allele1 Match_Norm_Validation_Allele2 Verification_Status Validation_Status
    Mutation_Status Sequencing_Phase Sequence_Source Validation_Method Score BAM_File Sequencer
    Tumor_Sample_UUID Matched_Norm_Sample_UUID Codon_Change Protein_Change Transcript_ID Exon_Number
);

# Parse through each variant in the VCF, pull out the snpEff effects, and choose one isoform per
# variant whose annotation will be used in the MAF
my ( $num_snvs, $num_indels, $num_svs ) = ( 0, 0, 0 );
$maf_fh->print( join( "\t", @mafHeader ), "\n" ); # Print MAF header
while( my $line = $vcf_fh->getline ) {

    # Skip all comment lines and the header cuz we know what we're doing!
    next if( $line =~ m/^#/ );

    # We don't need any of the data from the formatted columnms. Fetch everything else
    chomp( $line );
    my ( $chrom, $pos, $ids, $ref, $alt, $qual, $filter, $info_line ) = split( /\t/, $line );

    # Parse out the data in the info column, and store into a hash
    my %info = map {(m/=/ ? (split(/=/)) : ($_,1))} split( /\;/, $info_line );

    # Figure out the appropriate start/stop loci and var type/allele to report in the MAF
    my $start = my $stop = my $var = my $var_type = "";
    my ( $ref_length, $alt_length, $indel_size ) = ( length( $ref ), length( $alt ), 0 );
    if( $ref_length == $alt_length ) { # Handle SNP, DNP, TNP, or ONP
        ( $start, $stop ) = ( $pos, $pos + $alt_length - 1 );
        $var = $alt;
        $num_snvs++;
        my %np_type = qw( 1 SNP 2 DNP 3 TNP );
        $var_type = ( $alt_length > 3 ? "ONP" : $np_type{$alt_length} );
    }
    elsif( $ref_length < $alt_length ) { # Handle insertions
        $indel_size = length( $alt ) - $ref_length;
        ( $ref, $var ) = ( "-", substr( $alt, $ref_length, $indel_size ));
        ( $start, $stop ) = ( $pos, $pos + 1 );
        $num_indels++;
        $var_type = "INS";
    }
    elsif( $ref_length > $alt_length ) { # Handle deletions
        $indel_size = length( $ref ) - $alt_length;
        ( $ref, $var ) = ( substr( $ref, $alt_length, $indel_size ), "-" );
        ( $start, $stop ) = ( $pos + 1, $pos + $indel_size );
        $num_indels++;
        $var_type = "DEL";
    }
    else { # Poop-out on unknown variant types
        die "Unhandled variant type in VCF:\n$line\nPlease check/update this script!";
    }

    # Parse through each snpEff effect, which is stored in this format:
    # Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change | Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon  | GenotypeNum [ | ERRORS | WARNINGS ] )
    my $effLine = $info{EFF};
    my @effList = split( /,/, $effLine ) if( $effLine );

    # Load all effects & details into lists, skipping those with errors/warnings, and find the highest rank of their effect types, as defined in sub GetEffectPriority
    my $bestEffectRank;
    my @allEffects = ();
    my @allEffectDetails = ();
    foreach my $eff ( @effList ) {
        if( $eff =~ /^(\w+)\((.+)\)$/ ) {
            my $effect = $1;
            my @effectDetails = split( /\|/, $2 );
            my ( $warnings, $errors ) = @effectDetails[11,12];
            unless( $warnings or $errors ) {
                $bestEffectRank = GetEffectPriority( $effect ) if( !$bestEffectRank or $bestEffectRank > GetEffectPriority( $effect ));
                push( @allEffects, $effect );
                push( @allEffectDetails, \@effectDetails );
            }
        }
    }

    # Among effects with the bestEffectRank, find the highest rank of their biotypes, as defined in sub GetBiotypePriority
    my $bestBiotypeRank;
    for( my $i = 0; $i < scalar( @allEffectDetails ); $i++ ) {
        my $effect = $allEffects[$i];
        my @effectDetails = @{$allEffectDetails[$i]};
        my $biotype = $effectDetails[6];
        if( GetEffectPriority( $effect ) == $bestEffectRank ) {
            $bestBiotypeRank = GetBiotypePriority( $biotype ) if( !$bestBiotypeRank or $bestBiotypeRank > GetBiotypePriority( $biotype ));
        }
    }

    # Among effects with the bestEffectRank and bestBiotypeRank, choose the longest isoform as the one we annotate our variant to
    my $longest_aa_length = 0;
    my ( $bestEffect, @bestEffectDetails ) = ( "", () );
    for( my $i = 0; $i < scalar( @allEffectDetails ); $i++ ) {
        my $effect = $allEffects[$i];
        my @effectDetails = @{$allEffectDetails[$i]};
        my ( $aa_length, $biotype ) = @effectDetails[4,6];
        $aa_length = 0 unless( $aa_length );
        if( GetEffectPriority( $effect ) == $bestEffectRank and GetBiotypePriority( $biotype ) == $bestBiotypeRank and $longest_aa_length <= $aa_length ) {
            $longest_aa_length = $aa_length;
            $bestEffect = $effect;
            @bestEffectDetails = @effectDetails;
        }
    }

    # Construct the MAF columns and print to file
    my %mafLine = map{($_,'')} @mafHeader;
    $mafLine{'Hugo_Symbol'} = ( $bestEffectDetails[5] ? $bestEffectDetails[5] : 'Unknown' );
    $mafLine{'Entrez_Gene_Id'} = '0';
    $mafLine{'Center'} = '.';
    $mafLine{'NCBI_Build'} = 37;
    $mafLine{'Chromosome'} = $chrom;
    $mafLine{'Start_Position'} = $start;
    $mafLine{'End_Position'} = $stop;
    $mafLine{'Strand'} = '+';
    $mafLine{'Variant_Classification'} = GetVariantClassification($bestEffect, \@bestEffectDetails, $var_type);
    $mafLine{'Variant_Type'} = $var_type;
    $mafLine{'Reference_Allele'} = $ref;
    $mafLine{'Tumor_Seq_Allele1'} = $var;
    $mafLine{'Tumor_Seq_Allele2'} = $var;
    $mafLine{'dbSNP_RS'} = GetrsIDs( $ids );
    $mafLine{'Tumor_Sample_Barcode'} = $tumor_id;
    $mafLine{'Matched_Norm_Sample_Barcode'} = $normal_id;
    $mafLine{'Codon_Change'} = ( $bestEffectDetails[2] ? $bestEffectDetails[2] : "." );
    $mafLine{'Protein_Change'} = ( $bestEffectDetails[3] ? $bestEffectDetails[3] : "." );
    $mafLine{'Transcript_ID'} = ( $bestEffectDetails[8] ? $bestEffectDetails[8] : "." );
    $mafLine{'Exon_Number'} = ( $bestEffectDetails[9] ? $bestEffectDetails[9] : "." );

    foreach my $col (@mafHeader) {
        $maf_fh->print( "\t" ) if ( $col ne $mafHeader[0] );
        $maf_fh->print( $mafLine{$col} );
        $maf_fh->print( "\n" ) if ( $col eq $mafHeader[$#mafHeader] );
    }
}

# Prioritize snpEff effects, in the following generalized order:
#  1. Truncating
#  2. Missense
#  3. 5'UTR
#  4. 3'UTR
#  5. Upstream
#  6. Downstream
#  7. Synonymous
#  8. Non-protein
#  9. Intron - Conserved
# 10. Intron
# 11. Intergenic - Conserved
# 12. Intergenic/intragenic
sub GetEffectPriority {
    my ($effect) = @_;
    my %effectPriority = (
        'SPLICE_SITE_ACCEPTOR' => 1,
        'SPLICE_SITE_DONOR' => 1,
        'START_LOST' => 1,
        'EXON_DELETED' => 1,
        'FRAME_SHIFT' => 1,
        'STOP_GAINED' => 1,
        'STOP_LOST' => 1,
        'START_GAINED' => 2,
        'RARE_AMINO_ACID' => 2,
        'NON_SYNONYMOUS_START' => 2,
        'NON_SYNONYMOUS_CODING' => 2,
        'CODON_CHANGE' => 2,
        'CODON_INSERTION' => 2,
        'CODON_CHANGE_PLUS_CODON_INSERTION' => 2,
        'CODON_DELETION' => 2,
        'CODON_CHANGE_PLUS_CODON_DELETION' => 2,
        'UTR_5_PRIME' => 3,
        'UTR_5_DELETED' => 3,
        'UTR_3_PRIME' => 4,
        'UTR_3_DELETED' => 4,
        'UPSTREAM' => 5,
        'DOWNSTREAM' => 6,
        'SYNONYMOUS_START' => 7,
        'CDS' => 7,
        'SYNONYMOUS_CODING' => 7,
        'SYNONYMOUS_STOP' => 7,
        'GENE' => 8,
        'TRANSCRIPT' => 8,
        'EXON' => 8,
        'INTRON_CONSERVED' => 9,
        'INTRON' => 10,
        'INTERGENIC_CONSERVED' => 11,
        'INTRAGENIC' => 12,
        'INTERGENIC' => 12
    );
    $effectPriority{$effect} or die "ERROR: Unrecognized effect \"$effect\". Please update your hashes!";
    return $effectPriority{$effect};
}

# Prioritize the transcript biotypes that variants are annotated to
sub GetBiotypePriority {
    my ( $biotype ) = @_;
    my %biotypePriority = (
        'protein_coding' => 1,
        'IG_C_gene' => 2,
        'IG_D_gene' => 2,
        'IG_J_gene' => 2,
        'IG_V_gene' => 2,
        'TR_C_gene' => 2,
        'TR_D_gene' => 2,
        'TR_J_gene' => 2,
        'TR_V_gene' => 2,
        'miRNA' => 3,
        'snRNA' => 3,
        'snoRNA' => 3,
        'rRNA' => 3,
        'lincRNA' => 3,
        'Mt_tRNA' => 4,
        'Mt_rRNA' => 4,
        'antisense' => 5,
        'sense_intronic' => 5,
        'sense_overlapping' => 5,
        '3prime_overlapping_ncrna' => 5,
        'misc_RNA' => 5,
        'processed_transcript' => 6,
        'TEC' => 6,
        'retained_intron' => 7,
        'nonsense_mediated_decay' => 7,
        'non_stop_decay' => 7,
        'pseudogene' => 8,
        'processed_pseudogene' => 8,
        'polymorphic_pseudogene' => 8,
        'transcribed_processed_pseudogene' => 8,
        'transcribed_unprocessed_pseudogene' => 8,
        'unitary_pseudogene' => 8,
        'unprocessed_pseudogene' => 8,
        'IG_C_pseudogene' => 8,
        'IG_J_pseudogene' => 8,
        'IG_V_pseudogene' => 8,
        'TR_J_pseudogene' => 8,
        'TR_V_pseudogene' => 8,
        '' => 9
    );
    $biotypePriority{$biotype} or die "ERROR: Unrecognized biotype \"$biotype\". Please update your hashes!";
    return $biotypePriority{$biotype};
}

# Converts snpEff effect types to MAF variant classifications
sub GetVariantClassification {
    my ( $effect, $ref_bestEffectDetails, $var_type ) = @_;
    return "5'UTR" if( $effect =~ /UTR_5/ );
    return "3'UTR" if( $effect =~ /UTR_3/ );
    return "IGR" if( $effect =~ /INTERGENIC/ );
    return "Intron" if ( $effect =~ /INTRON/ );
    return "Splice_Site" if( $effect =~ /^SPLICE_SITE/ );
    return "5'Flank" if( $effect eq 'UPSTREAM' );
    return "3'Flank" if ( $effect eq 'DOWNSTREAM' );
    return "Missense_Mutation" if( $effect eq 'NON_SYNONYMOUS_CODING' or $effect eq 'CODON_CHANGE' or $effect eq 'RARE_AMINO_ACID' );
    return "Translation_Start_Site" if( $effect eq 'START_LOST' );
    return "De_novo_Start_InFrame" if( $effect eq 'START_GAINED' );
    return "Nonsense_Mutation" if( $effect eq 'STOP_GAINED' );
    return "Nonstop_Mutation" if( $effect eq 'STOP_LOST' );
    return "Frame_Shift_Ins" if( $effect eq 'FRAME_SHIFT' and $var_type eq 'INS' );
    return "Frame_Shift_Del" if ( $effect eq 'FRAME_SHIFT' and $var_type eq 'DEL' );
    return "In_Frame_Ins" if( $effect eq 'CODON_INSERTION' or $effect eq 'CODON_CHANGE_PLUS_CODON_INSERTION' );
    return "In_Frame_Del" if( $effect eq 'CODON_DELETION' or $effect eq 'CODON_CHANGE_PLUS_CODON_DELETION' );
    return "Silent" if( $effect =~ /^SYNONYMOUS_/ );
    return "IGR" if( $effect eq 'INTRAGENIC' );

    # All non-coding RNA genes are grouped into one classification
    return "RNA" if( $effect eq 'EXON' and ${$ref_bestEffectDetails}[6]=~m/RNA/ and ${$ref_bestEffectDetails}[7]=~m/^NON_CODING$/ );

    # Annotate everything else simply as a targeted region
    return "Targeted_Region";
}

# If the VCF variant ID contains multiple entries, return the ones that look like dbSNP rsIDs
sub GetrsIDs {

    my ( $id_line ) = @_;
    return '' if( !$id_line or $id_line eq '.' or $id_line=~m/^\s*$/ ); # Handle null data
    my @rs_ids = grep( /^rs\d+$/, split( /;/, $id_line )); # Pull out rsIDs into a list

    # If no rsIDs were found, return null. Else return a semicolon delimited list
    return '' unless( scalar( @rs_ids ) > 0 );    
    return join( ";", @rs_ids );
}

__DATA__

=head1 NAME

 vcf2maf.pl - Convert VCF to MAF. Use snpEff to annotate calls to all possible isoforms and choose one to report.

=head1 SYNOPSIS

 perl vcf2maf.pl \
   --input-vcf WD1309_vs_NB1308.vcf --output-maf WD1309_vs_NB1308.maf \
   --tumor-id WD1309 --normal-id NB1308 \
   --snpeff-cmd "java -Xmx2g -jar ~/snpEff/snpEff.jar eff -config ~/snpEff/snpEff.config -noStats -hgvs GRCh37.74"

=head1 OPTIONS

=over 8

=item B<--tumor-id>

 The tumor sample ID to report in the 16th column of the MAF

=item B<--normal-id>

 The normal sample ID to report in the 17th column of the MAF

=item B<--snpeff-cmd>

 Command to run snpEff annotation on a VCF file

=item B<--help>

 Prints a brief help message and quits

=item B<--man>

 Prints the detailed manual

=back

=head1 DESCRIPTION

 To convert a VCF into a MAF, each variant must be annotated to only one of all possible gene transcripts/isoforms that it might affect. This selection of a single affected transcript/isoform per variant, is often subjective. For now, we try to follow best-practices, but over time, the selection process will be made smarter and more configurable.

 This script needs snpEff (snpeff.sourceforge.net), a variant annotator that can quickly map each variant to all possible transcripts in a database. snpEff is downloadable as a java archive, so make sure you also have Java installed.

=head1 AUTHORS

 Cyriac Kandoth (ckandoth@gmail.com)
 William Lee (leew1@cbio.mskcc.org)

=head1 LICENSE

 LGPLv3, Memorial Sloan-Kettering Cancer Center, New York, NY 10065, USA

=cut
