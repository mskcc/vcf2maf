#!/usr/bin/env perl

use strict;
use warnings;
use IO::File;
use Getopt::Long qw( GetOptions );
use Pod::Usage qw( pod2usage );

# Set any default paths and constants
my ( $vep_path, $vep_data ) = ( "~/vep", "~/.vep" );
my ( $snpeff_path, $snpeff_data ) = ( "~/snpEff", "~/snpEff/data" );
my ( $tumor_id, $normal_id ) = ( "TUMOR", "NORMAL" );
my ( $ncbi_build, $maf_center ) = ( 37, "mskcc.org" );

# Check for missing or crappy arguments
unless( @ARGV and $ARGV[0]=~m/^-/ ) {
    pod2usage( -verbose => 0, -message => "$0: Missing or invalid arguments!\n", -exitval => 2 );
}

# Parse options and print usage if there is a syntax error, or if usage was explicitly requested
my ( $man, $help, $use_snpeff ) = ( 0, 0, 0 );
my ( $input_vcf, $vep_anno, $snpeff_anno, $output_maf );
my ( $vcf_tumor_id, $vcf_normal_id );
GetOptions(
    'help!' => \$help,
    'man!' => \$man,
    'input-vcf=s' => \$input_vcf,
    'use-snpeff!' => \$use_snpeff,
    'input-vep=s' => \$vep_anno,
    'input-snpeff=s' => \$snpeff_anno,
    'output-maf=s' => \$output_maf,
    'tumor-id=s' => \$tumor_id,
    'normal-id=s' => \$normal_id,
    'vcf-tumor-id=s' => \$vcf_tumor_id,
    'vcf-normal-id=s' => \$vcf_normal_id,
    'vep-path=s' => \$vep_path,
    'vep-data=s' => \$vep_data,
    'snpeff-path=s' => \$snpeff_path,
    'snpeff-data=s' => \$snpeff_data
) or pod2usage( -verbose => 1, -input => \*DATA, -exitval => 2 );
pod2usage( -verbose => 1, -input => \*DATA, -exitval => 0 ) if( $help );
pod2usage( -verbose => 2, -input => \*DATA, -exitval => 0 ) if( $man );

# Check for valid input combinations
if(( $input_vcf and $vep_anno ) or ( $input_vcf and $snpeff_anno ) or ( $vep_anno and $snpeff_anno )) {
    die "Please specify only one input file: input-vcf, input-vep, or input-snpeff\n";
}
elsif( $use_snpeff and ( $vep_anno or $snpeff_anno )) {
    die "The use-snpeff option can only be used with input-vcf\n";
}

# Unless specified, assume that the VCF uses the same sample IDs that the MAF will contain
$vcf_tumor_id = $tumor_id unless( $vcf_tumor_id );
$vcf_normal_id = $normal_id unless( $vcf_normal_id );

# Annotate variants in given VCF to all possible transcripts, unless an annotated VCF was provided
if( $vep_anno ) {
    ( -s $vep_anno ) or die "Provided VEP-annotated VCF file is missing or empty!\nPath: $vep_anno\n";
    ( $vep_anno !~ m/\.(gz|bz2|bcf)$/ ) or die "Compressed or binary VCFs are not supported\n";
}
elsif( $snpeff_anno ) {
    ( -s $snpeff_anno ) or die "Provided snpEff-annotated VCF file is missing or empty!\nPath: $snpeff_anno\n";
    ( $snpeff_anno !~ m/\.(gz|bz2|bcf)$/ ) or die "Compressed or binary VCFs are not supported\n";
}
elsif( $input_vcf ) {
    ( -s $input_vcf ) or die "Provided VCF file is missing or empty!\nPath: $input_vcf\n";
    ( $input_vcf !~ m/\.(gz|bz2|bcf)$/ ) or die "Compressed or binary VCFs are not supported\n";

    # Run snpEff if user specifically asks for it. Otherwise, run VEP by default
    if( $use_snpeff ) {
        $snpeff_anno = $input_vcf;
        $snpeff_anno =~ s/(\.vcf)*$/.snpeff.vcf/;
        print STDERR "Running snpEff on VCF, and writing output to: $snpeff_anno\n";

        # Make sure we can find the snpEff jar file and config
        unless( -e "$snpeff_path/snpEff.jar" and -e "$snpeff_path/snpEff.config" ) {
            die "Cannot find snpEff jar or config in path: $snpeff_path";
        }

        # Contruct snpEff command using our chosen defaults and run it
        my $snpeff_cmd = "java -Xmx4g -jar $snpeff_path/snpEff.jar eff -config $snpeff_path/snpEff.config -dataDir $snpeff_data -noStats -sequenceOntolgy -hgvs GRCh37.75 $input_vcf > $snpeff_anno";
        print STDERR `$snpeff_cmd`;
        ( -s $snpeff_anno ) or die "snpEff-annotated VCF file is missing or empty!\nPath: $snpeff_anno\n";
    }
    else {
        $vep_anno = $input_vcf;
        $vep_anno =~ s/(\.vcf)*$/.vep.vcf/;
        print STDERR "Running VEP on VCF, and writing output to: $vep_anno\n";

        # Contruct VEP command using our chosen defaults and run it
        my $vep_cmd = "perl $vep_path/variant_effect_predictor.pl --offline --no_stats --everything --xref_refseq --check_existing --total_length --allele_number --no_escape --fork 4 --dir $vep_data --fasta $vep_data --vcf --input_file $input_vcf --output_file $vep_anno";
        print STDERR `$vep_cmd`;
        ( -s $vep_anno ) or die "VEP-annotated VCF file is missing or empty!\nPath: $vep_anno\n";
    }
}
else {
    die "Please specify an input file: input-vcf, input-vep, or input-snpeff. STDIN is not supported.\n";
}

# Define default MAF Header (https://wiki.nci.nih.gov/x/eJaPAQ) with our vcf2maf additions
my @maf_header = qw(
    Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome Start_Position End_Position Strand
    Variant_Classification Variant_Type Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2
    dbSNP_RS dbSNP_Val_Status Tumor_Sample_Barcode Matched_Norm_Sample_Barcode
    Match_Norm_Seq_Allele1 Match_Norm_Seq_Allele2 Tumor_Validation_Allele1 Tumor_Validation_Allele2
    Match_Norm_Validation_Allele1 Match_Norm_Validation_Allele2 Verification_Status Validation_Status
    Mutation_Status Sequencing_Phase Sequence_Source Validation_Method Score BAM_File Sequencer
    Tumor_Sample_UUID Matched_Norm_Sample_UUID HGVSc HGVSp Transcript_ID Exon_Number
    t_depth t_ref_count t_alt_count n_depth n_ref_count n_alt_count
);

# Add extra columns to the MAF depending on whether we used VEP or snpEff
my @vepcsq_cols = qw( Allele Gene Feature Feature_type Consequence cDNA_position CDS_position Protein_position Amino_acids Codons Existing_variation AA_MAF EA_MAF ALLELE_NUM RefSeq EXON INTRON MOTIF_NAME MOTIF_POS HIGH_INF_POS MOTIF_SCORE_CHANGE DISTANCE STRAND CLIN_SIG CANONICAL SYMBOL SYMBOL_SOURCE SIFT PolyPhen GMAF BIOTYPE ENSP DOMAINS CCDS HGVSc HGVSp AFR_MAF AMR_MAF ASN_MAF EUR_MAF PUBMED );
my @snpeff_cols = qw( Effect Effect_Impact Functional_Class Codon_Change Amino_Acid_Change Amino_Acid_Length Gene_Name Transcript_BioType Gene_Coding Transcript_ID Exon_Rank Genotype_Number ERRORS WARNINGS );
push( @maf_header, ( $vep_anno ? @vepcsq_cols : @snpeff_cols ));

# Parse through each variant in the annotated VCF, pull out CSQ/EFF from the INFO column, and choose
# one transcript per variant whose annotation will be used in the MAF
my $vcf_file = ( $vep_anno ? $vep_anno : $snpeff_anno );
my $vcf_fh = IO::File->new( $vcf_file ) or die "Couldn't open annotated VCF: $vcf_file! $!";
my $maf_fh = *STDOUT; # Use STDOUT if an output MAF file was not defined
$maf_fh = IO::File->new( $output_maf, ">" ) or die "Couldn't open output file: $output_maf! $!" if( $output_maf );
$maf_fh->print( "#version 2.4\n" . join( "\t", @maf_header ), "\n" ); # Print MAF header
my ( $vcf_tumor_idx, $vcf_normal_idx );
while( my $line = $vcf_fh->getline ) {

    # Skip all VCF header descriptions cuz we're l33t!
    next if( $line =~ m/^##/ );

    chomp( $line );
    my ( $chrom, $pos, $ids, $ref, $alt, $qual, $filter, $info_line, $format_line, @rest ) = split( /\t/, $line );

    # If FORMATted genotype fields are available, find the sample with the variant, and matched normal
    if( $line =~ m/^#CHROM/ ) {
        if( $format_line and scalar( @rest ) > 0 ) {
            for( my $i = 0; $i <= $#rest; ++$i ) {
                $vcf_tumor_idx = $i if( $rest[$i] eq $vcf_tumor_id );
                $vcf_normal_idx = $i if( $rest[$i] eq $vcf_normal_id );
            }
            ( defined $vcf_tumor_idx ) or die "No genotypes for $vcf_tumor_id in VCF!\n";
        }
        next;
    }

    # Parse out the data in the info column, and store into a hash
    my %info = map {( m/=/ ? ( split( /=/, $_, 2 )) : ( $_, 1 ))} split( /\;/, $info_line );

    # If there are multiple ALT alleles, choose the one with the most supporting reads
    # Or if allelic depths are unavailable, choose the first one listed under ALT
    my @alleles = ( $ref, split( /,/, $alt ));
    my $alt_allele_idx = 1;
    my ( %var_sample_info, %nrm_sample_info );
    if( $format_line and scalar( @rest ) > 0 ) {

        # Load into a hash, the FORMATted genotype info for the sample containing the variant
        my @format_keys = split( /\:/, $format_line );
        my ( $t_idx, $n_idx ) = ( 0, 0 );
        %var_sample_info = map {( $format_keys[$t_idx++], $_ )} split( /\:/, $rest[$vcf_tumor_idx] );
        %nrm_sample_info = map {( $format_keys[$n_idx++], $_ )} split( /\:/, $rest[$vcf_normal_idx] ) if( defined $vcf_normal_idx );
        if( defined $var_sample_info{AD} ) {
            my @allelic_depth = split( /,/, $var_sample_info{AD} );
            # The first depth listed belongs to the reference allele. Of the rest, find the largest
            for( my $i = 1; $i <= $#allelic_depth; ++$i ) {
                $alt_allele_idx = $i if( $allelic_depth[$i] > $allelic_depth[$alt_allele_idx] );
            }
        }
    }
    my ( $var ) = $alleles[$alt_allele_idx];

    # Figure out the appropriate start/stop loci and variant type/allele to report in the MAF
    my $start = my $stop = my $var_type = "";
    my ( $ref_length, $var_length ) = ( length( $ref ), length( $var ));
    # Handle SNPs, DNPs, TNPs, or anything larger (ONP)
    if( $ref_length == $var_length ) {
        ( $start, $stop ) = ( $pos, $pos + $var_length - 1 );
        my %np_type = qw( 1 SNP 2 DNP 3 TNP );
        $var_type = ( $var_length > 3 ? "ONP" : $np_type{$var_length} );
    }
    # Handle all indels, including those complex ones which contain substitutions
    elsif( $ref_length != $var_length ) {
        # Remove the prefixed reference bp from all alleles, using "-" for simple indels
        ( $ref, $var, @alleles ) = map{$_ = substr( $_, 1 ); ( $_ ? $_ : "-" )} ( $ref, $var, @alleles );
        --$ref_length; --$var_length;
        if( $ref_length < $var_length ) { # Handle insertions, and the special case for complex ones
            ( $start, $stop ) = ( $pos, ( $ref eq "-" ? $pos + 1 : $pos + $ref_length ));
            $var_type = "INS";
        }
        else { # Handle deletions
            ( $start, $stop ) = ( $pos + 1, $pos + $ref_length );
            $var_type = "DEL";
        }
    }

    my @all_effects = (); # A list of all effects per variant that can be reported in extra MAF columns
    my $maf_effect = {}; # A single effect per variant to report in the standard MAF columns

    ### Parsing VEP consequences
    # INFO:CSQ is a comma-delimited list of VEP consequences, with pipe-delim details per consequence
    # VEP replaces commas in details with '&'. We'll assume that all '&'s we see, were formerly commas
    # Consequence often reports multiple effects on the same transcript e.g. missense, splice_region
    # CSQ = Allele | Gene | Feature | Feature_type | Consequence | cDNA_position | CDS_position | Protein_position | Amino_acids | Codons | Existing_variation | AA_MAF | EA_MAF | ALLELE_NUM | RefSeq | EXON | INTRON | MOTIF_NAME | MOTIF_POS | HIGH_INF_POS | MOTIF_SCORE_CHANGE | DISTANCE | STRAND | CLIN_SIG | CANONICAL | SYMBOL | SYMBOL_SOURCE | SIFT | PolyPhen | GMAF | BIOTYPE | ENSP | DOMAINS | CCDS | HGVSc | HGVSp | AFR_MAF | AMR_MAF | ASN_MAF | EUR_MAF | PUBMED , ...
    if( $info{CSQ} ) {

        foreach my $csq_line ( split( /,/, $info{CSQ} )) {
            my $idx = 0;
            my %effect = map{s/\&/,/g; ( $vepcsq_cols[$idx++], $_ )} split( /\|/, $csq_line );

            # Skip effects on other ALT alleles
            if( $effect{ALLELE_NUM} == $alt_allele_idx ) {

                # Fix potential warnings about undefined variables
                $effect{BIOTYPE} = '' unless( $effect{BIOTYPE} );
                $effect{Consequence} = '' unless( $effect{Consequence} );
                $effect{CANONICAL} = '' unless( $effect{CANONICAL} );
                $effect{SYMBOL} = '' unless( $effect{SYMBOL} );

                # Remove transcript ID from HGVS codon/protein changes, to make it easier on the eye
                $effect{HGVSc} =~ s/^.*:// if( $effect{HGVSc} );
                $effect{HGVSp} =~ s/^.*:// if( $effect{HGVSp} );

                # Transcript length isn't directly reported, but can be parsed out from another field
                ( $effect{Transcript_Length} ) = $effect{cDNA_position} =~ m/\/(\d+)$/;
                $effect{Transcript_Length} = 0 unless( $effect{Transcript_Length} );

                # If there are many possible consequences on a transcript, choose the most severe one
                ( $effect{Consequence} ) = sort { GetEffectPriority($a) <=> GetEffectPriority($b) } split( /,/, $effect{Consequence} );

                push( @all_effects, \%effect );
            }
        }

        # Sort effects first by transcript biotype, then by severity, and then by longest transcript
        @all_effects = sort {
            GetBiotypePriority( $a->{BIOTYPE} ) <=> GetBiotypePriority( $b->{BIOTYPE} ) ||
            GetEffectPriority( $a->{Consequence} ) <=> GetEffectPriority( $b->{Consequence} ) ||
            $b->{Transcript_Length} <=> $a->{Transcript_Length}
        } @all_effects;

        # For the MAF, we will report the effect on the canonical transcript of the first priority gene
        my $maf_gene = $all_effects[0]->{SYMBOL};
        if( $maf_gene and $maf_gene ne '' ) {
            ( $maf_effect ) = grep { $_->{SYMBOL} eq $maf_gene and $_->{CANONICAL} eq "YES" } @all_effects;

            # If that gene had no canonical transcript tagged, choose the longest transcript instead
            unless( $maf_effect ) {
                ( $maf_effect ) = sort { $b->{Transcript_Length} <=> $a->{Transcript_Length} } grep { $_->{SYMBOL} eq $maf_gene } @all_effects;
            }
        }
        # If the top priority effect doesn't have a gene name, then use the first one that does
        else {
            ( $maf_effect ) = grep { defined $_->{SYMBOL} } @all_effects;
            # If VEP still fails to provide any annotation, it's intergenic
            $maf_effect->{Consequence} = "intergenic_variant" unless( $maf_effect->{Consequence} );
        }
    }

    ### Parsing snpEff effects
    # INFO:EFF is a comma-delimited list of snpEff effects, with pipe-delim details per effect
    # But note the parentheses, the Effect defined separately, and the last two columns being optional
    # Only AA lengths for coding transcripts are provided. So we're SOL for non-coding transcripts
    # EFF = Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change | Amino_Acid_Length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype_Number [ | ERRORS | WARNINGS ] ) , ...
    elsif( $info{EFF} ) {

        foreach my $eff_line ( split( /,/, $info{EFF} )) {
            if( $eff_line =~ /^(\w+)\((.+)\)$/ ) {
                my $idx = 0;
                my %effect = map{( $snpeff_cols[$idx++], $_ )} ( $1, split( /\|/, $2 ));

                # Skip transcripts with errors/warnings, or effects on other ALT alleles
                unless( $effect{ERRORS} or $effect{WARNINGS} or $effect{Genotype_Number} != $alt_allele_idx ) {

                    # Fix potential warnings about undefined variables
                    $effect{Transcript_BioType} = '' unless( $effect{Transcript_BioType} );
                    $effect{Effect} = '' unless( $effect{Effect} );
                    $effect{Gene_Name} = '' unless( $effect{Gene_Name} );

                    # HGVS formatted codon/protein changes need to be parsed out of Amino_Acid_Change
                    ( $effect{HGVSp}, $effect{HGVSc} ) = $effect{Amino_Acid_Change} =~ m/^(.*)\/(.*)$/;

                    # Transcript length isn't reported, so we have to use AA length, where available
                    $effect{Amino_Acid_Length} = 0 unless( $effect{Amino_Acid_Length} );

                    push( @all_effects, \%effect );
                }
            }
        }

        # Sort effects first by transcript biotype, then by severity, and then by longest transcript
        @all_effects = sort {
            GetBiotypePriority( $a->{Transcript_BioType} ) <=> GetBiotypePriority( $b->{Transcript_BioType} ) ||
            GetEffectPriority( $a->{Effect} ) <=> GetEffectPriority( $b->{Effect} ) ||
            $b->{Amino_Acid_Length} <=> $a->{Amino_Acid_Length}
        } @all_effects;

        # For the MAF, we will report the effect on the longest transcript of the first priority gene
        my $maf_gene = $all_effects[0]->{Gene_Name};
        if( $maf_gene and $maf_gene ne '' ) {
            ( $maf_effect ) = sort { $b->{Amino_Acid_Length} <=> $a->{Amino_Acid_Length} } grep { $_->{Gene_Name} eq $maf_gene } @all_effects;
        }
        # If the top priority effect doesn't have a gene name, then use the first one that does
        else {
            ( $maf_effect ) = grep { defined $_->{Gene_Name} } @all_effects;
            # If snpEff still fails to provide any annotation, it's intergenic
            $maf_effect->{Effect} = "intergenic_variant" unless( $maf_effect->{Effect} );
        }
    }

    # Construct the MAF columns from the $maf_effect hash, and print to output
    my %maf_line = map{ ( $_, ( $maf_effect->{$_} ? $maf_effect->{$_} : '' )) } @maf_header;
    $maf_line{Hugo_Symbol} = ( $maf_effect->{SYMBOL} ? $maf_effect->{SYMBOL} : ( $maf_effect->{Gene_Name} ? $maf_effect->{Gene_Name} : 'Unknown' ));
    $maf_line{Entrez_Gene_Id} = '0';
    $maf_line{Center} = $maf_center;
    $maf_line{NCBI_Build} = $ncbi_build;
    $maf_line{Chromosome} = $chrom;
    $maf_line{Start_Position} = $start;
    $maf_line{End_Position} = $stop;
    $maf_line{Strand} = '+';
    my $so_effect = ( $maf_effect->{Consequence} ? $maf_effect->{Consequence} : $maf_effect->{Effect} );
    $maf_line{Variant_Classification} = GetVariantClassification( $so_effect, $var_type);
    $maf_line{Variant_Type} = $var_type;
    $maf_line{Reference_Allele} = $ref;
    # If the genotypes are unavailable, then we'll assume it's ref/var heterozygous
    $maf_line{Tumor_Seq_Allele1} = $ref;
    if( defined $var_sample_info{GT} and $var_sample_info{GT} ne "." ) {
        # ::NOTE:: MAF format doesn't support triallelic genotypes. So break it down as necessary
        my ( $idx1, $idx2 ) = split( /\//, $var_sample_info{GT} );
        $maf_line{Tumor_Seq_Allele1} = ( $alleles[$idx1] eq $var ? $alleles[$idx2] : $alleles[$idx1] );
    }
    $maf_line{Tumor_Seq_Allele2} = $var;
    $maf_line{dbSNP_RS} = GetrsIDs( $ids );
    $maf_line{Tumor_Sample_Barcode} = $tumor_id;
    $maf_line{Matched_Norm_Sample_Barcode} = $normal_id;
    $maf_line{HGVSc} = ( $maf_effect->{HGVSc} ? $maf_effect->{HGVSc} : '' );
    $maf_line{HGVSp} = ( $maf_effect->{HGVSp} ? $maf_effect->{HGVSp} : '' );
    $maf_line{Transcript_ID} = ( $maf_effect->{RefSeq} ? $maf_effect->{RefSeq} : ( $maf_effect->{Transcript_ID} ? $maf_effect->{Transcript_ID} : '' ));
    $maf_line{Exon_Number} = ( $maf_effect->{EXON} ? $maf_effect->{EXON} : ( $maf_effect->{Exon_Rank} ? $maf_effect->{Exon_Rank} : '' ));
    $maf_line{t_depth} = $var_sample_info{DP} if( defined $var_sample_info{DP} );
    ( $maf_line{t_ref_count}, $maf_line{t_alt_count} ) = split( /,/, $var_sample_info{AD} ) if( defined $var_sample_info{AD} );
    $maf_line{n_depth} = $nrm_sample_info{DP} if( defined $nrm_sample_info{DP} );
    ( $maf_line{n_ref_count}, $maf_line{n_alt_count} ) = split( /,/, $nrm_sample_info{AD} ) if( defined $nrm_sample_info{AD} );

    foreach my $col ( @maf_header ) {
        $maf_fh->print( "\t" ) if ( $col ne $maf_header[0] );
        $maf_fh->print( $maf_line{$col} );
        $maf_fh->print( "\n" ) if ( $col eq $maf_header[$#maf_header] );
    }

}
$maf_fh->close if( $output_maf );
$vcf_fh->close;

# Prioritize Sequence Ontology terms from VEP/snpEff in order of severity, as estimated by Ensembl:
# http://useast.ensembl.org/info/genome/variation/predicted_data.html#consequences
# ::NOTE:: snpEff conversion to SO terms has caveats, so handle exceptions as necessary
sub GetEffectPriority {
    my ( $effect ) = @_;
    my %effectPriority = (
        'transcript_ablation' => 1, # A feature ablation whereby the deleted region includes a transcript feature
        'splice_donor_variant' => 2, # A splice variant that changes the 2 base region at the 5' end of an intron
        'splice_acceptor_variant' => 2, # A splice variant that changes the 2 base region at the 3' end of an intron
        'stop_gained' => 3, # A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript
        'frameshift_variant' => 3, # A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three
        'stop_lost' => 3, # A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript
        'initiator_codon_variant' => 4, # A codon variant that changes at least one base of the first codon of a transcript
        'inframe_insertion' => 5, # An inframe non synonymous variant that inserts bases into in the coding sequence
        'inframe_deletion' => 5, # An inframe non synonymous variant that deletes bases from the coding sequence
        'missense_variant' => 6, # A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved
        'transcript_amplification' => 7, # A feature amplification of a region containing a transcript
        'splice_region_variant' => 8, # A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron
        'incomplete_terminal_codon_variant' => 9, # A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed
        'synonymous_variant' => 10, # A sequence variant where there is no resulting change to the encoded amino acid
        'stop_retained_variant' => 10, # A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
        'coding_sequence_variant' => 11, # A sequence variant that changes the coding sequence
        'mature_miRNA_variant' => 11, # A transcript variant located with the sequence of the mature miRNA
        '5_prime_UTR_variant' => 12, # A UTR variant of the 5' UTR
        '5_prime_UTR_premature_start_codon_gain_variant' => 12, # snpEff-specific effect, creating a start codon in 5' UTR
        '3_prime_UTR_variant' => 12, # A UTR variant of the 3' UTR
        'non_coding_exon_variant' => 13, # A sequence variant that changes non-coding exon sequence
        'non_coding_transcript_exon_variant' => 13, # snpEff-specific synonym for non_coding_exon_variant
        'nc_transcript_variant' => 14, # A transcript variant of a non coding RNA
        'intron_variant' => 14, # A transcript variant occurring within an intron
        'INTRAGENIC' => 14, # snpEff-specific effect where the variant hits a gene without transcripts??
        'NMD_transcript_variant' => 15, # A variant in a transcript that is the target of NMD
        'upstream_gene_variant' => 16, # A sequence variant located 5' of a gene
        'downstream_gene_variant' => 16, # A sequence variant located 3' of a gene
        'TFBS_ablation' => 17, # A feature ablation whereby the deleted region includes a transcription factor binding site
        'TFBS_amplification' => 17, # A feature amplification of a region containing a transcription factor binding site
        'TF_binding_site_variant' => 17, # A sequence variant located within a transcription factor binding site
        'regulatory_region_variant' => 17, # A sequence variant located within a regulatory region
        'regulatory_region_ablation' => 17, # A feature ablation whereby the deleted region includes a regulatory region
        'regulatory_region_amplification' => 17, # A feature amplification of a region containing a regulatory region
        'feature_elongation' => 18, # A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence
        'feature_truncation' => 18, # A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence
        'intergenic_variant' => 19, # A sequence variant located in the intergenic region, between genes
        '' => 20
    );
    $effectPriority{$effect} or die "ERROR: Unrecognized effect \"$effect\". Please update your hashes!";
    return $effectPriority{$effect};
}

# Prioritize the transcript biotypes that variants are annotated to, based on disease significance:
# All possible biotypes are defined here: http://www.gencodegenes.org/gencode_biotypes.html
sub GetBiotypePriority {
    my ( $biotype ) = @_;
    my %biotype_priority = (
        'protein_coding' => 1, # Contains an open reading frame (ORF)
        'LRG_gene' => 2, # Gene in a "Locus Reference Genomic" region known to have disease-related sequence variations
        'IG_C_gene' => 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
        'IG_D_gene' => 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
        'IG_J_gene' => 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
        'IG_V_gene' => 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
        'TR_C_gene' => 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
        'TR_D_gene' => 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
        'TR_J_gene' => 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
        'TR_V_gene' => 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
        'miRNA' => 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'snRNA' => 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'snoRNA' => 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'rRNA' => 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'lincRNA' => 3, # Long, intervening noncoding (linc) RNAs, that can be found in evolutionarily conserved, intergenic regions
        'Mt_tRNA' => 4, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'Mt_rRNA' => 4, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'antisense' => 5, # Has transcripts that overlap the genomic span (i.e. exon or introns) of a protein-coding locus on the opposite strand
        'sense_intronic' => 5, # Long non-coding transcript in introns of a coding gene that does not overlap any exons
        'sense_overlapping' => 5, # Long non-coding transcript that contains a coding gene in its intron on the same strand
        '3prime_overlapping_ncrna' => 5, # Transcripts where ditag and/or published experimental data strongly supports the existence of short non-coding transcripts transcribed from the 3'UTR
        'misc_RNA' => 5, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'non_coding' => 5, # Transcript which is known from the literature to not be protein coding
        'disrupted_domain' => 6, # Otherwise viable coding region omitted from this alternatively spliced transcript because the splice variation affects a region coding for a protein domain
        'processed_transcript' => 6, # Doesn't contain an ORF
        'TEC' => 6, # To be Experimentally Confirmed. This is used for non-spliced EST clusters that have polyA features. This category has been specifically created for the ENCODE project to highlight regions that could indicate the presence of protein coding genes that require experimental validation, either by 5' RACE or RT-PCR to extend the transcripts, or by confirming expression of the putatively-encoded peptide with specific antibodies
        'retained_intron' => 7, # Alternatively spliced transcript believed to contain intronic sequence relative to other, coding, variants
        'nonsense_mediated_decay' => 7, # If the coding sequence (following the appropriate reference) of a transcript finishes >50bp from a downstream splice site then it is tagged as NMD. If the variant does not cover the full reference coding sequence then it is annotated as NMD if NMD is unavoidable i.e. no matter what the exon structure of the missing portion is the transcript will be subject to NMD
        'non_stop_decay' => 7, # Transcripts that have polyA features (including signal) without a prior stop codon in the CDS, i.e. a non-genomic polyA tail attached directly to the CDS without 3' UTR. These transcripts are subject to degradation
        'ambiguous_orf' => 7, # Transcript believed to be protein coding, but with more than one possible open reading frame
        'pseudogene' => 8, # Have homology to proteins but generally suffer from a disrupted coding sequence and an active homologous gene can be found at another locus. Sometimes these entries have an intact coding sequence or an open but truncated ORF, in which case there is other evidence used (for example genomic polyA stretches at the 3' end) to classify them as a pseudogene. Can be further classified as one of the following
        'processed_pseudogene' => 8, # Pseudogene that lack introns and is thought to arise from reverse transcription of mRNA followed by reinsertion of DNA into the genome
        'polymorphic_pseudogene' => 8, # Pseudogene owing to a SNP/DIP but in other individuals/haplotypes/strains the gene is translated
        'retrotransposed' => 8, # Pseudogene owing to a reverse transcribed and re-inserted sequence
        'translated_processed_pseudogene' => 8, # Pseudogenes that have mass spec data suggesting that they are also translated
        'translated_unprocessed_pseudogene' => 8, # Pseudogenes that have mass spec data suggesting that they are also translated
        'transcribed_processed_pseudogene' => 8, # Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression
        'transcribed_unprocessed_pseudogene' => 8, # Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression
        'unitary_pseudogene' => 8, # A species specific unprocessed pseudogene without a parent gene, as it has an active orthologue in another species
        'unprocessed_pseudogene' => 8, # Pseudogene that can contain introns since produced by gene duplication
        'Mt_tRNA_pseudogene' => 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
        'tRNA_pseudogene' => 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
        'snoRNA_pseudogene' => 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
        'snRNA_pseudogene' => 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
        'scRNA_pseudogene' => 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
        'rRNA_pseudogene' => 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
        'misc_RNA_pseudogene' => 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
        'miRNA_pseudogene' => 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
        'IG_C_pseudogene' => 8, # Inactivated immunoglobulin gene
        'IG_J_pseudogene' => 8, # Inactivated immunoglobulin gene
        'IG_V_pseudogene' => 8, # Inactivated immunoglobulin gene
        'TR_J_pseudogene' => 8, # Inactivated immunoglobulin gene
        'TR_V_pseudogene' => 8, # Inactivated immunoglobulin gene
        'artifact' => 9, # Used to tag mistakes in the public databases (Ensembl/SwissProt/Trembl)
        '' => 9
    );
    $biotype_priority{$biotype} or die "ERROR: Unrecognized biotype \"$biotype\". Please update your hashes!";
    return $biotype_priority{$biotype};
}

# Converts Sequence Ontology variant types to MAF variant classifications
sub GetVariantClassification {
    my ( $effect, $var_type ) = @_;
    return "Splice_Site" if( $effect =~ /^(splice_acceptor_variant|splice_donor_variant|transcript_ablation)$/ );
    return "Nonsense_Mutation" if( $effect eq 'stop_gained' );
    return "Frame_Shift_Del" if ( $effect eq 'frameshift_variant' and $var_type eq 'DEL' );
    return "Frame_Shift_Ins" if( $effect eq 'frameshift_variant' and $var_type eq 'INS' );
    return "Nonstop_Mutation" if( $effect eq 'stop_lost' );
    return "Translation_Start_Site" if( $effect eq 'initiator_codon_variant' );
    return "In_Frame_Ins" if( $effect eq 'inframe_insertion' );
    return "In_Frame_Del" if( $effect eq 'inframe_deletion' );
    return "Missense_Mutation" if( $effect =~ /^(missense_variant|coding_sequence_variant)$/ );
    return "Intron" if ( $effect =~ /^(transcript_amplification|splice_region_variant|intron_variant|INTRAGENIC)$/ );
    return "Silent" if( $effect =~ /^(incomplete_terminal_codon_variant|synonymous_variant|stop_retained_variant|NMD_transcript_variant)$/ );
    return "RNA" if( $effect =~ /^(mature_miRNA_variant|non_coding_exon_variant|non_coding_transcript_exon_variant|nc_transcript_variant)$/ );
    return "5'UTR" if( $effect =~ /^(5_prime_UTR_variant|5_prime_UTR_premature_start_codon_gain_variant)$/ );
    return "3'UTR" if( $effect eq '3_prime_UTR_variant' );
    return "IGR" if( $effect =~ /^(TF_binding_site_variant|regulatory_region_variant|intergenic_variant)$/ );
    return "5'Flank" if( $effect eq 'upstream_gene_variant' );
    return "3'Flank" if ( $effect eq 'downstream_gene_variant' );

    # Annotate everything else simply as a targeted region
    # TFBS_ablation, TFBS_amplification,regulatory_region_ablation, regulatory_region_amplification,
    # feature_elongation, feature_truncation
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

 vcf2maf.pl - Map effects of variants in a given VCF, and report them in a MAF

=head1 SYNOPSIS

 perl vcf2maf.pl --help
 perl vcf2maf.pl --input-vcf test.vcf --output-maf test.maf --tumor-id WD4086 --normal-id NB4086

=head1 OPTIONS

 --input-vcf      Path to input file in VCF format
 --input-vep      Path to VEP-annotated VCF file
 --input-snpeff   Path to snpEff-annotated VCF file
 --output-maf     Path to output MAF file [Default: STDOUT]
 --tumor-id       Tumor_Sample_Barcode to report in the MAF [TUMOR]
 --normal-id      Matched_Norm_Sample_Barcode to report in the MAF [NORMAL]
 --vcf-tumor-id   Tumor sample ID used in VCF's genotype columns [--tumor-id]
 --vcf-normal-id  Matched normal ID used in VCF's genotype columns [--normal-id]
 --use-snpeff     Use snpEff to annotate VCF, instead of the default VEP
 --vep-path       Folder containing variant_effect_predictor.pl [~/vep]
 --vep-data       VEP's base cache/plugin directory [~/.vep]
 --snpeff-path    Folder containing snpEff.jar and snpEff.config [~/snpEff]
 --snpeff-data    Override for data_dir in snpEff.config [~/snpEff/data]
 --help           Print a brief help message and quit
 --man            Print the detailed manual

=head1 DESCRIPTION

To convert a VCF into a MAF, each variant must be annotated to only one of all possible gene transcripts/isoforms that it might affect. This selection of a single effect on a single transcript, per variant, is often subjective. So this script tries to follow best practices, but makes the selection criteria smarter, reproducible, and more configurable.

This script needs Ensembl's VEP or snpEff - variant annotators that can map the effect of a variant on all possible gene transcripts in a database. For more info, see the README.

=head2 Relevant links:

 vcf2maf Homepage: https://github.com/ckandoth/vcf2maf
 VCF format: http://samtools.github.io/hts-specs/
 MAF format: https://wiki.nci.nih.gov/x/eJaPAQ
 Ensembl VEP: http://ensembl.org/info/docs/tools/vep/index.html
 VEP annotated VCF format: http://ensembl.org/info/docs/tools/vep/vep_formats.html#vcfout
 snpEff: http://snpeff.sourceforge.net
 snpEff annotated VCF format: http://snpeff.sourceforge.net/SnpEff_manual.html#output

=head1 AUTHORS

 Cyriac Kandoth (ckandoth@gmail.com)
 William Lee, Senior Research Scientist, Memorial Sloan Kettering Cancer Center

=head1 LICENSE

 LGPLv3, Memorial Sloan-Kettering Cancer Center, New York, NY 10065, USA

=cut
