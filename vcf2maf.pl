#!/usr/bin/env perl

# vcf2maf - Convert a VCF into a MAF by mapping each variant to only one of all possible gene isoforms

use strict;
use warnings;
use IO::File;
use Getopt::Long qw( GetOptions );
use Pod::Usage qw( pod2usage );
use File::Copy qw( move );
use File::Path qw( mkpath );
use Config;
use Text::Wrap;

# Set any default paths and constants
my ( $tumor_id, $normal_id ) = ( "TUMOR", "NORMAL" );
my ( $vep_path, $vep_data, $vep_forks, $buffer_size, $any_allele, $inhibit_vep, $online, $vep_custom, $vep_config, $vep_overwrite  ) = ( "$ENV{HOME}/miniconda3/bin", "$ENV{HOME}/.vep", 4, 5000, 0, 0, 0, "", "", 0 );
my ( $ref_fasta ) = ( "$ENV{HOME}/.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa.gz" );
my ( $species, $ncbi_build, $cache_version, $maf_center, $retain_info, $retain_fmt, $retain_ann, $min_hom_vaf, $max_subpop_af ) = ( "homo_sapiens", "GRCh37", "", ".", "", "", "", 0.7, 0.0004 );
my $perl_bin = $Config{perlpath};

# Set default formatting for any output command lines:
$Text::Wrap::huge = 'overflow';
$Text::Wrap::separator = " \\$/";

# Find out if samtools and tabix are properly installed, and warn the user if it's not
my ( $samtools ) = map{chomp; $_}`which samtools`;
( $samtools and -e $samtools ) or die "ERROR: Please install samtools, and make sure it's in your PATH\n";
my ( $tabix ) = map{chomp; $_}`which tabix`;
( $tabix and -e $tabix ) or die "ERROR: Please install tabix, and make sure it's in your PATH\n";

# Hash to convert 3-letter amino-acid codes to their 1-letter codes
my %aa3to1 = qw( Ala A Arg R Asn N Asp D Asx B Cys C Glu E Gln Q Glx Z Gly G His H Ile I Leu L
    Lys K Met M Phe F Pro P Ser S Thr T Trp W Tyr Y Val V Xxx X Ter * );

# Prioritize Sequence Ontology terms in order of severity, as estimated by Ensembl:
# https://ensembl.org/info/genome/variation/prediction/predicted_data.html
sub GetEffectPriority {
    my ( $effect ) = @_;
    $effect = '' unless( defined $effect );
    my %effectPriority = (
        'transcript_ablation' => 1, # A feature ablation whereby the deleted region includes a transcript feature
        'exon_loss_variant' => 1, # A sequence variant whereby an exon is lost from the transcript
        'splice_donor_variant' => 2, # A splice variant that changes the 2 base region at the 5' end of an intron
        'splice_acceptor_variant' => 2, # A splice variant that changes the 2 base region at the 3' end of an intron
        'stop_gained' => 3, # A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript
        'frameshift_variant' => 3, # A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three
        'stop_lost' => 3, # A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript
        'start_lost' => 4, # A codon variant that changes at least one base of the canonical start codon
        'initiator_codon_variant' => 4, # A codon variant that changes at least one base of the first codon of a transcript
        'disruptive_inframe_insertion' => 5, # An inframe increase in cds length that inserts one or more codons into the coding sequence within an existing codon
        'disruptive_inframe_deletion' => 5, # An inframe decrease in cds length that deletes bases from the coding sequence starting within an existing codon
        'conservative_inframe_insertion' => 5, # An inframe increase in cds length that inserts one or more codons into the coding sequence between existing codons
        'conservative_inframe_deletion' => 5, # An inframe decrease in cds length that deletes one or more entire codons from the coding sequence but does not change any remaining codons
        'inframe_insertion' => 5, # An inframe non synonymous variant that inserts bases into the coding sequence
        'inframe_deletion' => 5, # An inframe non synonymous variant that deletes bases from the coding sequence
        'protein_altering_variant' => 5, # A sequence variant which is predicted to change the protein encoded in the coding sequence
        'missense_variant' => 6, # A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved
        'conservative_missense_variant' => 6, # A sequence variant whereby at least one base of a codon is changed resulting in a codon that encodes for a different but similar amino acid. These variants may or may not be deleterious
        'rare_amino_acid_variant' => 6, # A sequence variant whereby at least one base of a codon encoding a rare amino acid is changed, resulting in a different encoded amino acid
        'transcript_amplification' => 7, # A feature amplification of a region containing a transcript
        'splice_region_variant' => 8, # A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron
        'start_retained_variant' => 9, # A sequence variant where at least one base in the start codon is changed, but the start remains
        'stop_retained_variant' => 9, # A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
        'synonymous_variant' => 9, # A sequence variant where there is no resulting change to the encoded amino acid
        'incomplete_terminal_codon_variant' => 10, # A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed
        'coding_sequence_variant' => 11, # A sequence variant that changes the coding sequence
        'mature_miRNA_variant' => 11, # A transcript variant located with the sequence of the mature miRNA
        'exon_variant' => 11, # A sequence variant that changes exon sequence
        '5_prime_UTR_variant' => 12, # A UTR variant of the 5' UTR
        '5_prime_UTR_premature_start_codon_gain_variant' => 12, # snpEff-specific effect, creating a start codon in 5' UTR
        '3_prime_UTR_variant' => 12, # A UTR variant of the 3' UTR
        'non_coding_exon_variant' => 13, # A sequence variant that changes non-coding exon sequence
        'non_coding_transcript_exon_variant' => 13, # snpEff-specific synonym for non_coding_exon_variant
        'non_coding_transcript_variant' => 14, # A transcript variant of a non coding RNA gene
        'nc_transcript_variant' => 14, # A transcript variant of a non coding RNA gene (older alias for non_coding_transcript_variant)
        'intron_variant' => 14, # A transcript variant occurring within an intron
        'intragenic_variant' => 14, # A variant that occurs within a gene but falls outside of all transcript features. This occurs when alternate transcripts of a gene do not share overlapping sequence
        'INTRAGENIC' => 14, # snpEff-specific synonym of intragenic_variant
        'NMD_transcript_variant' => 15, # A variant in a transcript that is the target of NMD
        'upstream_gene_variant' => 16, # A sequence variant located 5' of a gene
        'downstream_gene_variant' => 16, # A sequence variant located 3' of a gene
        'TFBS_ablation' => 17, # A feature ablation whereby the deleted region includes a transcription factor binding site
        'TFBS_amplification' => 17, # A feature amplification of a region containing a transcription factor binding site
        'TF_binding_site_variant' => 17, # A sequence variant located within a transcription factor binding site
        'regulatory_region_ablation' => 17, # A feature ablation whereby the deleted region includes a regulatory region
        'regulatory_region_amplification' => 17, # A feature amplification of a region containing a regulatory region
        'regulatory_region_variant' => 17, # A sequence variant located within a regulatory region
        'regulatory_region' =>17, # snpEff-specific effect that should really be regulatory_region_variant
        'feature_elongation' => 18, # A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence
        'feature_truncation' => 18, # A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence
        'intergenic_variant' => 19, # A sequence variant located in the intergenic region, between genes
        'intergenic_region' => 19, # snpEff-specific effect that should really be intergenic_variant
        '' => 20
    );
    unless( defined $effectPriority{$effect} ) {
        warn "WARNING: Unrecognized effect \"$effect\". Assigning lowest priority!\n";
        return 20;
    }
    return $effectPriority{$effect};
}

# Prioritize the transcript biotypes that variants are annotated to, based on disease significance:
# All possible biotypes are defined here: https://www.gencodegenes.org/pages/biotypes.html
sub GetBiotypePriority {
    my ( $biotype ) = @_;
    $biotype = '' unless( defined $biotype );
    my %biotype_priority = (
        'protein_coding' => 1, # Contains an open reading frame (ORF)
        'LRG_gene' => 2, # Gene in a "Locus Reference Genomic" region known to have disease-related sequence variations
        'IG_C_gene' => 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
        'IG_D_gene' => 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
        'IG_J_gene' => 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
        'IG_LV_gene' => 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
        'IG_V_gene' => 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
        'TR_C_gene' => 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
        'TR_D_gene' => 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
        'TR_J_gene' => 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
        'TR_V_gene' => 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
        'miRNA' => 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'snRNA' => 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'snoRNA' => 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'ribozyme' => 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'tRNA' => 3, #Added by Y. Boursin
        'sRNA' => 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'scaRNA' => 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'rRNA' => 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'scRNA' => 3, # Non-coding RNA predicted using sequences from Rfam and miRBase
        'lincRNA' => 3, # Long, intervening noncoding (linc) RNAs, that can be found in evolutionarily conserved, intergenic regions
        'lncRNA' => 3, # Replaces 3prime_overlapping_ncRNA, antisense, bidirectional_promoter_lncRNA, lincRNA, macro_lncRNA, non_coding, processed_transcript, sense_intronic and sense_overlapping
        'bidirectional_promoter_lncrna' => 3, # A non-coding locus that originates from within the promoter region of a protein-coding gene, with transcription proceeding in the opposite direction on the other strand
        'bidirectional_promoter_lncRNA' => 3, # A non-coding locus that originates from within the promoter region of a protein-coding gene, with transcription proceeding in the opposite direction on the other strand
        'known_ncrna' => 4,
        'vaultRNA' => 4, # Short non coding RNA genes that form part of the vault ribonucleoprotein complex
        'macro_lncRNA' => 4, # unspliced lncRNAs that are several kb in size
        'Mt_tRNA' => 4, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'Mt_rRNA' => 4, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'antisense' => 5, # Has transcripts that overlap the genomic span (i.e. exon or introns) of a protein-coding locus on the opposite strand
        'antisense_RNA' => 5, # Alias for antisense (Y. Boursin)
        'sense_intronic' => 5, # Long non-coding transcript in introns of a coding gene that does not overlap any exons
        'sense_overlapping' => 5, # Long non-coding transcript that contains a coding gene in its intron on the same strand
        '3prime_overlapping_ncrna' => 5, # Transcripts where ditag and/or published experimental data strongly supports the existence of short non-coding transcripts transcribed from the 3'UTR
        '3prime_overlapping_ncRNA' => 5, # Transcripts where ditag and/or published experimental data strongly supports the existence of short non-coding transcripts transcribed from the 3'UTR
        'misc_RNA' => 5, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'non_coding' => 5, # Transcript which is known from the literature to not be protein coding
        'regulatory_region' => 6, # A region of sequence that is involved in the control of a biological process
        'disrupted_domain' => 6, # Otherwise viable coding region omitted from this alternatively spliced transcript because the splice variation affects a region coding for a protein domain
        'processed_transcript' => 6, # Doesn't contain an ORF
        'TEC' => 6, # To be Experimentally Confirmed. This is used for non-spliced EST clusters that have polyA features. This category has been specifically created for the ENCODE project to highlight regions that could indicate the presence of protein coding genes that require experimental validation, either by 5' RACE or RT-PCR to extend the transcripts, or by confirming expression of the putatively-encoded peptide with specific antibodies
        'TF_binding_site' => 7, # A region of a nucleotide molecule that binds a Transcription Factor or Transcription Factor complex
        'CTCF_binding_site' =>7, # A transcription factor binding site with consensus sequence CCGCGNGGNGGCAG, bound by CCCTF-binding factor
        'promoter_flanking_region' => 7, # A region immediately adjacent to a promoter which may or may not contain transcription factor binding sites
        'enhancer' => 7, # A cis-acting sequence that increases the utilization of (some) eukaryotic promoters, and can function in either orientation and in any location (upstream or downstream) relative to the promoter
        'promoter' => 7, # A regulatory_region composed of the TSS(s) and binding sites for TF_complexes of the basal transcription machinery
        'open_chromatin_region' => 7, # A DNA sequence that in the normal state of the chromosome corresponds to an unfolded, un-complexed stretch of double-stranded DNA
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
        'transcribed_unitary_pseudogene' => 8, #Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression
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
        'IG_D_pseudogene' => 8, # Inactivated immunoglobulin gene
        'IG_J_pseudogene' => 8, # Inactivated immunoglobulin gene
        'IG_V_pseudogene' => 8, # Inactivated immunoglobulin gene
        'TR_J_pseudogene' => 8, # Inactivated immunoglobulin gene
        'TR_V_pseudogene' => 8, # Inactivated immunoglobulin gene
        'artifact' => 9, # Used to tag mistakes in the public databases (Ensembl/SwissProt/Trembl)
        '' => 10
    );
    unless( defined $biotype_priority{$biotype} ) {
        warn "WARNING: Unrecognized biotype \"$biotype\". Assigning lowest priority!\n";
        return 10;
    }
    return $biotype_priority{$biotype};
}

# Check for missing or crappy arguments
unless( @ARGV and $ARGV[0] =~ m/^-/ ) {
    pod2usage( -verbose => 0, -message => "$0: Missing or invalid arguments!\n", -exitval => 2 );
}

# Parse options and print usage if there is a syntax error, or if usage was explicitly requested
my ( $man, $help, $verbose ) = ( 0, 0, 0 );
my ( $input_vcf, $output_maf, $tmp_dir, $custom_enst_file );
my ( $vcf_tumor_id, $vcf_normal_id, $remap_chain );
GetOptions(
    'help!' => \$help,
    'man!' => \$man,
    'verbose!' => \$verbose,
    'input-vcf=s' => \$input_vcf,
    'output-maf=s' => \$output_maf,
    'tmp-dir=s' => \$tmp_dir,
    'tumor-id=s' => \$tumor_id,
    'normal-id=s' => \$normal_id,
    'vcf-tumor-id=s' => \$vcf_tumor_id,
    'vcf-normal-id=s' => \$vcf_normal_id,
    'custom-enst=s' => \$custom_enst_file,
    'vep-path=s' => \$vep_path,
    'vep-data=s' => \$vep_data,
    'vep-forks=s' => \$vep_forks,
    'vep-custom=s' => \$vep_custom,
    'vep-config=s' => \$vep_config,
    'vep-overwrite!' => \$vep_overwrite,
    'buffer-size=i' => \$buffer_size,
    'any-allele!' => \$any_allele,
    'inhibit-vep!' => \$inhibit_vep,
    'online!' => \$online,
    'ref-fasta=s' => \$ref_fasta,
    'species=s' => \$species,
    'ncbi-build=s' => \$ncbi_build,
    'cache-version=s' => \$cache_version,
    'maf-center=s' => \$maf_center,
    'retain-info=s' => \$retain_info,
    'retain-fmt=s' => \$retain_fmt,
    'retain-ann=s' => \$retain_ann,
    'min-hom-vaf=s' => \$min_hom_vaf,
    'remap-chain=s' => \$remap_chain,
    'max-subpop-af=f' => \$max_subpop_af
) or pod2usage( -verbose => 1, -input => \*DATA, -exitval => 2 );
pod2usage( -verbose => 1, -input => \*DATA, -exitval => 0 ) if( $help );
pod2usage( -verbose => 2, -input => \*DATA, -exitval => 0 ) if( $man );

# Check if required arguments are missing or problematic
( defined $input_vcf and defined $output_maf ) or die "ERROR: Both input-vcf and output-maf must be defined!\n";
( -s $input_vcf ) or die "ERROR: Provided --input-vcf is missing or empty: $input_vcf\n";
( -s $ref_fasta ) or die "ERROR: Provided --ref-fasta is missing or empty: $ref_fasta\n";
( $input_vcf !~ m/\.(gz|bz2|bcf)$/ ) or die "ERROR: Unfortunately, --input-vcf cannot be in a compressed format\n";

# Unless specified, assume that the VCF uses the same sample IDs that the MAF will contain
$vcf_tumor_id = $tumor_id unless( $vcf_tumor_id );
$vcf_normal_id = $normal_id unless( $vcf_normal_id );

# Load up the custom isoform overrides if provided:
my %custom_enst;
if( $custom_enst_file ) {
    ( -s $custom_enst_file ) or die "ERROR: Provided --custom-enst file is missing or empty: $custom_enst_file\n";
    warn "STATUS: Reading --custom-enst $custom_enst_file...\n" if( $verbose );
    %custom_enst = map{chomp; ( $_, 1 )}`grep -v ^# $custom_enst_file | cut -f1`;
}

# Create a folder for the intermediate VCFs if user-defined, or default to the input VCF's folder
if( defined $tmp_dir ) {
    mkpath( $tmp_dir ) unless( -d $tmp_dir );
}
else {
    $tmp_dir = substr( $input_vcf, 0, rindex( $input_vcf, "/" )) if( $input_vcf =~ m/\// );
    $tmp_dir = "." unless( $tmp_dir ); # In case the input VCF is in the current working directory
}

# Also figure out the base name of the input VCF, cuz we'll be naming a lot of files based on that
my $input_name = substr( $input_vcf, rindex( $input_vcf, "/" ) + 1 );
$input_name =~ s/(\.vcf)*$//;

# If the VCF contains SVs, split the breakpoints into separate lines before passing to VEP
my ( $split_svs, $var_count ) = ( 0, 0 );
my $orig_vcf_fh = IO::File->new( $input_vcf ) or die "ERROR: Couldn't open --input-vcf: $input_vcf!\n";
my $split_vcf_fh = IO::File->new( "$tmp_dir/$input_name.split.vcf", "w" ) or die "ERROR: Couldn't open VCF: $tmp_dir/$input_name.split.vcf!\n";

warn "STATUS: Preprocessing $input_vcf: split SV breakpoints before passing to VEP...\n" if( $verbose );

while( my $line = $orig_vcf_fh->getline ) {
    # If the file uses Mac OS 9 newlines, quit with an error
    ( $line !~ m/\r$/ ) or die "ERROR: Your VCF uses CR line breaks, which we can't support. Please use LF or CRLF.\n";

    if( $line =~ m/^#/ ) {
        $split_vcf_fh->print( $line ); # Write header lines unchanged
        next;
    }

    chomp( $line );
    ++$var_count;
    my @cols = split( "\t", $line );
    my %info = map {( m/=/ ? ( split( /=/, $_, 2 )) : ( $_, "1" ))} split( /\;/, $cols[7] );
    if( $info{SVTYPE} ){
        # Remove SVTYPE tag if REF/ALT alleles are defined, or VEP won't report transcript effects
        if( $cols[3]=~m/^[ACGTN]+$/i and $cols[4]=~m/^[ACGTN,]+$/i ) {
            $cols[7]=~s/(SVTYPE=\w+;|;SVTYPE=\w+|SVTYPE=\w+)//;
            $split_vcf_fh->print( join( "\t", @cols ), "\n" );
        }
        # For legit SVs except insertions, split them into two separate breakpoint events
        elsif( $info{SVTYPE}=~m/^(BND|TRA|DEL|DUP|INV)$/ ) {
            $split_svs = 1;
            # Don't tell VEP it's an SV, by removing the SVTYPE tag
            $cols[7]=~s/(SVTYPE=\w+;|;SVTYPE=\w+|SVTYPE=\w+)//;
            # Rename two SV specific INFO keys to something friendlier
            $cols[7]=~s/CT=([35]to[35])/Frame=$1/;
            $cols[7]=~s/SVMETHOD=([\w.]+)/Method=$1/;
            $cols[4] = "<" . $info{SVTYPE} . ">";
            # Fetch the REF allele at the second breakpoint using samtools faidx
            my $ref2 = `$samtools faidx $ref_fasta $info{CHR2}:$info{END}-$info{END} | grep -v ^\\>`;
            chomp( $ref2 );
            $split_vcf_fh->print( join( "\t", $info{CHR2}, $info{END}, $cols[2], ( $ref2 ? $ref2 : $cols[3] ), @cols[4..$#cols] ), "\n" );
            $split_vcf_fh->print( join( "\t", @cols ), "\n" );
        }
        $input_vcf = "$tmp_dir/$input_name.split.vcf";
    }
    else {
        $split_vcf_fh->print( join( "\t", @cols ), "\n" );
    }
}
$split_vcf_fh->close;
$orig_vcf_fh->close;

# Delete the split.vcf created above if we didn't find any variants with the SVTYPE tag
unlink( "$tmp_dir/$input_name.split.vcf" ) if( $input_vcf ne "$tmp_dir/$input_name.split.vcf" );

# Make sure the --online option is only used with small GRCh38 VCFs
if( $online ) {
    ( $var_count < 100 and $ncbi_build eq "GRCh38" ) or die "ERROR: Option --online can only be used with GRCh38 VCFs listing <100 events\n";
}

# If a liftOver chain was provided, remap and switch the input VCF before annotation
my ( %remap );
if( $remap_chain ) {
    warn "STATUS: Running liftOver...\n" if( $verbose );

    # Find out if liftOver is properly installed, and warn the user if it's not
    my $liftover = `which liftOver`;
    chomp( $liftover );
    ( $liftover and -e $liftover ) or die "ERROR: Please install liftOver, and make sure it's in your PATH\n";

    # Make a BED file from the VCF, run liftOver on it, and create a hash mapping old to new loci
    `grep -v ^# $input_vcf | cut -f1,2 | awk '{OFS="\\t"; print \$1,\$2-1,\$2,\$1":"\$2}' > $tmp_dir/$input_name.bed`;
    %remap = map{chomp; my @c=split("\t"); ($c[3], "$c[0]:$c[2]")}`$liftover $tmp_dir/$input_name.bed $remap_chain /dev/stdout /dev/null 2> /dev/null`;
    unlink( "$tmp_dir/$input_name.bed" );

    # Create a new VCF in the temp folder, with remapped loci on which we'll run annotation
    my $orig_vcf_fh = IO::File->new( $input_vcf ) or die "ERROR: Couldn't open --input-vcf: $input_vcf!\n";
    my $remap_vcf_fh = IO::File->new( "$tmp_dir/$input_name.remap.vcf", "w" ) or die "ERROR: Couldn't open VCF: $tmp_dir/$input_name.remap.vcf!\n";
    while( my $line = $orig_vcf_fh->getline ) {
        if( $line =~ m/^#/ ) {
            $remap_vcf_fh->print( $line ); # Write header lines unchanged
        }
        else {
            chomp( $line );
            my @cols = split( "\t", $line );
            my $locus = $cols[0] . ":" . $cols[1];
            if( defined $remap{$locus} ) {
                # Retain original variant under INFO, so we can append it later to the output MAF
                $cols[7] = ( !$cols[7] or $cols[7] eq "." ? "" : "$cols[7];" ) . "REMAPPED_POS=" . join( ":", @cols[0,1,3,4] );
                @cols[0,1] = split( ":", $remap{$locus} );
                $remap_vcf_fh->print( join( "\t", @cols ), "\n" );
            }
            else {
                warn "WARNING: Skipping variant at $locus; Unable to liftOver using $remap_chain\n";
            }
        }
    }
    $remap_vcf_fh->close;
    $orig_vcf_fh->close;
    $input_vcf = "$tmp_dir/$input_name.remap.vcf";
}

# Before running annotation, let's pull flanking reference bps for each variant to do some checks
warn "STATUS: Pulling flanking reference bps for checks...\n" if( $verbose );
my $vcf_fh = IO::File->new( $input_vcf ) or die "ERROR: Couldn't open --input-vcf: $input_vcf!\n";
my ( %ref_bps, @ref_regions, %uniq_regions, %flanking_bps );
while( my $line = $vcf_fh->getline ) {
    # Skip header lines, and pull variant loci to pass to samtools later
    next if( $line =~ m/^#/ );
    chomp( $line );
    my ( $chr, $pos, undef, $ref ) = split( "\t", $line );
    # Create a region that spans the length of the reference allele and 1bp flanks around it
    my $region = "$chr:" . ( $pos - 1 ) . "-" . ( $pos + length( $ref ));
    $ref_bps{$region} = $ref;
    push( @ref_regions, $region );
    $uniq_regions{$region} = 1;
}
$vcf_fh->close;

# samtools runs faster when passed many loci at a time, but limited to around 125k args, at least
# on CentOS 6. If there are too many loci, split them into smaller chunks and run separately
warn "STATUS: Splitting loci into smaller chunks to run separately...\n" if( $verbose );
my ( $lines, @regions_split ) = ( "", ());
my @regions = keys %uniq_regions;
my $chr_prefix_in_use = ( @regions and $regions[0] =~ m/^chr/ ? 1 : 0 );
push( @regions_split, [ splice( @regions, 0, $buffer_size ) ] ) while @regions;
map{ my $region = join( " ", sort @{$_} ); $lines .= `$samtools faidx $ref_fasta $region` } @regions_split;
foreach my $line ( grep( length, split( ">", $lines ))) {
    # Carefully split this FASTA entry, properly chomping newlines for long indels
    my ( $region, $bps ) = split( "\n", $line, 2 );
    $bps =~ s/\r|\n//g;
    if( $bps ){
        $bps = uc( $bps );
        $flanking_bps{$region} = $bps;
    }
}

# If flanking_bps is entirely empty, then it's most likely that the user chose the wrong ref-fasta
# Or it's also possible that an outdated samtools was unable to parse the gzipped FASTA files
# ::NOTE:: If input had no variants, don't break here, so we can continue to create an empty MAF
( !@regions_split or %flanking_bps ) or die "ERROR: You're either using an outdated samtools, or --ref-fasta is not the same genome build as your --input-vcf.";

# For each variant locus and reference allele in the input VCF, report any problems
warn "STATUS: Reporting any problems on variant loci and reference alleles...\n" if( $verbose );
foreach my $region ( @ref_regions ) {
    my $ref = $ref_bps{$region};
    my ( $locus ) = map{ my ( $chr, $pos ) = split( ":" ); ++$pos; "$chr:$pos" } split( "-", $region );
    if( !defined $flanking_bps{$region} ) {
        warn "WARNING: Couldn't retrieve bps around $locus from reference FASTA: $ref_fasta\n";
    }
    elsif( $flanking_bps{$region} !~ m/^[ACGTN]+$/ ) {
        warn "WARNING: Retrieved invalid bps " . $flanking_bps{$region} . " around $locus from reference FASTA: $ref_fasta\n";
    }
    elsif( $ref ne substr( $flanking_bps{$region}, 1, length( $ref ))) {
        warn "WARNING: Reference allele $ref at $locus doesn't match " .
            substr( $flanking_bps{$region}, 1, length( $ref )) . " (flanking bps: " .
            $flanking_bps{$region} . ") from reference FASTA: $ref_fasta\n";
    }
}

# Annotate variants in given VCF to all possible transcripts, unless user requested to skip VEP
my $output_vcf = $input_vcf;
unless( $inhibit_vep ) {
    $output_vcf = ( $remap_chain ? "$tmp_dir/$input_name.remap.vep.vcf" : "$tmp_dir/$input_name.vep.vcf" );
    warn "STATUS: Running VEP and writing to: $output_vcf\n";
    # Make sure we can find the VEP script
    my $vep_script = ( -s "$vep_path/vep" ? "$vep_path/vep" : "$vep_path/variant_effect_predictor.pl" );
    ( -s $vep_script ) or die "ERROR: Cannot find VEP script under: $vep_path\n";

    # Contruct VEP command using some default options and run it
    my $vep_cmd = "$perl_bin $vep_script --species $species --assembly $ncbi_build";
    $vep_cmd .= " --no_progress" unless( $verbose );
    $vep_cmd .= " --no_stats --buffer_size $buffer_size --sift b --ccds";
    $vep_cmd .= " --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical";
    $vep_cmd .= " --protein --biotype --uniprot --tsl --variant_class --shift_hgvs 1";
    $vep_cmd .= " --check_existing --total_length --allele_number --no_escape --xref_refseq";
    $vep_cmd .= " --failed 1 --vcf --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length";
    $vep_cmd .= " --dir $vep_data --fasta $ref_fasta --format vcf --input_file $input_vcf --output_file $output_vcf";
    $vep_cmd .= " --force_overwrite" if( $vep_overwrite );
    # Change options based on whether we are running in offline mode or not
    $vep_cmd .= ( $online ? " --database --host useastdb.ensembl.org" : " --offline --pubmed" );
    # VEP barks if --fork is set to 1. So don't use this argument unless it's >1
    $vep_cmd .= " --fork $vep_forks" if( $vep_forks > 1 );
    # Add --custom if requested at command line
    $vep_cmd .= " --custom $vep_custom" if ($vep_custom);
    # Add --config if requested at command line
    $vep_cmd .= " --config $vep_config" if ($vep_config);
    # Require allele match for co-located variants unless user-rejected or we're using a newer VEP
    $vep_cmd .= " --check_allele" unless( $any_allele or $vep_script =~ m/vep$/ );
    # Add --cache-version only if the user specifically asked for a version
    $vep_cmd .= " --cache_version $cache_version" if( $cache_version );
    # Add options that only work on human variants
    if( $species eq "homo_sapiens" ) {
        # Slight change in options if in offline mode, or if using the newer VEP
        $vep_cmd .= " --polyphen b" . ( $vep_script =~ m/vep$/ ? " --af" : " --gmaf" );
        $vep_cmd .= ( $vep_script =~ m/vep$/ ? " --af_1kg --af_esp --af_gnomad" : " --maf_1kg --maf_esp" ) unless( $online );
    }
    # Do not use the --regulatory option in situations where we know it will break
    $vep_cmd .= " --regulatory" unless( $species eq "canis_familiaris" or $online );

    warn "STATUS: Running this VEP command:  \n". wrap( "  ", "    ", $vep_cmd. "\n" ) if( $verbose );

    # Make sure it ran without error codes
    system( $vep_cmd ) == 0 or die "\nERROR: Failed to run the VEP annotator! Command: $vep_cmd\n";
    ( -s $output_vcf ) or warn "WARNING: VEP-annotated VCF file is missing or empty: $output_vcf\n";

    warn "STATUS: Finished with vep...\n" if( $verbose );
}

# Define default MAF Header (https://wiki.nci.nih.gov/x/eJaPAQ) with our vcf2maf additions
my @maf_header = qw(
    Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome Start_Position End_Position Strand
    Variant_Classification Variant_Type Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2
    dbSNP_RS dbSNP_Val_Status Tumor_Sample_Barcode Matched_Norm_Sample_Barcode
    Match_Norm_Seq_Allele1 Match_Norm_Seq_Allele2 Tumor_Validation_Allele1 Tumor_Validation_Allele2
    Match_Norm_Validation_Allele1 Match_Norm_Validation_Allele2 Verification_Status
    Validation_Status Mutation_Status Sequencing_Phase Sequence_Source Validation_Method Score
    BAM_File Sequencer Tumor_Sample_UUID Matched_Norm_Sample_UUID HGVSc HGVSp HGVSp_Short Transcript_ID
    Exon_Number t_depth t_ref_count t_alt_count n_depth n_ref_count n_alt_count all_effects
);

# Add extra annotation columns to the MAF in a consistent order
my @ann_cols = qw( Allele Gene Feature Feature_type Consequence cDNA_position CDS_position
    Protein_position Amino_acids Codons Existing_variation ALLELE_NUM DISTANCE STRAND_VEP SYMBOL
    SYMBOL_SOURCE HGNC_ID BIOTYPE CANONICAL CCDS ENSP SWISSPROT TREMBL UNIPARC RefSeq SIFT PolyPhen
    EXON INTRON DOMAINS AF AFR_AF AMR_AF ASN_AF EAS_AF EUR_AF SAS_AF AA_AF EA_AF CLIN_SIG SOMATIC
    PUBMED MOTIF_NAME MOTIF_POS HIGH_INF_POS MOTIF_SCORE_CHANGE IMPACT PICK VARIANT_CLASS TSL
    HGVS_OFFSET PHENO MINIMISED GENE_PHENO FILTER flanking_bps vcf_id vcf_qual gnomAD_AF gnomAD_AFR_AF
    gnomAD_AMR_AF gnomAD_ASJ_AF gnomAD_EAS_AF gnomAD_FIN_AF gnomAD_NFE_AF gnomAD_OTH_AF gnomAD_SAS_AF );

# push any requested custom VEP annotations from the CSQ/ANN section into @ann_cols
if ($retain_ann) {
    push @ann_cols, split(',',$retain_ann);
}
my @ann_cols_format; # To store the actual order of VEP data, that may differ between runs
push( @maf_header, @ann_cols );

# Add original VCF POS column header
push( @maf_header, "vcf_pos" );

# If the user has INFO fields they want to retain, create additional columns for those
my @addl_info_cols = ();
if( $retain_info or $remap_chain or $split_svs ) {
    # But let's not overwrite existing columns with the same name
    my %maf_cols = map{ my $c = lc; ( $c, 1 )} @maf_header;
    @addl_info_cols = grep{ my $c = lc; !$maf_cols{$c}} split( ",", $retain_info );
    # If a remap-chain was used, add a column to retain the original chr:pos:ref:alt
    push( @addl_info_cols, "REMAPPED_POS" ) if( $remap_chain );
    # If we had to split some SVs earlier, add some columns with some useful info about SVs
    push( @addl_info_cols, qw( Fusion Method Frame CONSENSUS )) if( $split_svs );
    push( @maf_header, @addl_info_cols );
}

# If the user has FORMAT fields they want to retain, create additional columns for those
my @addl_fmt_cols = ();
if( $retain_fmt ) {
    foreach my $fmt_tag ( split( ",", $retain_fmt )) {
        # Create 2 columns for tumor/normal, except if a column name is already used
        my ( $tc, $nc ) = ( "t_$fmt_tag", "n_$fmt_tag" );
        my %maf_cols = map{ my $c = lc; ( $c, 1 ) } @maf_header;
        if(!$maf_cols{lc($tc)}) { push (@addl_fmt_cols, $tc); }
        if(!$maf_cols{lc($nc)}) { push (@addl_fmt_cols, $nc); }
    }
    push( @maf_header, @addl_fmt_cols );
}

# Locate and load the file mapping ENSG IDs to Entrez IDs
my ( $script_dir ) = $0 =~ m/^(.*)\/vcf2maf/;
$script_dir = "." unless( $script_dir );

my $entrez_id_file = "$script_dir/data/ensg_to_entrez_id_map_ensembl_feb2014.tsv";
my %entrez_id_map = ();
if( -s $entrez_id_file ) {
    %entrez_id_map = map{chomp; split("\t")} `grep -hv ^# $entrez_id_file`;
}

# Parse through each variant in the annotated VCF, pull out CSQ/ANN from the INFO column, and choose
# one transcript per variant whose annotation will be used in the MAF
warn "STATUS: Parsing variants in annotated VCF...\n" if( $verbose );
my $maf_fh = IO::File->new( $output_maf, ">" ) or die "ERROR: Couldn't open --output-maf: $output_maf!\n";
$maf_fh->print( "#version 2.4\n" . join( "\t", @maf_header ), "\n" ); # Print MAF header
( -s $output_vcf ) or exit; # Warnings on this were printed earlier, but quit here, only after a blank MAF is created
my $annotated_vcf_fh = IO::File->new( $output_vcf ) or die "ERROR: Couldn't open annotated VCF: $output_vcf!\n";
my ( $vcf_tumor_idx, $vcf_normal_idx, %sv_pair );
while( my $line = $annotated_vcf_fh->getline ) {

    # Parse out the VEP CSQ/ANN format, which seems to differ between runs
    if( $line =~ m/^##INFO=<ID=(CSQ|ANN).*: '?([^"']+)["']/ ) {
        # Use this as the expected column order of VEP annotation, unless we already got it from CSQ
        @ann_cols_format = map{s/\s//g; $_} split( /\|/, $2 ) unless( @ann_cols_format and $1 eq "ANN" );
    }
    # Skip all other header lines
    next if( $line =~ m/^##/ );

    chomp( $line );
    my ( $chrom, $pos, $var_id, $ref, $alt, $var_qual, $filter, $info_line, $format_line, @rest ) = split( "\t", $line );

    # Set ID, QUAL, and FILTER to "." unless defined and non-empty
    $var_id = "." unless( defined $var_id and $var_id ne "" );
    $var_qual = "." unless( defined $var_qual and $var_qual ne "" );
    $filter = "." unless( defined $filter and $filter ne "" );

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
        next;
    }

    # Parse out the data in the info column, and store into a hash
    my %info = map {( $_, "1" )} grep { !m/=/ } split( /\;/, $info_line );
    map { my ( $key, $val ) = split( /=/, $_, 2 ); $info{$key} .= ( $info{$key} ? ",$val" : $val ) unless( $info{$key} and grep { /^$val$/ } split( ",", $info{$key} ))} grep { m/=/ } split( /\;/, $info_line );

    # By default, the variant allele is the first (usually the only) allele listed under ALT. If
    # there are >1 alleles in ALT, choose the first non-REF allele listed under tumor GT, that is
    # also not seen under normal GT. If tumor GT is undefined or ambiguous, choose the tumor allele
    # with the most supporting read depth, if available.
    my @alleles = ( $ref, split( /,/, $alt ));
    my $var_allele_idx = 1;

    # Parse out info from the normal genotype field
    my ( %nrm_info, @nrm_depths );
    if( defined $vcf_normal_idx ) {
        my @format_keys = split( /\:/, $format_line );
        my $idx = 0;
        %nrm_info = map {( $format_keys[$idx++], $_ )} split( /\:/, $rest[$vcf_normal_idx] );
    }

    # Parse out info from the tumor genotype field
    my ( %tum_info, @tum_depths );
    if( defined $vcf_tumor_idx ) {
        my @format_keys = split( /\:/, $format_line );
        my $idx = 0;
        %tum_info = map {( $format_keys[$idx++], $_ )} split( /\:/, $rest[$vcf_tumor_idx] );

        # If possible, parse the tumor genotype to identify the variant allele
        if( defined $tum_info{GT} and $tum_info{GT} ne "." and $tum_info{GT} ne "./." ) {
            my @tum_gt = split( /[\/|]/, $tum_info{GT} );
            # Default to the first non-REF allele seen in tumor GT
            ( $var_allele_idx ) = grep {$_ ne "0"} @tum_gt;
            # If possible, choose the first non-REF tumor allele that is also not in normal GT
            if( defined $nrm_info{GT} and $nrm_info{GT} ne "." and $nrm_info{GT} ne "./." ) {
                my %nrm_gt = map {( $_, 1 )} split( /[\/|]/, $nrm_info{GT} );
                ( $var_allele_idx ) = grep {$_ ne "0" and !$nrm_gt{$_}} @tum_gt;
            }
            # If GT was unhelpful, default to the first ALT allele and set GT to undefined
            if( !defined $var_allele_idx or $var_allele_idx !~ m/^\d+$/ or $var_allele_idx >= scalar( @alleles )) {
                $var_allele_idx = 1;
                $tum_info{GT} = "./.";
            }
        }

        # Standardize tumor AD and DP based on data in the genotype fields
        FixAlleleDepths( \@alleles, $var_allele_idx, \%tum_info );
        @tum_depths = split( ",", $tum_info{AD} );

        # If genotype is undefined, use the allele depths collected to choose the major variant allele
        unless( defined $tum_info{GT} and $tum_info{GT} ne '.' and $tum_info{GT} ne "./." ) {
            # The first depth listed belongs to the reference allele. Of the rest, find the largest
            for( my $i = 1; $i <= $#tum_depths; ++$i ) {
                $var_allele_idx = $i if( $tum_depths[$i] and $tum_depths[$i] > $tum_depths[$var_allele_idx] );
            }
            $tum_info{GT} = "./.";
            if( defined $tum_info{DP} and $tum_info{DP} ne '.' and $tum_info{DP} != 0 and defined $tum_depths[$var_allele_idx] ) {
                my $vaf = $tum_depths[$var_allele_idx] / $tum_info{DP};
                $tum_info{GT} = ( $vaf < $min_hom_vaf ? "0/1" : "1/1" );
            }
        }
    }

    # Set the variant allele to whatever we selected above
    my $var = $alleles[$var_allele_idx];

    # Standardize normal AD and DP based on data in the genotype fields
    if( defined $vcf_normal_idx ) {
        FixAlleleDepths( \@alleles, $var_allele_idx, \%nrm_info );
        @nrm_depths = split( ",", $nrm_info{AD} );
        $nrm_info{GT} = "./." unless( defined $nrm_info{GT} and $nrm_info{GT} ne '.' );
    }

    # Figure out the appropriate start/stop loci and variant type/allele to report in the MAF
    my $start = my $stop = my $var_type = my $inframe = "";
    my ( $ref_length, $var_length ) = ( length( $ref ), length( $var ));
    # Backup the VCF-style position and REF/ALT alleles, so we can use it later
    my ( $vcf_pos, $vcf_ref, $vcf_var ) = ( $pos, $ref, $var );
    # Remove any prefixed reference bps from all alleles, using "-" for simple indels
    while( $ref and $var and substr( $ref, 0, 1 ) eq substr( $var, 0, 1 ) and $ref ne $var ) {
        ( $ref, $var, @alleles ) = map{$_ = substr( $_, 1 ); ( $_ ? $_ : "-" )} ( $ref, $var, @alleles );
        --$ref_length; --$var_length; ++$pos;
    }
    # Handle SNPs, DNPs, TNPs, or anything larger (ONP)
    if( $ref_length == $var_length ) {
        ( $start, $stop ) = ( $pos, $pos + $var_length - 1 );
        my %np_type = qw( 1 SNP 2 DNP 3 TNP );
        $var_type = ( $var_length > 3 ? "ONP" : $np_type{$var_length} );
    }
    # Handle all indels, including those complex ones which contain substitutions
    elsif( $ref_length != $var_length ) {
        if( $ref_length < $var_length ) { # Handle insertions, and the special case for complex ones
            ( $start, $stop ) = (( $ref eq "-" ? $pos - 1 : $pos ), ( $ref eq "-" ? $pos : $pos + $ref_length - 1 ));
            $var_type = "INS";
        }
        else { # Handle deletions
            ( $start, $stop ) = ( $pos, $pos + $ref_length - 1 );
            $var_type = "DEL";
        }
        $inframe = ( abs( $ref_length - $var_length ) % 3 == 0 ? 1 : 0 );
    }

    my @all_effects; # A list of effects of this variant on all possible transcripts
    my $maf_effect; # A single effect per variant to report in the standard MAF columns
    my %maf_line = map{( $_, '' )} @maf_header; # Initialize MAF fields with blank strings

    # VEP provides a comma-delimited list of consequences, with pipe-delim details per consequence
    # It replaces ',' in details with '&'. We'll assume that all '&'s we see, were formerly commas
    # "Consequence" might list multiple effects on the same transcript e.g. missense,splice_region
    if( $info{CSQ} or $info{ANN} ) {

        my $ann_lines = ( $info{CSQ} ? $info{CSQ} : $info{ANN} );
        foreach my $ann_line ( split( /,/, $ann_lines )) {
            my $idx = 0;
            my %effect = map{s/\&/,/g; ( $ann_cols_format[$idx++], ( defined $_ ? $_ : '' ))} split( /\|/, $ann_line );

            # Remove transcript ID from HGVS codon/protein changes, to make it easier on the eye
            $effect{HGVSc} =~ s/^.*:// if( $effect{HGVSc} );
            $effect{HGVSp} =~ s/^.*:// if( $effect{HGVSp} );

            # Remove the prefixed HGVSc code in HGVSp, if found
            $effect{HGVSp} =~ s/^.*\((p\.\S+)\)/$1/ if( $effect{HGVSp} and $effect{HGVSp} =~ m/^c\./ );

            # If we find any snpEff fields, rename them to the corresponding VEP field names
            $effect{Consequence} = $effect{Annotation} if( $effect{Annotation} );
            $effect{IMPACT} = $effect{Annotation_Impact} if( $effect{Annotation_Impact} );
            $effect{SYMBOL} = $effect{Gene_Name} if( $effect{Gene_Name} );
            $effect{Gene} = $effect{Gene_ID} if( $effect{Gene_ID} );
            $effect{Feature_type} = $effect{Feature_Type} if( $effect{Feature_Type} );
            $effect{Feature} = $effect{Feature_ID} if( $effect{Feature_ID} );
            $effect{BIOTYPE} = $effect{Transcript_BioType} if( $effect{Transcript_BioType} );
            $effect{HGVSc} = $effect{'HGVS.c'} if( $effect{'HGVS.c'} );
            $effect{HGVSp} = $effect{'HGVS.p'} if( $effect{'HGVS.p'} );
            $effect{cDNA_position} = $effect{'cDNA.pos/cDNA.length'} if( $effect{'cDNA.pos/cDNA.length'} );
            $effect{CDS_position} = $effect{'CDS.pos/CDS.length'} if( $effect{'CDS.pos/CDS.length'} );
            $effect{Protein_position} = $effect{'AA.pos/AA.length'} if( $effect{'AA.pos/AA.length'} );
            $effect{DISTANCE} = $effect{Distance} if( $effect{Distance} );

            # Sort consequences by decreasing order of severity, and pick the most severe one
            $effect{Consequence} = join( ",", sort { GetEffectPriority($a) <=> GetEffectPriority($b) } split( ",", $effect{Consequence} ));
            ( $effect{One_Consequence} ) = split( ",", $effect{Consequence} );

            # When VEP fails to provide any value in Consequence, tag it as an intergenic variant
            $effect{One_Consequence} = "intergenic_variant" unless( $effect{Consequence} );

            # Create a shorter HGVS protein format using 1-letter codes
            if( $effect{HGVSp} ) {
                my $hgvs_p_short = $effect{HGVSp};
                while( $hgvs_p_short and my ( $find, $replace ) = each %aa3to1 ) {
                    eval "\$hgvs_p_short =~ s{$find}{$replace}g";
                }
                $effect{HGVSp_Short} = $hgvs_p_short;
            }

            # Fix HGVSp_Short, CDS_position, and Protein_position for splice acceptor/donor variants
            if( $effect{One_Consequence} =~ m/^(splice_acceptor_variant|splice_donor_variant)$/ ) {
                my ( $c_pos ) = $effect{HGVSc} =~ m/^c.(\d+)/;
                if( defined $c_pos ) {
                    $c_pos = 1 if( $c_pos < 1 ); # Handle negative cDNA positions used in 5' UTRs
                    my $p_pos = sprintf( "%.0f", ( $c_pos + $c_pos % 3 ) / 3 );
                    $effect{HGVSp_Short} = "p.X" . $p_pos . "_splice";
                    $effect{CDS_position} =~ s/^-(\/\d+)$/$c_pos$1/ if( $effect{CDS_position} );
                    $effect{Protein_position} =~ s/^-(\/\d+)$/$p_pos$1/ if( $effect{Protein_position} )
                }
            }

            # Fix HGVSp_Short for Silent mutations, so it mentions the amino-acid and position
            if( defined $effect{HGVSp_Short} and $effect{HGVSp_Short} eq "p.=" ) {
                my ( $p_pos ) = $effect{Protein_position} =~ m/^(\d+)(-\d+)?\/\d+$/;
                my $aa = $effect{Amino_acids};
                $effect{HGVSp_Short} = "p.$aa" . $p_pos . "=";
            }

            # Copy VEP data into MAF fields that don't share the same identifier
            $effect{Transcript_ID} = $effect{Feature};
            $effect{Exon_Number} = $effect{EXON};
            $effect{Hugo_Symbol} = ( $effect{SYMBOL} ? $effect{SYMBOL} : '' );

            # If AF columns from the older VEP are found, rename to the newer ones for consistency
            my %af_col = qw( GMAF AF AFR_MAF AFR_AF AMR_MAF AMR_AF ASN_MAF ASN_AF EAS_MAF EAS_AF
                EUR_MAF EUR_AF SAS_MAF SAS_AF AA_MAF AA_AF EA_MAF EA_AF );
            map { $effect{$af_col{$_}} = $effect{$_} if( defined $effect{$_} )} keys %af_col;

            # If VEP couldn't find this variant in dbSNP/COSMIC/etc., we'll say it's "novel"
            if( $effect{Existing_variation} ) {
                # ::NOTE:: If seen in a DB other than dbSNP, this field will remain blank
                $effect{dbSNP_RS} = join( ",", grep{m/^rs\d+$/} split( /,/, $effect{Existing_variation} ));
            }
            else {
                $effect{dbSNP_RS} = "novel";
            }

            # Transcript_Length isn't separately reported, but can be parsed out from cDNA_position
            ( $effect{Transcript_Length} ) = $effect{cDNA_position} =~ m/\/(\d+)$/ if( $effect{cDNA_position} );
            $effect{Transcript_Length} = 0 unless( defined $effect{Transcript_Length} );

            # Skip effects on other ALT alleles. If ALLELE_NUM is undefined (e.g. for INFO:SVTYPE), don't skip any
            push( @all_effects, \%effect ) unless( $effect{ALLELE_NUM} and $effect{ALLELE_NUM} != $var_allele_idx );
        }

        # Sort effects first by transcript biotype, then by severity, and then by longest transcript
        @all_effects = sort {
            GetBiotypePriority( $a->{BIOTYPE} ) <=> GetBiotypePriority( $b->{BIOTYPE} ) ||
            GetEffectPriority( $a->{One_Consequence} ) <=> GetEffectPriority( $b->{One_Consequence} ) ||
            $b->{Transcript_Length} <=> $a->{Transcript_Length}
        } @all_effects;

        # Find the highest priority effect with a gene symbol i.e. the worst affected gene
        my ( $effect_with_gene_name ) = grep { $_->{SYMBOL} } @all_effects;
        my $maf_gene = $effect_with_gene_name->{SYMBOL} if( $effect_with_gene_name );

        # If that gene has a user-preferred isoform, report the effect on that isoform
        ( $maf_effect ) = grep { $_->{SYMBOL} and $_->{SYMBOL} eq $maf_gene and $_->{Transcript_ID} and $custom_enst{$_->{Transcript_ID}} } @all_effects;

        # If that gene has no user-preferred isoform, then use the VEP-preferred (canonical) isoform
        ( $maf_effect ) = grep { $_->{SYMBOL} and $_->{SYMBOL} eq $maf_gene and $_->{CANONICAL} and $_->{CANONICAL} eq "YES" } @all_effects unless( $maf_effect );

        # If that gene has no VEP-preferred isoform either, then choose the worst affected user-preferred isoform with a gene symbol
        ( $maf_effect ) = grep { $_->{SYMBOL} and $_->{Transcript_ID} and $custom_enst{$_->{Transcript_ID}} } @all_effects unless( $maf_effect );

        # If none of the isoforms are user-preferred, then choose the worst affected VEP-preferred isoform with a gene symbol
        ( $maf_effect ) = grep { $_->{SYMBOL} and $_->{CANONICAL} and $_->{CANONICAL} eq "YES" } @all_effects unless( $maf_effect );

        # If we still have nothing selected, then just report the worst effect
        $maf_effect = $all_effects[0] unless( $maf_effect );
    }

    # Construct the MAF columns from the $maf_effect hash
    %maf_line = map{( $_, ( $maf_effect->{$_} ? $maf_effect->{$_} : '' ))} @maf_header;
    $maf_line{Hugo_Symbol} = $maf_effect->{Transcript_ID} unless( $maf_effect->{Hugo_Symbol} );
    $maf_line{Hugo_Symbol} = 'Unknown' unless( $maf_effect->{Transcript_ID} );
    $maf_line{Entrez_Gene_Id} = ( defined $maf_effect->{Gene} && defined $entrez_id_map{$maf_effect->{Gene}} ? $entrez_id_map{$maf_effect->{Gene}} : "0" );
    $maf_line{Center} = $maf_center;
    $maf_line{NCBI_Build} = $ncbi_build;
    $maf_line{Chromosome} = $chrom;
    $maf_line{Start_Position} = $start;
    $maf_line{End_Position} = $stop;
    $maf_line{Strand} = '+'; # Per MAF definition, only the positive strand is an accepted value
    $maf_line{STRAND_VEP} = $maf_effect->{STRAND}; # Renamed to avoid mixup with "Strand" above
    $maf_line{Variant_Classification} = GetVariantClassification( $maf_effect->{One_Consequence}, $var_type, $inframe );
    $maf_line{Variant_Type} = $var_type;
    $maf_line{Reference_Allele} = $ref;
    # ::NOTE:: If tumor genotype is unavailable, then we'll assume it's ref/var heterozygous
    $maf_line{Tumor_Seq_Allele1} = $ref;
    $maf_line{Tumor_Seq_Allele2} = $var;
    if( defined $tum_info{GT} and $tum_info{GT} ne "." and $tum_info{GT} ne "./." ) {
        # ::NOTE:: MAF only supports biallelic sites. Tumor_Seq_Allele2 must always be the $var
        # picked earlier. For Tumor_Seq_Allele1, pick the first non-var allele in GT (usually $ref)
        my ( $idx1, $idx2 ) = split( /[\/|]/, $tum_info{GT} );
        # If GT was monoploid, then $idx2 will be undefined, and we should set it equal to $idx1
        $idx2 = $idx1 unless( defined $idx2 );
        $maf_line{Tumor_Seq_Allele1} = ( $alleles[$idx1] ne $var ? $alleles[$idx1] : $alleles[$idx2] );
    }
    # ::NOTE:: If normal genotype is unavailable, then we'll assume it's ref/ref homozygous
    $maf_line{Match_Norm_Seq_Allele1} = $ref;
    $maf_line{Match_Norm_Seq_Allele2} = $ref;
    if( defined $nrm_info{GT} and $nrm_info{GT} ne "." and $nrm_info{GT} ne "./." ) {
        # ::NOTE:: MAF only supports biallelic sites. So choose the first two alleles listed in GT
        my ( $idx1, $idx2 ) = split( /[\/|]/, $nrm_info{GT} );
        # If GT was monoploid, then $idx2 will be undefined, and we should set it equal to $idx1
        $idx2 = $idx1 unless( defined $idx2 );
        $maf_line{Match_Norm_Seq_Allele1} = $alleles[$idx1];
        $maf_line{Match_Norm_Seq_Allele2} = $alleles[$idx2];
    }
    $maf_line{Tumor_Sample_Barcode} = $tumor_id;
    $maf_line{Matched_Norm_Sample_Barcode} = $normal_id;
    $maf_line{t_depth} = $tum_info{DP} if( defined $tum_info{DP} and $tum_info{DP} ne "." );
    ( $maf_line{t_ref_count}, $maf_line{t_alt_count} ) = @tum_depths[0,$var_allele_idx] if( @tum_depths );
    $maf_line{n_depth} = $nrm_info{DP} if( defined $nrm_info{DP} and $nrm_info{DP} ne "." );
    ( $maf_line{n_ref_count}, $maf_line{n_alt_count} ) = @nrm_depths[0,$var_allele_idx] if( @nrm_depths );

    # Create a semicolon delimited list summarizing the prioritized effects in @all_effects
    $maf_line{all_effects} = "";
    foreach my $effect ( @all_effects ) {
        my $gene_name = $effect->{Hugo_Symbol};
        my $effect_type = $effect->{One_Consequence};
        my $protein_change = ( $effect->{HGVSp} ? $effect->{HGVSp} : '' );
        my $transcript_id = ( $effect->{Transcript_ID} ? $effect->{Transcript_ID} : '' );
        my $refseq_ids = ( $effect->{RefSeq} ? $effect->{RefSeq} : '' );
        $maf_line{all_effects} .= "$gene_name,$effect_type,$protein_change,$transcript_id,$refseq_ids;" if( $effect_type and $transcript_id );
    }

    # Copy FILTER from input VCF, and tag calls with high allele freq in any gnomAD subpopulation
    my $subpop_count = 0;
    foreach my $subpop ( qw( AFR AMR ASJ EAS FIN NFE SAS )) {
        if( $maf_line{"gnomAD_$subpop\_AF"} ) {
            my ( $subpop_af ) = split( "/", $maf_line{"gnomAD_$subpop\_AF"} );
            $subpop_count++ if( $subpop_af > $max_subpop_af );
        }
    }
    # Remove existing common_variant tags from input, so it's redefined by this new criteria
    $filter = join( ";", grep{ $_ ne "common_variant" } split( /,|;/, $filter ));
    if( $subpop_count > 0 ) {
        $filter = (( $filter eq "PASS" or $filter eq "." or !$filter ) ? "common_variant" : "$filter;common_variant" );
    }
    $maf_line{FILTER} = $filter;

    # Also add the reference allele flanking bps that we generated earlier with samtools
    my $region = "$chrom:" . ( $vcf_pos - 1 ) . "-" . ( $vcf_pos + length( $vcf_ref ));
    $maf_line{flanking_bps} = $flanking_bps{$region};

    # Add ID and QUAL from the input VCF into respective MAF columns
    $maf_line{vcf_id} = $var_id;
    $maf_line{vcf_qual} = $var_qual;

    # Add original VCF POS column
    $maf_line{vcf_pos} = $vcf_pos;

    # If there are additional INFO data to add, then add those
    foreach my $info_col ( @addl_info_cols ) {
        $maf_line{$info_col} = ( defined $info{$info_col} ? $info{$info_col} : "" );
    }

    # If there are additional FORMAT data to add, then add those
    foreach my $fmt_col ( @addl_fmt_cols ) {
        my $fmt_key = $fmt_col;
        if ( $fmt_key =~ /^t_/ ) { $fmt_key =~ s/^t_//; $maf_line{$fmt_col} = ( defined $tum_info{$fmt_key} ? $tum_info{$fmt_key} : "" ); }
        if ( $fmt_key =~ /^n_/ ) { $fmt_key =~ s/^n_//; $maf_line{$fmt_col} = ( defined $nrm_info{$fmt_key} ? $nrm_info{$fmt_key} : "" ); }
    }

    # If this is an SV, pair up gene names from separate lines to backfill the Fusion column later
    if( $split_svs and $var=~m/^<BND|DEL|DUP|INV>$/ ) {
        my $sv_key = "$var_id-$tumor_id";
        if( $sv_pair{$sv_key} ) {
            $sv_pair{$sv_key} = $sv_pair{$sv_key} . "-" . $maf_line{Hugo_Symbol} . " fusion";
        }
        else {
            $sv_pair{$sv_key} = $maf_line{Hugo_Symbol};
        }
    }

    # At this point, we've generated all we can about this variant, so write it to the MAF
    $maf_fh->print( join( "\t", map{( defined $maf_line{$_} ? $maf_line{$_} : "" )} @maf_header ) . "\n" );
}
$maf_fh->close;
$annotated_vcf_fh->close;

# If the MAF lists SVs, backfill the Fusion column with gene-pair names
warn "STATUS: For any SVs, backfilling Fusion column with gene-pair names...\n" if( $verbose );
if( $split_svs ) {
    my $output_name = substr( $output_maf, rindex( $output_maf, "/" ) + 1 );
    $output_name =~ s/(\.maf)*$//;
    my $tmp_output_maf = "$tmp_dir/$output_name.tmp.maf";

    my $in_maf_fh = IO::File->new( $output_maf ) or die "ERROR: Couldn't open: $output_maf!\n";
    my $out_maf_fh = IO::File->new( $tmp_output_maf, ">" ) or die "ERROR: Couldn't open: $tmp_output_maf!\n";
    my ( $tid_idx, $fusion_idx, $var_id_idx ) = ( 0, 0, 0 );
    while( my $line = $in_maf_fh->getline ) {
        chomp( $line );
        if( $line =~ m/^#/ ) {
            $out_maf_fh->print( "$line\n" ); # Copy comments unchanged
        }
        elsif( $line =~ m/^Hugo_Symbol/ ) {
            # Copy the header unchanged, after figuring out necessary column indexes
            foreach( split( /\t/, $line )) { last if( $_ eq "Tumor_Sample_Barcode" ); ++$tid_idx; }
            foreach( split( /\t/, $line )) { last if( $_ eq "Fusion" ); ++$fusion_idx; }
            foreach( split( /\t/, $line )) { last if( $_ eq "vcf_id" ); ++$var_id_idx; }
            $out_maf_fh->print( "$line\n" ); # Copy header unchanged
        }
        else {
            # Write the gene-pair name into the Fusion column if it was backfilled earlier
            my @cols = split( /\t/, $line, -1 );
            my $sv_key = $cols[$var_id_idx] . "-" . $cols[$tid_idx];
            $cols[$fusion_idx] = $sv_pair{$sv_key} if( $sv_pair{$sv_key} );
            $out_maf_fh->print( join( "\t", @cols ) . "\n" );
        }
    }
    $out_maf_fh->close;
    $in_maf_fh->close;

    move( $tmp_output_maf, $output_maf );
}

warn "STATUS: Finished! Check results in $output_maf\n" if( $verbose );

# Converts Sequence Ontology variant types to MAF variant classifications
sub GetVariantClassification {
    my ( $effect, $var_type, $inframe ) = @_;
    return "Targeted_Region" if( not defined $effect or not $effect ); # In case VEP was skipped
    return "Splice_Site" if( $effect =~ /^(splice_acceptor_variant|splice_donor_variant|transcript_ablation|exon_loss_variant)$/ );
    return "Nonsense_Mutation" if( $effect eq 'stop_gained' );
    return "Frame_Shift_Del" if(( $effect eq 'frameshift_variant' or ( $effect eq 'protein_altering_variant' and !$inframe )) and $var_type eq 'DEL' );
    return "Frame_Shift_Ins" if(( $effect eq 'frameshift_variant' or ( $effect eq 'protein_altering_variant' and !$inframe )) and $var_type eq 'INS' );
    return "Nonstop_Mutation" if( $effect eq 'stop_lost' );
    return "Translation_Start_Site" if( $effect =~ /^(initiator_codon_variant|start_lost)$/ );
    return "In_Frame_Ins" if( $effect =~ /inframe_insertion$/ or ( $effect eq 'protein_altering_variant' and $inframe and $var_type eq 'INS' ));
    return "In_Frame_Del" if( $effect =~ /inframe_deletion$/ or ( $effect eq 'protein_altering_variant' and $inframe and $var_type eq 'DEL' ));
    return "Missense_Mutation" if( $effect =~ /^(missense_variant|coding_sequence_variant|conservative_missense_variant|rare_amino_acid_variant)$/ );
    return "Intron" if ( $effect =~ /^(transcript_amplification|intron_variant|INTRAGENIC|intragenic_variant)$/ );
    return "Splice_Region" if( $effect eq 'splice_region_variant' );
    return "Silent" if( $effect =~ /^(incomplete_terminal_codon_variant|synonymous_variant|stop_retained_variant|NMD_transcript_variant)$/ );
    return "RNA" if( $effect =~ /^(mature_miRNA_variant|exon_variant|non_coding_exon_variant|non_coding_transcript_exon_variant|non_coding_transcript_variant|nc_transcript_variant)$/ );
    return "5'UTR" if( $effect =~ /^(5_prime_UTR_variant|5_prime_UTR_premature_start_codon_gain_variant)$/ );
    return "3'UTR" if( $effect eq '3_prime_UTR_variant' );
    return "IGR" if( $effect =~ /^(TF_binding_site_variant|regulatory_region_variant|regulatory_region|intergenic_variant|intergenic_region)$/ );
    return "5'Flank" if( $effect eq 'upstream_gene_variant' );
    return "3'Flank" if ( $effect eq 'downstream_gene_variant' );

    # Annotate everything else simply as a targeted region
    # TFBS_ablation, TFBS_amplification,regulatory_region_ablation, regulatory_region_amplification,
    # feature_elongation, feature_truncation
    return "Targeted_Region";
}

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
            my %acgt_depths = map{( defined $fmt_info{$_.'U'} ? ( $_, $fmt_info{$_.'U'} ) : ( $_, "" ))} qw( A C G T );
            my @deepest = sort {$acgt_depths{$b} <=> $acgt_depths{$a}} keys %acgt_depths;
            ( $alleles[1] ) = ( $deepest[0] ne $alleles[0] ? $deepest[0] : $deepest[1] );
        }
        @depths = map{( defined $fmt_info{$_.'U'} ? $fmt_info{$_.'U'} : "" )} @alleles;
    }
    # Handle VCF Indel lines by Strelka, where variant allele depth is in TIR and reference allele depth in TAR
    elsif( !defined $fmt_info{AD} and defined $fmt_info{TIR} and defined $fmt_info{TAR}) {
        @depths = map{""} @alleles;
        $depths[0] = ( split /,/, $fmt_info{TAR} )[0];
        $depths[$var_allele_idx] = ( split /,/, $fmt_info{TIR} )[0];
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
    # Handle VCF lines where REF/ALT allele counts must be extracted from DP4
    elsif( !defined $fmt_info{AD} and defined $fmt_info{DP4} and scalar( split( /,/, $fmt_info{DP4} )) == 4 ) {
        # Reference allele depth and depths for any other ALT alleles must be left undefined
        @depths = map{""} @alleles;
        # DP4 is usually a comma-delimited list for ref-forward, ref-reverse, alt-forward and alt-reverse read counts
        my @count = split( /,/, $fmt_info{DP4} );
        $depths[0] = $count[0] + $count[1];
        $depths[$var_allele_idx] = $count[2] + $count[3];
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

 vcf2maf.pl - Convert a VCF into a MAF by mapping each variant to only one of all possible gene isoforms

=head1 SYNOPSIS

 perl vcf2maf.pl --help
 perl vcf2maf.pl --input-vcf WD4086.vcf --output-maf WD4086.maf --tumor-id WD4086 --normal-id NB4086

=head1 OPTIONS

 --input-vcf      Path to input file in VCF format
 --output-maf     Path to output MAF file
 --tmp-dir        Folder to retain intermediate VCFs after runtime [Default: Folder containing input VCF]
 --tumor-id       Tumor_Sample_Barcode to report in the MAF [TUMOR]
 --normal-id      Matched_Norm_Sample_Barcode to report in the MAF [NORMAL]
 --vcf-tumor-id   Tumor sample ID used in VCF's genotype columns [--tumor-id]
 --vcf-normal-id  Matched normal ID used in VCF's genotype columns [--normal-id]
 --custom-enst    List of custom ENST IDs that override canonical selection
 --vep-path       Folder containing the vep script [~/miniconda3/bin]
 --vep-data       VEP's base cache/plugin directory [~/.vep]
 --vep-forks      Number of forked processes to use when running VEP [4]
 --vep-custom     String to pass into VEP's --custom option []
 --vep-config     Config file to pass into VEP's --config option []
 --vep-overwrite  Allow VEP to overwrite output VCF if it exists
 --buffer-size    Number of variants VEP loads at a time; Reduce this for low memory systems [5000]
 --any-allele     When reporting co-located variants, allow mismatched variant alleles too
 --inhibit-vep    Skip running VEP, but extract VEP annotation in VCF if found
 --online         Use useastdb.ensembl.org instead of local cache (supports only GRCh38 VCFs listing <100 events)
 --ref-fasta      Reference FASTA file [~/.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa.gz]
 --max-subpop-af  Add FILTER tag common_variant if gnomAD reports any subpopulation AFs greater than this [0.0004]
 --species        Ensembl-friendly name of species (e.g. mus_musculus for mouse) [homo_sapiens]
 --ncbi-build     NCBI reference assembly of variants MAF (e.g. GRCm38 for mouse) [GRCh37]
 --cache-version  Version of offline cache to use with VEP (e.g. 75, 91, 102) [Default: Installed version]
 --maf-center     Variant calling center to report in MAF [.]
 --retain-info    Comma-delimited names of INFO fields to retain as extra columns in MAF []
 --retain-fmt     Comma-delimited names of FORMAT fields to retain as extra columns in MAF []
 --retain-ann     Comma-delimited names of annotations (within the VEP CSQ/ANN) to retain as extra columns in MAF []
 --min-hom-vaf    If GT undefined in VCF, minimum allele fraction to call a variant homozygous [0.7]
 --remap-chain    Chain file to remap variants to a different assembly before running VEP
 --verbose        Print more things to log progress
 --help           Print a brief help message and quit
 --man            Print the detailed manual

=head1 DESCRIPTION

To convert a VCF into a MAF, each variant must be mapped to only one of all possible gene transcripts/isoforms that it might affect. This selection of a single effect per variant, is often subjective. So this project is an attempt to make the selection criteria smarter, reproducible, and more configurable.

This script uses VEP, a variant annotator that maps effects of a variant on all possible genes and transcripts. For more info, see the README.

=head2 Relevant links:

 Homepage: https://github.com/ckandoth/vcf2maf
 VCF format: http://samtools.github.io/hts-specs/
 MAF format: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format
 VEP: http://ensembl.org/info/docs/tools/vep/index.html
 VEP annotated VCF format: http://ensembl.org/info/docs/tools/vep/vep_formats.html#vcfout

=head1 AUTHORS

 Cyriac Kandoth (ckandoth@gmail.com)
 Shweta Chavan (chavan.shweta@gmail.com)
 Zuojian Tang (zuojian.tang@gmail.com)

=head1 LICENSE

 Apache-2.0 | Apache License, Version 2.0 | https://www.apache.org/licenses/LICENSE-2.0

=cut
