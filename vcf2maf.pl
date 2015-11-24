#!/usr/bin/env perl

# vcf2maf - Convert a VCF into a MAF by mapping each variant to only one of all possible gene isoforms

use strict;
use warnings;
use IO::File;
use Getopt::Long qw( GetOptions );
use Pod::Usage qw( pod2usage );
use Config;

# Set any default paths and constants
my ( $tumor_id, $normal_id ) = ( "TUMOR", "NORMAL" );
my ( $vep_path, $vep_data, $vep_forks, $ref_fasta ) = ( "$ENV{HOME}/vep", "$ENV{HOME}/.vep", 4, "$ENV{HOME}/.vep/homo_sapiens/82_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz" );
my ( $species, $ncbi_build, $maf_center, $min_hom_vaf ) = ( "homo_sapiens", "GRCh37", ".", 0.7 );
my $perl_bin = $Config{perlpath};

# Hash to convert 3-letter amino-acid codes to their 1-letter codes
my %aa3to1 = qw( Ala A Arg R Asn N Asp D Asx B Cys C Glu E Gln Q Glx Z Gly G His H Ile I Leu L
    Lys K Met M Phe F Pro P Ser S Thr T Trp W Tyr Y Val V Xxx X Ter * );

# Prioritize Sequence Ontology terms in order of severity, as estimated by Ensembl:
# http://useast.ensembl.org/info/genome/variation/predicted_data.html#consequences
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
        'inframe_insertion' => 5, # An inframe non synonymous variant that inserts bases into the coding sequence
        'inframe_deletion' => 5, # An inframe non synonymous variant that deletes bases from the coding sequence
        'missense_variant' => 6, # A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved
        'conservative_missense_variant' => 6, # A sequence variant whereby at least one base of a codon is changed resulting in a codon that encodes for a different but similar amino acid. These variants may or may not be deleterious
        'rare_amino_acid_variant' => 6, # A sequence variant whereby at least one base of a codon encoding a rare amino acid is changed, resulting in a different encoded amino acid
        'transcript_amplification' => 7, # A feature amplification of a region containing a transcript
        'stop_retained_variant' => 8, # A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
        'synonymous_variant' => 8, # A sequence variant where there is no resulting change to the encoded amino acid
        'splice_region_variant' => 9, # A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron
        'incomplete_terminal_codon_variant' => 10, # A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed
        'protein_altering_variant' => 11, # A sequence variant which is predicted to change the protein encoded in the coding sequence
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
    $effectPriority{$effect} or die "ERROR: Unrecognized effect \"$effect\". Please update your hashes!\n";
    return $effectPriority{$effect};
}

# Prioritize the transcript biotypes that variants are annotated to, based on disease significance:
# All possible biotypes are defined here: http://www.gencodegenes.org/gencode_biotypes.html
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
        'sRNA' => 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'scaRNA' => 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'rRNA' => 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'lincRNA' => 3, # Long, intervening noncoding (linc) RNAs, that can be found in evolutionarily conserved, intergenic regions
        'known_ncrna' => 4,
        'vaultRNA' => 4, # Short non coding RNA genes that form part of the vault ribonucleoprotein complex
        'macro_lncRNA' => 4, # unspliced lncRNAs that are several kb in size
        'Mt_tRNA' => 4, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'Mt_rRNA' => 4, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'antisense' => 5, # Has transcripts that overlap the genomic span (i.e. exon or introns) of a protein-coding locus on the opposite strand
        'sense_intronic' => 5, # Long non-coding transcript in introns of a coding gene that does not overlap any exons
        'sense_overlapping' => 5, # Long non-coding transcript that contains a coding gene in its intron on the same strand
        '3prime_overlapping_ncrna' => 5, # Transcripts where ditag and/or published experimental data strongly supports the existence of short non-coding transcripts transcribed from the 3'UTR
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
        '' => 9
    );
    $biotype_priority{$biotype} or die "ERROR: Unrecognized biotype \"$biotype\". Please update your hashes!\n";
    return $biotype_priority{$biotype};
}

# Check for missing or crappy arguments
unless( @ARGV and $ARGV[0] =~ m/^-/ ) {
    pod2usage( -verbose => 0, -message => "$0: Missing or invalid arguments!\n", -exitval => 2 );
}

# Parse options and print usage if there is a syntax error, or if usage was explicitly requested
my ( $man, $help ) = ( 0, 0 );
my ( $input_vcf, $output_maf, $custom_enst_file );
my ( $vcf_tumor_id, $vcf_normal_id );
GetOptions(
    'help!' => \$help,
    'man!' => \$man,
    'input-vcf=s' => \$input_vcf,
    'output-maf=s' => \$output_maf,
    'tumor-id=s' => \$tumor_id,
    'normal-id=s' => \$normal_id,
    'vcf-tumor-id=s' => \$vcf_tumor_id,
    'vcf-normal-id=s' => \$vcf_normal_id,
    'custom-enst=s' => \$custom_enst_file,
    'vep-path=s' => \$vep_path,
    'vep-data=s' => \$vep_data,
    'vep-forks=s' => \$vep_forks,
    'ref-fasta=s' => \$ref_fasta,
    'species=s' => \$species,
    'ncbi-build=s' => \$ncbi_build,
    'maf-center=s' => \$maf_center,
    'min-hom-vaf=s' => \$min_hom_vaf
) or pod2usage( -verbose => 1, -input => \*DATA, -exitval => 2 );
pod2usage( -verbose => 1, -input => \*DATA, -exitval => 0 ) if( $help );
pod2usage( -verbose => 2, -input => \*DATA, -exitval => 0 ) if( $man );

# Unless specified, assume that the VCF uses the same sample IDs that the MAF will contain
$vcf_tumor_id = $tumor_id unless( $vcf_tumor_id );
$vcf_normal_id = $normal_id unless( $vcf_normal_id );

# Load up the custom isoform overrides if provided:
my %custom_enst;
if( $custom_enst_file ) {
	( -s $custom_enst_file ) or die "ERROR: The provided custom ENST file is missing or empty!\nPath: $custom_enst_file\n";
    %custom_enst = map{chomp; ( $_, 1 )}`grep -v ^# $custom_enst_file | cut -f1`;
}

# Annotate variants in given VCF to all possible transcripts
my $output_vcf;
if( $input_vcf ) {
    ( -s $input_vcf ) or die "ERROR: Provided VCF file is missing or empty!\nPath: $input_vcf\n";
    ( $input_vcf !~ m/\.(gz|bz2|bcf)$/ ) or die "ERROR: Compressed or binary VCFs are not supported\n";

    $output_vcf = $input_vcf;
    $output_vcf =~ s/(\.vcf)*$/.vep.vcf/;

    # Skip running VEP if an annotated VCF already exists
    unless( -s $output_vcf ) {
        warn "STATUS: Running VEP and writing to: $output_vcf\n";
        # Make sure we can find the VEP script and the reference FASTA
        ( -s "$vep_path/variant_effect_predictor.pl" ) or die "ERROR: Cannot find VEP script variant_effect_predictor.pl in path: $vep_path\n";
        ( -s $ref_fasta ) or die "ERROR: Reference FASTA not found: $ref_fasta\n";

        # Contruct VEP command using some default options and run it
        my $vep_cmd = "$perl_bin $vep_path/variant_effect_predictor.pl --species $species --assembly $ncbi_build --offline --no_progress --no_stats --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --regulatory --canonical --protein --biotype --uniprot --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --check_alleles --check_ref --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --minimal --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --dir $vep_data --fasta $ref_fasta --input_file $input_vcf --output_file $output_vcf";
        $vep_cmd .= " --fork $vep_forks" if( $vep_forks > 1 ); # VEP barks if it's set to 1
        # Add options that only work on human variants
        $vep_cmd .= " --polyphen b --gmaf --maf_1kg --maf_esp" if( $species eq "homo_sapiens" );
        # Add options that only work on human variants mapped to the GRCh37 reference genome
        $vep_cmd .= " --plugin ExAC,$vep_data/ExAC.r0.3.sites.minus_somatic.vcf.gz" if( $species eq "homo_sapiens" and $ncbi_build eq "GRCh37" );

        # Make sure it ran without error codes
        system( $vep_cmd ) == 0 or die "\nERROR: Failed to run the VEP annotator!\nCommand: $vep_cmd\n";
        ( -s $output_vcf ) or warn "WARNING: VEP-annotated VCF file is missing or empty!\nPath: $output_vcf\n";
    }
}
else {
    die "ERROR: Please specify an input file: input-vcf. STDIN is not supported.\n";
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
    Protein_position Amino_acids Codons Existing_variation ALLELE_NUM DISTANCE STRAND SYMBOL
    SYMBOL_SOURCE HGNC_ID BIOTYPE CANONICAL CCDS ENSP SWISSPROT TREMBL UNIPARC RefSeq SIFT PolyPhen
    EXON INTRON DOMAINS GMAF AFR_MAF AMR_MAF ASN_MAF EAS_MAF EUR_MAF SAS_MAF AA_MAF EA_MAF CLIN_SIG
    SOMATIC PUBMED MOTIF_NAME MOTIF_POS HIGH_INF_POS MOTIF_SCORE_CHANGE IMPACT PICK VARIANT_CLASS
    TSL HGVS_OFFSET PHENO MINIMISED ExAC_AF ExAC_AF_AFR ExAC_AF_AMR ExAC_AF_EAS ExAC_AF_FIN
    ExAC_AF_NFE ExAC_AF_OTH ExAC_AF_SAS GENE_PHENO FILTER );
my @ann_cols_format; # To store the actual order of VEP data, that may differ between runs
push( @maf_header, @ann_cols );

# Parse through each variant in the annotated VCF, pull out ANN/CSQ from the INFO column, and choose
# one transcript per variant whose annotation will be used in the MAF
my $maf_fh = *STDOUT; # Use STDOUT if an output MAF file was not defined
$maf_fh = IO::File->new( $output_maf, ">" ) or die "ERROR: Couldn't open output file: $output_maf!\n" if( $output_maf );
$maf_fh->print( "#version 2.4\n" . join( "\t", @maf_header ), "\n" ); # Print MAF header
( -s $output_vcf ) or exit; # Warnings on this were printed earlier, but quit here, only after a blank MAF is created
my $vcf_fh = IO::File->new( $output_vcf ) or die "ERROR: Couldn't open annotated VCF: $output_vcf!\n";
my ( $vcf_tumor_idx, $vcf_normal_idx );
while( my $line = $vcf_fh->getline ) {

    # Parse out the VEP ANN/CSQ format, which seems to differ between runs
    if( $line =~ m/^##INFO=<ID=(ANN|CSQ).*Format: (\S+)">$/ ) {
        @ann_cols_format = split( /\|/, $2 );
    }
    # Skip all other header lines
    next if( $line =~ m/^##/ );

    chomp( $line );
    my ( $chrom, $pos, $ids, $ref, $alt, $qual, $filter, $info_line, $format_line, @rest ) = split( /\t/, $line );

    # Set QUAL and FILTER to "." unless defined and non-empty
    $qual = "." unless( defined $qual and $qual ne "" );
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
    my %info = map {( m/=/ ? ( split( /=/, $_, 2 )) : ( $_, 1 ))} split( /\;/, $info_line );

    # By default, the variant allele is the first (usually the only) allele listed under ALT
    # If there are multiple ALT alleles, choose the allele specified in the tumor GT field
    # If tumor GT is undefined or ambiguous, choose the one with the most supporting read depth
    my @alleles = ( $ref, split( /,/, $alt ));
    my $var_allele_idx = 1;

    # Parse out info from the tumor genotype field
    my ( %tum_info, @tum_depths );
    if( defined $vcf_tumor_idx ) {
        my @format_keys = split( /\:/, $format_line );
        my $idx = 0;
        %tum_info = map {( $format_keys[$idx++], $_ )} split( /\:/, $rest[$vcf_tumor_idx] );

        # If possible, parse the tumor genotype to identify the variant allele
        if( defined $tum_info{GT} and $tum_info{GT} ne "." and $tum_info{GT} ne "./." ) {
            my @genotype = split( /[\/|]/, $tum_info{GT} );
            # In case of polyploid calls, choose the first non-REF allele, if any
            ( $var_allele_idx ) = grep {$_ ne "0"} @genotype;
            $var_allele_idx = 1 unless( defined $var_allele_idx and $var_allele_idx =~ m/^\d+$/ );
        }

        # If AD is defined, then parse out all REF/ALT allele depths, or whatever is in it
        if( defined $tum_info{AD} and $tum_info{AD} ne "." ) {
            @tum_depths = map{( m/^\d+$/ ? $_ : "" )}split( /,/, $tum_info{AD} );
        }

        # Handle VarScan VCF lines where AD contains only 1 depth, and REF allele depth is in RD
        if( scalar( @tum_depths ) == 1 and defined $tum_info{RD} ) {
            @tum_depths = map{""} @alleles;
            $tum_depths[0] = $tum_info{RD};
            $tum_depths[$var_allele_idx] = $tum_info{AD};
        }
        # Handle SomaticSniper VCF lines, where allele depths must be extracted from BCOUNT
        elsif( defined $tum_info{BCOUNT} ) {
            my %b_idx = ( A=>0, C=>1, G=>2, T=>3 );
            my @bcount = split( /,/, $tum_info{BCOUNT} );
            @tum_depths = map{(( defined $b_idx{$_} and defined $bcount[$b_idx{$_}] ) ? $bcount[$b_idx{$_}] : "" )} @alleles;
        }
        # Handle VCF SNV lines by Strelka, where allele depths are in AU:CU:GU:TU
        elsif( scalar( grep{defined $tum_info{$_}} qw/AU CU GU TU/ ) == 4 ) {
            # Strelka allele depths come in tiers 1,2. We'll use tier1 cuz it's stricter, and DP already is
            map{( $tum_info{$_.'U'} ) = split( ",", $tum_info{$_.'U'} )} qw( A C G T );

            # If the ALT allele is just N, then set it to the allele with the highest non-ref readcount
            if( scalar( @alleles ) == 2 and $alleles[1] eq "N" ) {
                my %acgt_depths = map{( defined $tum_info{$_.'U'} ? ( $_, $tum_info{$_.'U'} ) : ( $_, "" ))} qw( A C G T );
                my @deepest = sort {$acgt_depths{$b} <=> $acgt_depths{$a}} keys %acgt_depths;
                ( $alleles[1] ) = ( $deepest[0] ne $ref ? $deepest[0] : $deepest[1] );
            }
            @tum_depths = map{( defined $tum_info{$_.'U'} ? $tum_info{$_.'U'} : "" )} @alleles;
        }
        # Handle VCF Indel lines by Strelka, where variant allele depth is in TIR
        elsif( $tum_info{TIR} ) {
            # Reference allele depth is not provided by Strelka for indels, so we have to skip it
            @tum_depths = ( "", ( split /,/, $tum_info{TIR} )[0] );
        }
        # Handle VCF lines from the Ion Torrent Suite where ALT depths are in AO and REF depths are in RO
        elsif( defined $tum_info{AO} and defined $tum_info{RO} ) {
            @tum_depths = ( $tum_info{RO}, map{( m/^\d+$/ ? $_ : "" )}split( /,/, $tum_info{AO} ));
        }
        # Handle VCF lines with ALT allele-frac in FA, which needs to be multiplied by DP to get AD
        elsif( defined $tum_info{FA} and defined $tum_info{DP} and $tum_info{DP} ne '.' ) {
            # Reference allele depth and depths for any other ALT alleles must be left undefined
            @tum_depths = map{""} @alleles;
            $tum_depths[$var_allele_idx] = sprintf( "%.0f", $tum_info{FA} * $tum_info{DP} );
        }
        # Handle VCF lines where AD contains only 1 value, that we can assume is the variant allele
        elsif( defined $tum_info{AD} and @tum_depths and scalar( @tum_depths ) == 1 ) {
            # Reference allele depth and depths for any other ALT alleles must be left undefined
            @tum_depths = map{""} @alleles;
            $tum_depths[$var_allele_idx] = $tum_info{AD};
        }
        # Handle VCF lines from mpileup/bcftools where DV contains the ALT allele depth
        elsif( defined $tum_info{DV} and defined $tum_info{DP} ) {
            # Reference allele depth and depths for any other ALT alleles must be left undefined
            @tum_depths = map{""} @alleles;
            $tum_depths[$var_allele_idx] = $tum_info{DV};
        }
        # For all other lines where #depths is not equal to #alleles, blank out the depths
        elsif( @tum_depths and $#tum_depths != $#alleles ) {
            warn "WARNING: Unusual AD format for alleles $ref,$alt in $format_line = " . $rest[$vcf_tumor_idx] . "\n";
            @tum_depths = map{""} @alleles;
        }

        # Sanity check that REF/ALT allele depths are lower than the total depth
        if( defined $tum_info{DP} and $tum_info{DP} ne '.' and (( $tum_depths[0] and $tum_depths[0] > $tum_info{DP} ) or
            ( $tum_depths[$var_allele_idx] and $tum_depths[$var_allele_idx] > $tum_info{DP} ) or
            ( $tum_depths[0] and $tum_depths[$var_allele_idx] and $tum_depths[0] + $tum_depths[$var_allele_idx] > $tum_info{DP} ))) {
            warn "WARNING: AD of alleles $ref,$alt exceed DP. Setting DP to sum of allele depths in $format_line = " . $rest[$vcf_tumor_idx] . "\n";
            $tum_info{DP} = 0;
            map{$tum_info{DP} += $_ if($_ and $_ ne '.')} @tum_depths;
        }

        # If we have REF/ALT allele depths, but no DP, then set DP equal to the sum of all ADs
        if(( defined $tum_depths[0] and defined $tum_depths[$var_allele_idx] ) and ( !defined $tum_info{DP} or $tum_info{DP} eq '.' )) {
            $tum_info{DP} = 0;
            map{$tum_info{DP} += $_ if($_ and $_ ne '.')} @tum_depths;
        }

        # If genotype is undefined, use the allele depths collected to choose the major variant allele
        unless( defined $tum_info{GT} and $tum_info{GT} ne '.' and $tum_info{GT} ne "./." ) {
            # The first depth listed belongs to the reference allele. Of the rest, find the largest
            for( my $i = 1; $i <= $#tum_depths; ++$i ) {
                $var_allele_idx = $i if( $tum_depths[$i] and $tum_depths[$i] > $tum_depths[$var_allele_idx] );
            }
        }
    }

    # Set the variant allele to whatever we selected above
    my $var = $alleles[$var_allele_idx];

    # Unless GT is defined, assign het/hom status based on a simple variant allele fraction cutoff
    unless( defined $tum_info{GT} and $tum_info{GT} ne "." and $tum_info{GT} ne "./." ) {
        $tum_info{GT} = "./.";
        if( defined $tum_info{DP} and $tum_info{DP} != 0 and defined $tum_depths[$var_allele_idx] ) {
            my $vaf = $tum_depths[$var_allele_idx] / $tum_info{DP};
            $tum_info{GT} = ( $vaf < $min_hom_vaf ? "0/1" : "1/1" );
        }
    }

    # Same as above, parse out info from the normal genotype field
    my ( %nrm_info, @nrm_depths );
    if( defined $vcf_normal_idx ) {
        my @format_keys = split( /\:/, $format_line );
        my $idx = 0;
        %nrm_info = map {( $format_keys[$idx++], $_ )} split( /\:/, $rest[$vcf_normal_idx] );

        # If AD is defined, then parse out all REF/ALT allele depths, or whatever is in it
        if( defined $nrm_info{AD} and $nrm_info{AD} ne "." ) {
            @nrm_depths = map{( m/^\d+$/ ? $_ : "" )}split( /,/, $nrm_info{AD} );
        }

        # Handle VCF lines by VarScan where REF allele depth is stored separately in an RD tag
        if( scalar( @nrm_depths ) == 1 and defined $nrm_info{RD} ) {
            @nrm_depths = map{""} @alleles;
            $nrm_depths[0] = $nrm_info{RD};
            $nrm_depths[$var_allele_idx] = $nrm_info{AD};
        }
        # Handle VCF lines by SomaticSniper, where allele depths must be extracted from BCOUNT
        elsif( defined $nrm_info{BCOUNT} ) {
            my %b_idx = ( A=>0, C=>1, G=>2, T=>3 );
            my @bcount = split( /,/, $nrm_info{BCOUNT} );
            @nrm_depths = map{(( defined $b_idx{$_} and defined $bcount[$b_idx{$_}] ) ? $bcount[$b_idx{$_}] : "" )} @alleles;
        }
        # Handle VCF SNV lines by Strelka, where allele depths are in AU:CU:GU:TU
        elsif( scalar( grep{defined $nrm_info{$_}} qw/AU CU GU TU/ ) == 4 ) {
            # Strelka allele depths come in tiers 1,2. We'll use tier1 cuz it's stricter, and DP already is
            map{( $nrm_info{$_.'U'} ) = split( ",", $nrm_info{$_.'U'} )} qw( A C G T );

            # If the ALT allele is just N, then set it to the allele with the highest non-ref readcount
            if( scalar( @alleles ) == 2 and $alleles[1] eq "N" ) {
                my %acgt_depths = map{( defined $nrm_info{$_.'U'} ? ( $_, $nrm_info{$_.'U'} ) : ( $_, "" ))} qw( A C G T );
                my @deepest = sort {$acgt_depths{$b} <=> $acgt_depths{$a}} keys %acgt_depths;
                ( $alleles[1] ) = ( $deepest[0] ne $ref ? $deepest[0] : $deepest[1] );
            }
            @nrm_depths = map{( defined $nrm_info{$_.'U'} ? $nrm_info{$_.'U'} : "" )} @alleles;
        }
        # Handle VCF Indel lines by Strelka, where variant allele depth is in TIR
        elsif( $nrm_info{TIR} ) {
            # Reference depth is not explicitly defined by Strelka for indels, so we have to skip it
            @nrm_depths = ( "", ( split /,/, $nrm_info{TIR} )[0] );
        }
        # Handle VCF lines from the Ion Torrent Suite where ALT depths are in AO and REF depths are in RO
        elsif( defined $nrm_info{AO} and defined $nrm_info{RO} ) {
            @nrm_depths = ( $nrm_info{RO}, map{( m/^\d+$/ ? $_ : "" )}split( /,/, $nrm_info{AO} ));
        }
        # Handle VCF lines with ALT allele-frac in FA, which needs to be multiplied by DP to get AD
        elsif( defined $nrm_info{FA} and defined $nrm_info{DP} and $nrm_info{DP} ne '.' ) {
            # Reference allele depth and depths for any other ALT alleles must be left undefined
            @nrm_depths = map{""} @alleles;
            $nrm_depths[$var_allele_idx] = sprintf( "%.0f", $nrm_info{FA} * $nrm_info{DP} );
        }
        # Handle VCF lines where AD contains only 1 value, that we can assume is the variant allele
        elsif( defined $nrm_info{AD} and @nrm_depths and scalar( @nrm_depths ) == 1 ) {
            # Reference allele depth and depths for any other ALT alleles must be left undefined
            @nrm_depths = map{""} @alleles;
            $nrm_depths[$var_allele_idx] = $nrm_info{AD};
        }
        # Handle VCF lines from mpileup/bcftools where DV contains the ALT allele depth
        elsif( defined $nrm_info{DV} and defined $nrm_info{DP} ) {
            # Reference allele depth and depths for any other ALT alleles must be left undefined
            @nrm_depths = map{""} @alleles;
            $nrm_depths[$var_allele_idx] = $nrm_info{DV};
        }
        # For all other lines where #depths is not equal to #alleles, blank out the depths
        elsif( @nrm_depths and $#nrm_depths != $#alleles ) {
            warn "WARNING: Unusual AD format for alleles $ref,$alt in $format_line = " . $rest[$vcf_normal_idx] . "\n";
            @nrm_depths = map{""} @alleles;
        }

        # Sanity check that REF/ALT allele depths are lower than the total depth
        if( defined $nrm_info{DP} and $nrm_info{DP} ne '.' and (( $nrm_depths[0] and $nrm_depths[0] > $nrm_info{DP} ) or
            ( $nrm_depths[$var_allele_idx] and $nrm_depths[$var_allele_idx] > $nrm_info{DP} ) or
            ( $nrm_depths[0] and $nrm_depths[$var_allele_idx] and $nrm_depths[0] + $nrm_depths[$var_allele_idx] > $nrm_info{DP} ))) {
            warn "WARNING: AD of alleles $ref,$alt exceed DP. Setting DP to sum of allele depths in $format_line = " . $rest[$vcf_normal_idx] . "\n";
            $nrm_info{DP} = 0;
            map{$nrm_info{DP} += $_ if($_ and $_ ne '.')} @nrm_depths;
        }

        # If we have REF/ALT allele depths, but no DP, then set DP equal to the sum of all ADs
        if(( defined $nrm_depths[0] and defined $nrm_depths[$var_allele_idx] ) and ( !defined $nrm_info{DP} or $nrm_info{DP} eq '.' )) {
            $nrm_info{DP} = 0;
            map{$nrm_info{DP} += $_ if($_ and $_ ne '.')} @nrm_depths;
        }
    }

    # Figure out the appropriate start/stop loci and variant type/allele to report in the MAF
    my $start = my $stop = my $var_type = "";
    my ( $ref_length, $var_length ) = ( length( $ref ), length( $var ));
    # Remove any prefixed reference bps from all alleles, using "-" for simple indels
    while( substr( $ref, 0, 1 ) eq substr( $var, 0, 1 )) {
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
            ( $start, $stop ) = ( $pos - 1, ( $ref eq "-" ? $pos : $pos + $ref_length - 1 ));
            $var_type = "INS";
        }
        else { # Handle deletions
            ( $start, $stop ) = ( $pos, $pos + $ref_length - 1 );
            $var_type = "DEL";
        }
    }

    my @all_effects; # A list of effects of this variant on all possible transcripts
    my $maf_effect; # A single effect per variant to report in the standard MAF columns
    my %maf_line = map{( $_, '' )} @maf_header; # Initialize MAF fields with blank strings

    # VEP provides a comma-delimited list of consequences, with pipe-delim details per consequence
    # It replaces ',' in details with '&'. We'll assume that all '&'s we see, were formerly commas
    # "Consequence" might list multiple effects on the same transcript e.g. missense,splice_region
    if( $info{ANN} or $info{CSQ} ) {

        my $ann_lines = ( $info{ANN} ? $info{ANN} : $info{CSQ} );
        foreach my $ann_line ( split( /,/, $ann_lines )) {
            my $idx = 0;
            my %effect = map{s/\&/,/g; ( $ann_cols_format[$idx++], ( defined $_ ? $_ : '' ))} split( /\|/, $ann_line );

            # Remove transcript ID from HGVS codon/protein changes, to make it easier on the eye
            $effect{HGVSc} =~ s/^.*:// if( $effect{HGVSc} );
            $effect{HGVSp} =~ s/^.*:// if( $effect{HGVSp} );

            # Remove the prefixed HGVSc code in HGVSp, if found
            $effect{HGVSp} =~ s/^.*\((p\.\S+)\)/$1/ if( $effect{HGVSp} and $effect{HGVSp} =~ m/^c\./ );

            # If there are several consequences listed for a transcript, choose the most severe one
            ( $effect{Consequence} ) = sort { GetEffectPriority($a) <=> GetEffectPriority($b) } split( /,/, $effect{Consequence} );

            # When VEP fails to provide any value in Consequence, tag it as an intergenic variant
            $effect{Consequence} = "intergenic_variant" unless( $effect{Consequence} );

            # Create a shorter HGVS protein format using 1-letter codes
            if( $effect{HGVSp} ) {
                my $hgvs_p_short = $effect{HGVSp};
                while( $hgvs_p_short and my ( $find, $replace ) = each %aa3to1 ) {
                    eval "\$hgvs_p_short =~ s{$find}{$replace}g";
                }
                $effect{HGVSp_Short} = $hgvs_p_short;
            }

            # Fix HGVSp_Short, CDS_position, and Protein_position for splice acceptor/donor variants
            if( $effect{Consequence} =~ m/^(splice_acceptor_variant|splice_donor_variant)$/ ) {
                my ( $c_pos ) = $effect{HGVSc} =~ m/^c.(\d+)/;
                if( defined $c_pos ) {
                    $c_pos = 1 if( $c_pos < 1 ); # Handle negative cDNA positions used in 5' UTRs
                    my $p_pos = sprintf( "%.0f", ( $c_pos + $c_pos % 3 ) / 3 );
                    $effect{HGVSp_Short} = "p.X" . $p_pos . "_splice";
                    $effect{CDS_position} =~ s/^-(\/\d+)$/$c_pos$1/;
                    $effect{Protein_position} =~ s/^-(\/\d+)$/$p_pos$1/;
                }
            }

            # Fix HGVSp_Short for Silent mutations, so it mentions the amino-acid and position
            if( $effect{Consequence} eq "synonymous_variant" ) {
                my ( $p_pos ) = $effect{Protein_position} =~ m/^(\d+)(-\d+)?\/\d+$/;
                my $aa = $effect{Amino_acids};
                $effect{HGVSp_Short} = "p.$aa" . $p_pos . $aa;
            }

            # Copy VEP data into MAF fields that don't share the same identifier
            $effect{Transcript_ID} = $effect{Feature};
            $effect{Exon_Number} = $effect{EXON};
            $effect{Hugo_Symbol} = ( $effect{SYMBOL} ? $effect{SYMBOL} : '' );
            # If VEP couldn't find this variant in dbSNP/COSMIC/etc., we'll say it's "novel"
            if( $effect{Existing_variation} ) {
                # ::NOTE:: If seen in a DB other than dbSNP, this field will remain blank
                $effect{dbSNP_RS} = join( ",", grep{m/^rs\d+$/} split( /,/, $effect{Existing_variation} ));
            }
            else {
                $effect{dbSNP_RS} = "novel";
            }

            # Transcript_Length isn't separately reported, but can be parsed out from cDNA_position
            ( $effect{Transcript_Length} ) = $effect{cDNA_position} =~ m/\/(\d+)$/;
            $effect{Transcript_Length} = 0 unless( defined $effect{Transcript_Length} );

            # Skip effects on other ALT alleles. If ALLELE_NUM is undefined (e.g. for INFO:SVTYPE), don't skip any
            push( @all_effects, \%effect ) unless( $effect{ALLELE_NUM} and $effect{ALLELE_NUM} != $var_allele_idx );
        }

        # Sort effects first by transcript biotype, then by severity, and then by longest transcript
        @all_effects = sort {
            GetBiotypePriority( $a->{BIOTYPE} ) <=> GetBiotypePriority( $b->{BIOTYPE} ) ||
            GetEffectPriority( $a->{Consequence} ) <=> GetEffectPriority( $b->{Consequence} ) ||
            $b->{Transcript_Length} <=> $a->{Transcript_Length}
        } @all_effects;

        # Find the highest priority effect with a gene symbol (usually the first one)
        my ( $effect_with_gene_name ) = grep { $_->{SYMBOL} } @all_effects;
        my $maf_gene = $effect_with_gene_name->{SYMBOL} if( $effect_with_gene_name );

        # If the gene has user-defined custom isoform overrides, choose that instead
        ( $maf_effect ) = grep { $_->{SYMBOL} and $_->{SYMBOL} eq $maf_gene and $_->{Transcript_ID} and $custom_enst{$_->{Transcript_ID}} } @all_effects;

        # Find the effect on the canonical transcript of that highest priority gene
        ( $maf_effect ) = grep { $_->{SYMBOL} and $_->{SYMBOL} eq $maf_gene and $_->{CANONICAL} and $_->{CANONICAL} eq "YES" } @all_effects unless( $maf_effect );

        # If that gene has no canonical transcript tagged, choose the highest priority canonical effect on any gene
        ( $maf_effect ) = grep { $_->{CANONICAL} and $_->{CANONICAL} eq "YES" } @all_effects unless( $maf_effect );

        # If none of the effects are tagged as canonical, then just report the top priority effect
        $maf_effect = $all_effects[0] unless( $maf_effect );
    }

    # Construct the MAF columns from the $maf_effect hash, and print to output
    %maf_line = map{( $_, ( $maf_effect->{$_} ? $maf_effect->{$_} : '' ))} @maf_header;
    $maf_line{Hugo_Symbol} = $maf_effect->{Transcript_ID} unless( $maf_effect->{Hugo_Symbol} );
    $maf_line{Entrez_Gene_Id} = '0';
    $maf_line{Center} = $maf_center;
    $maf_line{NCBI_Build} = $ncbi_build;
    $maf_line{Chromosome} = $chrom;
    $maf_line{Start_Position} = $start;
    $maf_line{End_Position} = $stop;
    $maf_line{Strand} = '+'; # Per MAF definition, only the positive strand is an accepted value
    my $so_effect = ( $maf_effect->{Effect} ? $maf_effect->{Effect} : $maf_effect->{Consequence} );
    $maf_line{Variant_Classification} = GetVariantClassification( $so_effect, $var_type);
    $maf_line{Variant_Type} = $var_type;
    $maf_line{Reference_Allele} = $ref;
    # ::NOTE:: If tumor genotype is unavailable, then we'll assume it's ref/var heterozygous
    $maf_line{Tumor_Seq_Allele1} = $ref;
    $maf_line{Tumor_Seq_Allele2} = $var;
    if( defined $tum_info{GT} and $tum_info{GT} ne "." and $tum_info{GT} ne "./." ) {
        # ::NOTE:: MAF only supports biallelic sites. Tumor_Seq_Allele2 must always be the $var
        # picked earlier. For Tumor_Seq_Allele1, pick the first non-var allele in GT (usually $ref)
        my ( $idx1, $idx2 ) = split( /[\/|]/, $tum_info{GT} );
        # If GT was monoploid, then $idx2 will undefined, and we should set it equal to $idx1
        $idx2 = $idx1 unless( defined $idx2 );
        $maf_line{Tumor_Seq_Allele1} = ( $alleles[$idx1] ne $var ? $alleles[$idx1] : $alleles[$idx2] );
    }
    # ::NOTE:: If normal genotype is unavailable, then we'll assume it's ref/ref homozygous
    $maf_line{Match_Norm_Seq_Allele1} = $ref;
    $maf_line{Match_Norm_Seq_Allele2} = $ref;
    if( defined $nrm_info{GT} and $nrm_info{GT} ne "." and $nrm_info{GT} ne "./." ) {
        # ::NOTE:: MAF only supports biallelic sites. So choose the first two alleles listed in GT
        my ( $idx1, $idx2 ) = split( /[\/|]/, $nrm_info{GT} );
        # If GT was monoploid, then $idx2 will undefined, and we should set it equal to $idx1
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
        my $effect_type = ( $effect->{Effect} ? $effect->{Effect} : $effect->{Consequence} );
        my $protein_change = ( $effect->{HGVSp} ? $effect->{HGVSp} : '' );
        my $transcript_id = ( $effect->{Transcript_ID} ? $effect->{Transcript_ID} : '' );
        my $refseq_ids = ( $effect->{RefSeq} ? $effect->{RefSeq} : '' );
        $maf_line{all_effects} .= "$gene_name,$effect_type,$protein_change,$transcript_id,$refseq_ids;" if( $gene_name and $effect_type and $transcript_id );
    }

    # Apply FILTER from the input VCF, and also tag calls with high minor allele frequency in
    # any EVS/1000G/ExAC subpopulation, unless ClinVar says pathogenic, risk_factor, or protective
    my ( $max_subpop_maf, $subpop_count ) = ( 0.0004, 0 );
    foreach my $subpop ( qw( GMAF AFR_MAF AMR_MAF ASN_MAF EAS_MAF EUR_MAF SAS_MAF AA_MAF EA_MAF )) {
        my ( undef, $freq ) = split( /:|,/, $maf_line{$subpop} );
        $subpop_count++ if( defined $freq and $freq > $max_subpop_maf );
    }
    foreach my $subpop ( qw( ExAC_AF ExAC_AF_AFR ExAC_AF_AMR ExAC_AF_EAS ExAC_AF_FIN ExAC_AF_NFE ExAC_AF_OTH ExAC_AF_SAS )) {
        $subpop_count++ if( $maf_line{$subpop} ne "" and $maf_line{$subpop} > $max_subpop_maf );
    }
    if( $subpop_count > 0 and $maf_line{CLIN_SIG} !~ /pathogenic|risk_factor|protective/ ) {
        $filter = (( $filter eq "PASS" or $filter eq "." ) ? "common_variant" : "$filter,common_variant" );
    }
    $maf_line{FILTER} = $filter;

    # At this point, we've generated all we can about this variant, so write it to the MAF
    $maf_fh->print( join( "\t", map{ $maf_line{$_} } @maf_header ) . "\n" );
}
$maf_fh->close if( $output_maf );
$vcf_fh->close;

# Converts Sequence Ontology variant types to MAF variant classifications
sub GetVariantClassification {
    my ( $effect, $var_type ) = @_;
    return "Splice_Site" if( $effect =~ /^(splice_acceptor_variant|splice_donor_variant|transcript_ablation)$/ );
    return "Nonsense_Mutation" if( $effect eq 'stop_gained' );
    return "Frame_Shift_Del" if ( $effect eq 'frameshift_variant' and $var_type eq 'DEL' );
    return "Frame_Shift_Ins" if( $effect eq 'frameshift_variant' and $var_type eq 'INS' );
    return "Nonstop_Mutation" if( $effect eq 'stop_lost' );
    return "Translation_Start_Site" if( $effect =~ /^(initiator_codon_variant|start_lost)$/ );
    return "In_Frame_Ins" if( $effect eq 'inframe_insertion' );
    return "In_Frame_Del" if( $effect eq 'inframe_deletion' );
    return "Missense_Mutation" if( $effect =~ /^(missense_variant|coding_sequence_variant|conservative_missense_variant|rare_amino_acid_variant)$/ );
    return "Intron" if ( $effect =~ /^(transcript_amplification|splice_region_variant|intron_variant|INTRAGENIC|intragenic_variant)$/ );
    return "Silent" if( $effect =~ /^(incomplete_terminal_codon_variant|synonymous_variant|stop_retained_variant|NMD_transcript_variant)$/ );
    return "RNA" if( $effect =~ /^(mature_miRNA_variant|non_coding_exon_variant|non_coding_transcript_exon_variant|non_coding_transcript_variant|nc_transcript_variant)$/ );
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

__DATA__

=head1 NAME

 vcf2maf.pl - Convert a VCF into a MAF by mapping each variant to only one of all possible gene isoforms

=head1 SYNOPSIS

 perl vcf2maf.pl --help
 perl vcf2maf.pl --input-vcf WD4086.vcf --output-maf WD4086.maf --tumor-id WD4086 --normal-id NB4086

=head1 OPTIONS

 --input-vcf      Path to input file in VCF format
 --output-maf     Path to output MAF file [Default: STDOUT]
 --tumor-id       Tumor_Sample_Barcode to report in the MAF [TUMOR]
 --normal-id      Matched_Norm_Sample_Barcode to report in the MAF [NORMAL]
 --vcf-tumor-id   Tumor sample ID used in VCF's genotype columns [--tumor-id]
 --vcf-normal-id  Matched normal ID used in VCF's genotype columns [--normal-id]
 --custom-enst    List of custom ENST IDs that override canonical selection
 --vep-path       Folder containing variant_effect_predictor.pl [~/vep]
 --vep-data       VEP's base cache/plugin directory [~/.vep]
 --vep-forks      Number of forked processes to use when running VEP [4]
 --ref-fasta      Reference FASTA file [~/.vep/homo_sapiens/82_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz]
 --species        Ensembl-friendly name of species (e.g. mus_musculus for mouse) [homo_sapiens]
 --ncbi-build     NCBI reference assembly of variants MAF (e.g. GRCm38 for mouse) [GRCh37]
 --maf-center     Variant calling center to report in MAF [.]
 --min-hom-vaf    If GT undefined in VCF, minimum allele fraction to call a variant homozygous [0.7]
 --help           Print a brief help message and quit
 --man            Print the detailed manual

=head1 DESCRIPTION

To convert a VCF into a MAF, each variant must be mapped to only one of all possible gene transcripts/isoforms that it might affect. This selection of a single effect per variant, is often subjective. So this project is an attempt to make the selection criteria smarter, reproducible, and more configurable.

This script needs VEP, a variant annotator that maps effects of a variant on all possible genes and transcripts. For more info, see the README.

=head2 Relevant links:

 Homepage: https://github.com/ckandoth/vcf2maf
 VCF format: http://samtools.github.io/hts-specs/
 MAF format: https://wiki.nci.nih.gov/x/eJaPAQ
 VEP: http://ensembl.org/info/docs/tools/vep/index.html
 VEP annotated VCF format: http://ensembl.org/info/docs/tools/vep/vep_formats.html#vcfout

=head1 AUTHORS

 Cyriac Kandoth (ckandoth@gmail.com)

=head1 LICENSE

 Apache-2.0 | Apache License, Version 2.0 | https://www.apache.org/licenses/LICENSE-2.0

=cut
