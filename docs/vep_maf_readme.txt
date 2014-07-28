The first 34 columns are standard MAF format, and described here:
https://wiki.nci.nih.gov/x/eJaPAQ

The subsequent 10 columns are relevant to most analyses.

35. HGVSc - the coding sequence of the variant in HGVS recommended format
36. HGVSp - the protein sequence of the variant in HGVS recommended format
37. Transcript_ID - transcript onto which the consequence of the variant has been mapped
38. Exon_Number - the exon number (out of total number)
39. t_depth - read depth across this locus in tumor BAM
40. t_ref_count - read depth supporting the reference allele in tumor BAM
41. t_alt_count - read depth supporting the variant allele in tumor BAM
42. n_depth - read depth across this locus in normal BAM
43. n_ref_count - read depth supporting the reference allele in normal BAM
44. n_alt_count - read depth supporting the variant allele in normal BAM

The next column is relevant to analyses that consider the effect of the variant on all alternate
isoforms of the gene, or on non-coding/regulatory transcripts. The effects are sorted first by
transcript biotype priority, then by effect severity, and finally by decreasing order of transcript
length. Each effect in the list is in the format [SYMBOL,Consequence,HGVSp,Transcript_ID].

45. all_effects - a semicolon delimited list of all possible variant effects, sorted by priority

All remaining columns are straight out of Ensembl's VEP annotator, as described here:
http://useast.ensembl.org/info/docs/tools/vep/vep_formats.html#output

46. Allele - the variant allele used to calculate the consequence
47. Gene - stable Ensembl ID of affected gene
48. Feature - stable Ensembl ID of feature
49. Feature_type - type of feature. Currently one of Transcript, RegulatoryFeature, MotifFeature
50. Consequence - consequence type of this variation
51. cDNA_position - relative position of base pair in cDNA sequence
52. CDS_position - relative position of base pair in coding sequence
53. Protein_position - relative position of amino acid in protein
54. Amino_acids - only given if the variation affects the protein-coding sequence
55. Codons - the alternative codons with the variant base in upper case
56. Existing_variation - known identifier of existing variation
57. AA_MAF - minor allele and frequency of existing variant in NHLBI-ESP African American population
58. EA_MAF - minor allele and frequency of existing variant in NHLBI-ESP European American population
59. ALLELE_NUM - allele number from input; 0 is reference, 1 is first alternate etc
60. RefSeq - RefSeq identifier for this transcript
61. EXON - the exon number (out of total number)
62. INTRON - the intron number (out of total number)
63. MOTIF_NAME - the source and identifier of a transcription factor binding profile aligned at this position
64. MOTIF_POS - the relative position of the variation in the aligned TFBP
65. HIGH_INF_POS - a flag indicating if the variant falls in a high information position of a transcription factor binding profile (TFBP)
66. MOTIF_SCORE_CHANGE - the difference in motif score of the reference and variant sequences for the TFBP
67. DISTANCE - shortest distance from variant to transcript
68. STRAND - the DNA strand (1 or -1) on which the transcript/feature lies
69. CLIN_SIG - clinical significance of variant from dbSNP
70. CANONICAL - a flag indicating if the transcript is denoted as the canonical transcript for this gene
71. SYMBOL - the gene symbol
72. SYMBOL_SOURCE - the source of the gene symbol
73. SIFT - the SIFT prediction and/or score, with both given as prediction (score)
74. PolyPhen - the PolyPhen prediction and/or score
75. GMAF - minor allele and frequency of existing variation in 1000 Genomes Phase 1
76. BIOTYPE - biotype of transcript
77. ENSP - the Ensembl protein identifier of the affected transcript
78. DOMAINS - the source and identifer of any overlapping protein domains
79. CCDS - the CCDS identifer for this transcript, where applicable
80. HGVSc - the coding sequence of the variant in HGVS recommended format
81. HGVSp - the protein sequence of the variant in HGVS recommended format
82. AFR_MAF - minor allele and frequency of existing variation in 1000 Genomes Phase 1 combined African population
83. AMR_MAF - minor allele and frequency of existing variation in 1000 Genomes Phase 1 combined American population
84. ASN_MAF - minor allele and frequency of existing variation in 1000 Genomes Phase 1 combined Asian population
85. EUR_MAF - minor allele and frequency of existing variation in 1000 Genomes Phase 1 combined European population
86. PUBMED - pubmed ID(s) of publications that cite existing variant
