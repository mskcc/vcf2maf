The first 34 columns are standard MAF format, and described here:
https://wiki.nci.nih.gov/x/eJaPAQ

The subsequent 10 columns are relevant to most analyses:

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

All remaining columns are straight out of Ensembl's VEP annotator, as described here:
http://useast.ensembl.org/info/docs/tools/vep/vep_formats.html#output

45. Allele - the variant allele used to calculate the consequence
46. Gene - stable Ensembl ID of affected gene
47. Feature - stable Ensembl ID of feature
48. Feature_type - type of feature. Currently one of Transcript, RegulatoryFeature, MotifFeature
49. Consequence - consequence type of this variation
50. cDNA_position - relative position of base pair in cDNA sequence
51. CDS_position - relative position of base pair in coding sequence
52. Protein_position - relative position of amino acid in protein
53. Amino_acids - only given if the variation affects the protein-coding sequence
54. Codons - the alternative codons with the variant base in upper case
55. Existing_variation - known identifier of existing variation
56. AA_MAF - minor allele and frequency of existing variant in NHLBI-ESP African American population
57. EA_MAF - minor allele and frequency of existing variant in NHLBI-ESP European American population
58. ALLELE_NUM - allele number from input; 0 is reference, 1 is first alternate etc
59. RefSeq - RefSeq identifier for this transcript
60. EXON - the exon number (out of total number)
61. INTRON - the intron number (out of total number)
62. MOTIF_NAME - the source and identifier of a transcription factor binding profile aligned at this position
63. MOTIF_POS - the relative position of the variation in the aligned TFBP
64. HIGH_INF_POS - a flag indicating if the variant falls in a high information position of a transcription factor binding profile (TFBP)
65. MOTIF_SCORE_CHANGE - the difference in motif score of the reference and variant sequences for the TFBP
66. DISTANCE - shortest distance from variant to transcript
67. STRAND - the DNA strand (1 or -1) on which the transcript/feature lies
68. CLIN_SIG - clinical significance of variant from dbSNP
69. CANONICAL - a flag indicating if the transcript is denoted as the canonical transcript for this gene
70. SYMBOL - the gene symbol
71. SYMBOL_SOURCE - the source of the gene symbol
72. SIFT - the SIFT prediction and/or score, with both given as prediction (score)
73. PolyPhen - the PolyPhen prediction and/or score
74. GMAF - minor allele and frequency of existing variation in 1000 Genomes Phase 1
75. BIOTYPE - biotype of transcript
76. ENSP - the Ensembl protein identifier of the affected transcript
77. DOMAINS - the source and identifer of any overlapping protein domains
78. CCDS - the CCDS identifer for this transcript, where applicable
79. HGVSc - the coding sequence of the variant in HGVS recommended format
80. HGVSp - the protein sequence of the variant in HGVS recommended format
81. AFR_MAF - minor allele and frequency of existing variation in 1000 Genomes Phase 1 combined African population
82. AMR_MAF - minor allele and frequency of existing variation in 1000 Genomes Phase 1 combined American population
83. ASN_MAF - minor allele and frequency of existing variation in 1000 Genomes Phase 1 combined Asian population
84. EUR_MAF - minor allele and frequency of existing variation in 1000 Genomes Phase 1 combined European population
85. PUBMED - pubmed ID(s) of publications that cite existing variant
