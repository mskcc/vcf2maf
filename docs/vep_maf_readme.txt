The first 34 columns are standard MAF format, and described here:
https://wiki.nci.nih.gov/x/eJaPAQ

The subsequent 10 columns are relevant to most analyses.

35. HGVSc - the coding sequence of the variant in HGVS recommended format
36. HGVSp - the protein sequence of the variant in HGVS recommended format
37. HGVSp_Short - Same as HGVSp, but using 1-letter amino-acid codes
38. Transcript_ID - transcript onto which the consequence of the variant has been mapped
39. Exon_Number - the exon number (out of total number)
40. t_depth - read depth across this locus in tumor BAM
41. t_ref_count - read depth supporting the reference allele in tumor BAM
42. t_alt_count - read depth supporting the variant allele in tumor BAM
43. n_depth - read depth across this locus in normal BAM
44. n_ref_count - read depth supporting the reference allele in normal BAM
45. n_alt_count - read depth supporting the variant allele in normal BAM

The next column is relevant to analyses that consider the effect of the variant on all alternate
isoforms of the gene, or on non-coding/regulatory transcripts. The effects are sorted first by
transcript biotype priority, then by effect severity, and finally by decreasing order of transcript
length. Each effect in the list is in the format [SYMBOL,Consequence,HGVSp,Transcript_ID,RefSeq].

46. all_effects - a semicolon delimited list of all possible variant effects, sorted by priority

All remaining columns are straight out of Ensembl's VEP annotator, as described here:
http://useast.ensembl.org/info/docs/tools/vep/vep_formats.html#output

47. Allele - the variant allele used to calculate the consequence
48. Gene - stable Ensembl ID of affected gene
49. Feature - stable Ensembl ID of feature
50. Feature_type - type of feature. Currently one of Transcript, RegulatoryFeature, MotifFeature
51. Consequence - consequence type of this variation
52. cDNA_position - relative position of base pair in cDNA sequence
53. CDS_position - relative position of base pair in coding sequence
54. Protein_position - relative position of amino acid in protein
55. Amino_acids - only given if the variation affects the protein-coding sequence
56. Codons - the alternative codons with the variant base in upper case
57. Existing_variation - known identifier of existing variation
58. ALLELE_NUM - allele number from input; 0 is reference, 1 is first alternate etc
59. DISTANCE - shortest distance from variant to transcript
60. STRAND - the DNA strand (1 or -1) on which the transcript/feature lies
61. SYMBOL - the gene symbol
62. SYMBOL_SOURCE - the source of the gene symbol
63. HGNC_ID - gene identifier from the HUGO Gene Nomenclature Committee
64. BIOTYPE - biotype of transcript
65. CANONICAL - a flag indicating if the transcript is denoted as the canonical transcript for this gene
66. CCDS - the CCDS identifer for this transcript, where applicable
67. ENSP - the Ensembl protein identifier of the affected transcript
68. SWISSPROT - UniProtKB/Swiss-Prot identifier of protein product
69. TREMBL - UniProtKB/TrEMBL identifier of protein product
70. UNIPARC - UniParc identifier of protein product
71. RefSeq - RefSeq identifier for this transcript
72. SIFT - the SIFT prediction and/or score, with both given as prediction (score)
73. PolyPhen - the PolyPhen prediction and/or score
74. EXON - the exon number (out of total number)
75. INTRON - the intron number (out of total number)
76. DOMAINS - the source and identifer of any overlapping protein domains
77. GMAF - Non-reference allele and frequency of existing variant in 1000 Genomes
78. AFR_MAF - Non-reference allele and frequency of existing variant in 1000 Genomes combined African population
79. AMR_MAF - Non-reference allele and frequency of existing variant in 1000 Genomes combined American population
80. ASN_MAF - Non-reference allele and frequency of existing variant in 1000 Genomes combined Asian population
81. EAS_MAF - Non-reference allele and frequency of existing variant in 1000 Genomes combined East Asian population
82. EUR_MAF - Non-reference allele and frequency of existing variant in 1000 Genomes combined European population
83. SAS_MAF - Non-reference allele and frequency of existing variant in 1000 Genomes combined South Asian population
84. AA_MAF - Non-reference allele and frequency of existing variant in NHLBI-ESP African American population
85. EA_MAF - Non-reference allele and frequency of existing variant in NHLBI-ESP European American population
86. CLIN_SIG - clinical significance of variant from dbSNP
87. SOMATIC - somatic status of existing variation(s)
88. PUBMED - pubmed ID(s) of publications that cite existing variant
89. MOTIF_NAME - the source and identifier of a transcription factor binding profile aligned at this position
90. MOTIF_POS - the relative position of the variation in the aligned TFBP
91. HIGH_INF_POS - a flag indicating if the variant falls in a high information position of a transcription factor binding profile (TFBP)
92. MOTIF_SCORE_CHANGE - the difference in motif score of the reference and variant sequences for the TFBP
93. IMPACT - the impact modifier for the consequence type
94. PICK - indicates if this block of consequence data was picked by VEP's pick feature
95. VARIANT_CLASS - Sequence Ontology variant class
96. TSL - Transcript support level
97. HGVS_OFFSET - Indicates by how many bases the HGVS notations for this variant have been shifted
98. PHENO - Indicates if existing variant is associated with a phenotype, disease or trait
99. MINIMISED - Alleles in this variant have been converted to minimal representation before consequence calculation
100. ExAC_AF - Global Allele Frequency from ExAC
101. ExAC_AF_AFR - African/African American Allele Frequency from ExAC
102. ExAC_AF_AMR - American Allele Frequency from ExAC
103. ExAC_AF_EAS - East Asian Allele Frequency from ExAC
104. ExAC_AF_FIN - Finnish Allele Frequency from ExAC
105. ExAC_AF_NFE - Non-Finnish European Allele Frequency from ExAC
106. ExAC_AF_OTH - Other Allele Frequency from ExAC
107. ExAC_AF_SAS - South Asian Allele Frequency from ExAC
108. GENE_PHENO - Indicates if gene that the variant maps to is associated with a phenotype, disease or trait
109. FILTER - False-positive filtering status, first borrowed from the input MAF/VCF, and then sprinkled with magic
