
### AGP_genus_table.csv
- Genus-level (taxonomic level 6) abundance table for the American Gut Project (AGP) cohort.
- Sample IDs were reassigned because the original IDs caused errors during analysis.
- Contains 3,678 samples.

### AGP_HC_index_and_metadata.csv
- HC-index values (calculated using `hc_index_calc_new17taxa.py`) and metadata for samples in `AGP_genus_table.csv`.

### AGP_metadata_10317_20260402-160702.txt
- Metadata for the AGP cohort, obtained from [Qiita Study 10317](https://qiita.ucsd.edu/study/description/10317).
- Direct download link: https://qiita.ucsd.edu/download/3715213

### AGP_sample_ids_candidates_diet.txt
- List of sample IDs generated from `AGP_metadata_10317_20260402-160702.txt` using `id_generation4download_from_AGP.R`.
- Used to download data from redbiom.

### Animals_30JMIC_30Hadza_genus_table.csv
- Genus-level (taxonomic level 6) abundance table for 858 animal samples, 30 individuals from the J-MIC cohort, and 30 individuals from the Hadza cohort.

### Animals_alpha_HC_index.csv
- Alpha diversity and HC-index values for 858 animal samples.
- Alpha diversity was not calculated for samples with fewer than 5,000 reads; missing values are indicated as "na".

### Animals_diet_HC_index_by_species.csv
- Mean HC-index values calculated for each animal species.
- Proportions of plant- and animal-based foods in the diet of each species, calculated using the EltonTraits database (Wilman et al., 2014).
- Species without dietary information in the EltonTraits database are indicated as "na".

### Animals_genus_table.csv
- Genus-level (taxonomic level 6) abundance table for 858 animal samples.

### DF_intervention_HC_index_alpha_feremntation.csv
- HC-index values, alpha diversity, and fermentation parameters from a dietary intervention study (Inoue et al., 2025).

### Hadza_genus_table.csv
- Genus-level (taxonomic level 6) abundance table for the Hadza cohort.
- Raw FASTQ files were downloaded from the European Nucleotide Archive (ENA; accession: PRJEB27517; PMID: 28839072) and processed using our analysis pipeline.

### Hadza_metadata.tsv
- Metadata for samples in `Hadza_genus_table.csv`.

### Human_GreatApe_HC_index.csv
- HC-index values for urbanized individuals from the J-MIC cohort, individuals from the Hadza cohort, and great apes.

### JMIC_genus_table.csv
- Genus-level (taxonomic level 6) abundance table for the J-MIC cohort, filtered to include only samples with calculated alpha diversity values (i.e., samples with more than 5,000 reads).
- Contains 26,783 samples.

### JMIC_HC_index_and_metadata.csv
- HC-index values (calculated using `hc_index_calc_new17taxa.py`) and metadata for samples in `JMIC_genus_table.csv`.

### JMIC_metagenome_CAZY_tpm.tsv
- TPM-normalized abundance table for carbohydrate-active enzyme (CAZyme) genes based on metagenomic analysis of the J-MIC cohort.

### JMIC_metagenome_HC_index.csv
- HC-index values for participants included in the metagenomic analysis of the J-MIC cohort.

### JMIC_metagenome_KO_tpm.tsv
- TPM-normalized abundance table for KEGG Orthology (KO) groups based on metagenomic analysis of the J-MIC cohort.

### Kyotango_genus_table.csv
- Genus-level (taxonomic level 6) abundance table for the Kyotango cohort.
- Contains 770 samples.

### Kyotango_HC_index_and_metadata.csv
- HC-index values (calculated using `hc_index_calc_new17taxa.py`) and metadata for samples in `Kyotango_genus_table.csv`.

### KPUM_genus_table.csv
- Genus-level (taxonomic level 6) abundance table for the KPUM dataset.
- Contains 1803 samples.

### Kyotango_HC_index_and_metadata.csv
- HC-index values (calculated using `hc_index_calc_new17taxa.py`) and metadata for samples in `KPUM_genus_table.csv`.