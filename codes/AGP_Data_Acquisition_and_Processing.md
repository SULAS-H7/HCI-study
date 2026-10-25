# Download AGP data from redbiom
（sample_ids_candidates_diet.txt was generated using "id_generation4download_from_AGP.R")

redbiom fetch samples \
  --from sample_ids_candidates_diet.txt \
  --context Deblur_2021.09-Illumina-16S-V4-90nt-dd6875 \
  --resolve-ambiguities most-reads \
  --output redbiom_diet/Deblur_2021.09-Illumina-16S-V4-90nt-dd6875.biom


redbiom fetch samples \
  --from sample_ids_candidates_diet.txt \
  --context Deblur_2021.09-Illumina-16S-V4-125nt-92f954 \
  --resolve-ambiguities most-reads \
  --output redbiom_diet/Deblur_2021.09-Illumina-16S-V4-125nt-92f954.biom


# Data Processing for QIIME2 analysis
biom convert \
  -i redbiom_diet/Deblur_2021.09-Illumina-16S-V4-90nt-dd6875.biom \
  -o redbiom_diet/90nt.tsv \
  --to-tsv

biom convert \
  -i redbiom_diet/Deblur_2021.09-Illumina-16S-V4-125nt-92f954.biom \
  -o redbiom_diet/125nt.tsv \
  --to-tsv

grep -v '^#' redbiom_diet/90nt.tsv | tail -n +2 | cut -f1 | \
awk '{print ">"$1"\n"$1}' > redbiom_diet/90nt_seqid.fasta

grep -v '^#' redbiom_diet/125nt.tsv | tail -n +2 | cut -f1 | \
awk '{print ">"$1"\n"$1}' > redbiom_diet/125nt_seqid.fasta


# Analyze by QIIME2 (ver 2024.5)
qiime tools import \
  --input-path redbiom_diet/Deblur_2021.09-Illumina-16S-V4-90nt-dd6875.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path diet90_table.qza

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path redbiom_diet/90nt_seqid.fasta \
  --output-path diet90_features.qza

qiime tools import \
  --input-path redbiom_diet/Deblur_2021.09-Illumina-16S-V4-125nt-92f954.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path diet125_table.qza

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path redbiom_diet/125nt_seqid.fasta \
  --output-path diet125_features.qza

qiime feature-classifier classify-sklearn \
  --i-classifier ~/Files4analysis/QIIME2.2024.5/silva-138-99-nb-classifier2.qza \
  --i-reads diet90_features.qza \
  --o-classification diet90_taxonomy.qza

qiime feature-classifier classify-sklearn \
  --i-classifier ~/Files4analysis/QIIME2.2024.5/silva-138-99-nb-classifier2.qza \
  --i-reads diet125_features.qza \
  --o-classification diet125_taxonomy.qza

cat > missing_feature_diet125.txt << 'EOF'
feature-id
TACGTAGGTGGCAAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGTGTAGGCGATATCAAGTCAGAAGTGAAAATTACGGGCTCAACTCGTAACCTGCTTTTGAAACTGACATTCTTGAGTGAA
EOF

qiime feature-table filter-features \
  --i-table diet125_table.qza \
  --m-metadata-file missing_feature_diet125.txt \
  --p-exclude-ids \
  --o-filtered-table diet125_table_nomiss.qza


qiime taxa filter-table \
  --i-table diet125_table_nomiss.qza \
  --i-taxonomy diet125_taxonomy.qza \
  --p-exclude "mitochondria","chloroplast" \
  --o-filtered-table diet125_table_tmp.qza

qiime taxa filter-table \
  --i-table diet125_table_tmp.qza \
  --i-taxonomy diet125_taxonomy.qza \
  --p-include "p__" \
  --o-filtered-table diet125_table_filtered.qza

cat > missing_feature_diet90.txt << 'EOF'
feature-id
TACGTAGGGAGCAAGCGTTGTCCGGAATCATTGGGCGTAAAGCGCGCGTAGGTGGCCACATAAGTCCGTTGTGAAAGTCAAAGGCTCAAC
EOF



qiime feature-table filter-features \
  --i-table diet90_table.qza \
  --m-metadata-file missing_feature_diet90.txt \
  --p-exclude-ids \
  --o-filtered-table diet90_table_nomiss.qza


qiime taxa filter-table \
  --i-table diet90_table_nomiss.qza \
  --i-taxonomy diet90_taxonomy.qza \
  --p-exclude "mitochondria","chloroplast" \
  --o-filtered-table diet90_table_tmp.qza

qiime taxa filter-table \
  --i-table diet90_table_tmp.qza \
  --i-taxonomy diet90_taxonomy.qza \
  --p-include "p__" \
  --o-filtered-table diet90_table_filtered.qza


qiime taxa collapse \
  --i-table diet125_table_filtered_singleton.qza \
  --i-taxonomy diet125_taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table diet125_L6.qza

qiime taxa collapse \
  --i-table diet90_table_filtered_singleton.qza \
  --i-taxonomy diet90_taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table diet90_L6.qza


qiime feature-table relative-frequency \
  --i-table diet125_L6.qza \
  --o-relative-frequency-table diet125_L6_rel.qza

qiime feature-table relative-frequency \
  --i-table diet90_L6.qza \
  --o-relative-frequency-table diet90_L6_rel.qza

qiime tools export \
  --input-path diet125_L6_rel.qza \
  --output-path diet125_export

qiime tools export \
  --input-path diet90_L6_rel.qza \
  --output-path diet90_export

biom convert \
  -i diet125_export/feature-table.biom \
  -o diet125_L6_rel.tsv \
  --to-tsv

biom convert \
  -i diet90_export/feature-table.biom \
  -o diet90_L6_rel.tsv \
  --to-tsv

qiime diversity alpha \
  --i-table diet125_table_filtered_singleton.qza \
  --p-metric shannon \
  --o-alpha-diversity shannon_vector.qza

qiime diversity alpha \
  --i-table diet125_table_filtered_singleton.qza \
  --p-metric chao1 \
  --o-alpha-diversity chao1_vector.qza

qiime tools export \
  --input-path shannon_vector.qza \
  --output-path shannon_export125

qiime tools export \
  --input-path chao1_vector.qza \
  --output-path chao1_export125


qiime diversity alpha \
  --i-table diet90_table_filtered_singleton.qza \
  --p-metric shannon \
  --o-alpha-diversity shannon_vector.qza

qiime diversity alpha \
  --i-table diet90_table_filtered_singleton.qza \
  --p-metric chao1 \
  --o-alpha-diversity chao1_vector.qza

qiime tools export \
  --input-path shannon_vector.qza \
  --output-path shannon_export90

qiime tools export \
  --input-path chao1_vector.qza \
  --output-path chao1_export90



# Further processing and integration with metadata by R codes "AfterQIIMEProcessing2tsvExport_ver2.R"