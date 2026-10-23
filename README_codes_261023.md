# Codes

## JMIC_filterling_based_on_AlphaDiversity.R 
  - R code for filtering J-MIC sample based on alpha diversity value (i.e. >5000 reads).


## hc_index_calc_new17taxa.py
  - Python code for calculating HC-index based on final set of 17 taxa.

Usage:

```
python hc_index_calc_new17taxa.py --input input_file --output output_file
```


## HMC_hc_index_calc_new17taxa.py
  - Python code for calculating HC-index based on final set of 17 taxa for Human Microbiome Compendium dataset.

Usage:  

```
python HMC_hc_index_calc_new17taxa.py -i taxonomic_table.csv.gz -o HC_index_compendium.csv
```


## Model_compariosn.R
  - R script for evaluating the additional contribution of microbiome indices to models adjusted for age and sex.
  
  - For categorical outcomes, likelihood ratio tests (LRT) were performed to assess the contribution of microbiome indices beyond covariates (Age and Sex). 
  
  - Full models including the microbiome index were compared to reduced models including only covariates, and the LRT statistic was calculated as:
  
    `LRT = 2 × (log-likelihood of full model − log-likelihood of reduced model)`
  
  - This evaluates whether inclusion of the microbiome index significantly improves model fit.
  
  - While ΔR² quantifies the increase in explained variance for numeric outcomes, LRT assesses the statistical significance of model improvement for categorical outcomes.


## Model_compariosn_for_AGP.R
  - R script for performing the same model comparison analyses as `Model_compariosn.R`, adapted for the American Gut Project (AGP) dataset.
   

## AGP_id_generation_for_download.R
  - R code for sample id generation from metafile.


## AGP_Data_Acquisition_and_Processing.md
  - This markdown file describes the codes used to download biom data from redbiom and to process the downloaded data for QIIME2 (taxonomy assignment by SILVA 138 and alpha diversity calculation).


## AGP_AfterQIIME_Processing.R
  - R code for making abundance table combined with meta data from files processed by above codes. 


## cor_stat_new17taxa_with_DetectionLimit.py
  - Python codes for evaluating correaltion between HC-taxa and alpha diversity indices in 5 cohorts.
  - Performs correlation analyses and generates forest plots to visualize the results.


## HistGramGeneration_ver2.R
  - R code for generating histgram of HC-index from 5 cohorts and data from Human Microbiome Compendium.
  
Usage:  

```
python hc_index_compendium.py -i taxonomic_table.csv.gz -o HC_index_compendium.csv
```


## calc_hc_alpha_region_correlation.py
  - Python code for calculating HC-index and alpha diversity indices, and evaluating correlation between HC-index and alpha diversity indices among region.
  - Only samples from adults were included in the analysis.
  - Projects in which more than 20% of samples had near-zero HC-index values (absolute HC-index < 0.01) were excluded to minimize potential methodological bias.

Usage:  

```
    python calc_hc_alpha_region_correlation.py
    -i taxonomic_table.csv.gz
    -m sample_metadata.tsv
    -p projects.csv
    --sample_col auto
    --region_col region
    --eps 5e-5
    --detection_limit 1e-4
    --min_n 30
```

## plot_hc_by_region.py
  - Python script for examining the distribution of HC-index values across geographic regions, analytical methods, and clinical conditions (e.g., inflammatory bowel disease [IBD]).
  - Only samples from adults were included in the analysis.
  - Projects in which more than 20% of samples had near-zero HC-index values (absolute HC-index < 0.01) were excluded to minimize potential methodological bias.

Usage: 

```
python plot_hc_by_region.py -i HC_index_compendium.csv -m sample_metadata.tsv -p projects.csv -o HC_region_result
```


## analyze_hc_change.R
  - R scripts for calculating baseline-adjusted changes in HC-index in dietary fiber intervention study.



## NOTE:
- The HMC data (`taxonomic_table.csv.gz`, `sample_metadata.tsv`, and `projects.csv`) are available from [Zenodo](https://zenodo.org/records/15122187).
