# mouse-gut-brain

#### Repository for the analysis of primate GM effects on mouse brain gene expression

This repository contains scripts used in our analysis of primate GM effects on mouse brain gene expression.

We ran most analysis steps using [R](https://cran.r-project.org/) (v4.1). We recommend the following utility or visualization packages to extend base R's functionality.

# Inputs

The following files are expected:

Mouse experiment - brain RNAseq data:
* A .csv file containing count data ```data/counts.csv```
* A .csv file containing metadata data ```data/mouse_meta.csv```
* A .csv file containing technical metadata data ```data/mouse_tech_meta.csv```
* A .csv file containing normalized expression data ```data/exp-data-matched-to-pathways.csv```

Mouse experiment - GM data:
* A .csv file containing pathway data ```data/pathway-data-matched-to-expression.csv```

Published datasets:
* An .rds file containing primate expression data (from Zhu et al. 2018) ```data/human_v_macaque_zhu_MFC.rds```
* A .csv file containing the results from Dumas et al. 2021 ```data/dumas.csv```
* A .tsv file containing disease risk genes from the DISEASES database ```data/human_disease_associations.tsv```
* A .txt file containing network parameters for MiMeNet analysis ```data/network_parameters.txt```
  
# Pipeline

* **Key libraries:** variancePartition, BRETIGEA, biomaRt, ggplot2, corrplot

```
scripts/mouse-code-final.R
```

MiMeNet analysis

* **Key libraries:** MiMeNet_train.py

```
scripts/microbiome.sh
```
