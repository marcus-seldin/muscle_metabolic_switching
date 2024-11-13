# muscle_metabolic_switching
#### This repository contains all raw data and reproducible scripts used to analyze Tabula Sapiens skeletal muscle gene expression for metabolic gene correlations, as well as a filtered GTEx dataset used to compare covariance of relationships between enzyme expression and pathways
#### Figures produced here contributed to Campos et al. 2024

### Raw data used can be found here: 
#### https://drive.google.com/drive/folders/1Dif_0wEpJrYFzFSzPUpQijxYz5kV-gjr?usp=sharing
#### TS_Muscle.h5ad - contains muscle single-cell sequencing data available here: https://cellxgene.cziscience.com/e/1c9eb291-6d31-47e1-96b2-129b5e1ae64f.cxg/
#### GTEx NA included env.RData - contains a filtered and full GTEx v8 environment to analyze muscle variation

### scripts included 
#### single-cell muscle metabolic expression analyses.R - analyze and plot gene correlations in muslce cell types in single-cell seq
#### mediation scripts for metabolic gene relationships in bulk muscle.R - scripts used to run regression mediation analyses
