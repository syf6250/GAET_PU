# GAET_PU

This repository contains the R code and data used in the paper on GAET-based methods for positive-unlabeled (PU) data analysis.

The repository includes simulation studies, real-data analyses, and semi-synthetic data experiments.

---

## Repository Structure

---

## simulation/

R scripts for simulation experiments corresponding to the tables and figures in the paper.

- `Tab1_Set1.R` – Results for Setting 1 in Table 1  
- `Tab1_Set2.R` – Results for Setting 2 in Table 1  
- `Tab1_Set3.R` – Results for Setting 3 in Table 1  
- `Tab1_lin.R` – Linear method results for Table 1  
- `Tab1_Bayes.R` – Bayes classifier results for Table 1  
- `FigS1S2_Set1.R` – Data for Figures S1 and S2 (Setting 1)  
- `FigS1S2_Set2.R` – Data for Figures S1 and S2 (Setting 2)  
- `FigS1S2_Set3.R` – Data for Figures S1 and S2 (Setting 3)

---

## results/

Scripts generating final tables and figures reported in the paper.

- `Fig1-FigS3.R` – Figures 1 and S3  
- `FigS1S2.R` – Figures S1 and S2  
- `Fig2.R` – Figure 2  
- `Tab2_Set1.R` – Setting 1 in Table 2  
- `Tab2_Set2.R` – Setting 2 in Table 2  
- `Tab2_Set3.R` – Setting 3 in Table 2  
- `Tab3_Set1.R` – Setting 1* in Table 3  
- `Tab3_Set2.R` – Setting 2* in Table 3  
- `TabS1.R` – Table S1  
- `TabS2-FigS4.R` – Table S2 and Figure S4  

Other auxiliary files in this folder contain necessary supporting functions.

---

## real_data/

R scripts for real-data applications.

- `survey_analysis.R` – Survey data analysis and Figure S5  
- `presence-only_analysis.R` – Presence-only data analysis and Figure S6  

Survey data are publicly available at:  
https://search.gesis.org/research_data/ZA5688

---

## semi_synthetic_data/

R scripts and datasets for semi-synthetic experiments.

Datasets:
- `wilt_training.csv`, `wilt_testing.csv`
- `spambase.data`

Scripts and outputs:
- `TabS4_wilt.R` – Table S4 (wilt data)  
- `TabS4_spam.R` – Table S4 (spambase data)  
- `TabS5.R` – Table S5  
- `TabS6.R` – Table S6  

Some files in this directory contain auxiliary functions required for computation.

---

## Reproducibility

All simulation and empirical results in the paper can be reproduced by running the corresponding R scripts in each folder.

Users are encouraged to run scripts sequentially according to the tables and figures reported in the manuscript.

---

## License

This repository is released under the MIT License.
