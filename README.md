# Estimating effective reproduction numbers using wastewater data from multiple sewersheds for SARS-CoV-2 in California counties

This repository contains the code for the manuscript ["Estimating effective reproduction numbers using wastewater data from multiple sewersheds for SARS-CoV-2 in California counties"](https://www.sciencedirect.com/science/article/pii/S1755436524000641) (Ravuri et al. 2024). This manuscript details a pipeline to estimate county-aggregated, sewershed-restricted wastewater-based effective reproduction numbers. This repository can inform other jurisdictions’ and public health agencies’ construction of estimation pipelines at distinct geographic scales (i.e., sewershed or county). 

## Code
_ATTN: This repository is designed to run with Rstudio._ 

**Note(s)**: While this code is broadly designed to illustrate the pipeline detailed within the manuscript, it has been modified by the authors to use publicly available wastewater and cases data as much as possible. Given this, analyses relying on sewershed-restricted case data (not publicly available) are commented out; however, all related code is included. Please contact the corresponding author for any queries.

* [01_ww_re_script.Rmd](https://github.com/cdphmodeling/wastewaterReff/blob/main/01_ww_re_script.Rmd): R Markdown file which stores:
    * The pipeline to estimate county-aggregated, sewershed-restricted wastewater effective reproduction numbers. Each step of the pipeline is described, including data import; calculation of cohort and instantaneous effective reproduction numbers at the sewershed level; county-aggregation of sewershed-level effective reproduction numbers; and comparison of wastewater-based estimates against the case-based estimates of [CalCAT](calcat.cdph.ca.gov) (the California Communicable diseases Assessment Tool, a publicly available, county-level ensemble of case-based effective reproduction number estimates).
    * Code to calculate all metrics included within the manuscript (e.g., cross-correlation coefficients, confusion matrix classification accuracy, mean absolute error, spearman rank correlation coefficients, etc.)

* [00_ww_re_script_functions.R](https://github.com/cdphmodeling/wastewaterReff/blob/main/00_ww_re_script_functions.R): A numbered list of all functions necessary to run [01_ww_re_script.Rmd](https://github.com/cdphmodeling/wastewaterReff/blob/main/01_ww_re_script.Rmd)

* [00_gaussian_process_functions.R](https://github.com/cdphmodeling/wastewaterReff/blob/main/00_gaussian_process_functions.R): Functions necessary to execute Step 8 of [01_ww_re_script.Rmd](https://github.com/cdphmodeling/wastewaterReff/blob/main/01_ww_re_script.Rmd) (calculation of county-level effective reproduction number confidence intervals using Gaussian interpolation)

## Data
* The [data](https://github.com/cdphmodeling/wastewaterReff/tree/main/data) folder hosts all all necessary data files to run [01_ww_re_script.Rmd](https://github.com/cdphmodeling/wastewaterReff/blob/main/01_ww_re_script.Rmd). Currently, the only data file in this repository is 'CalCATEnsemble_timeseries.csv' - a historical time series of CalCAT Ensemble data, necessary for comparison of wastewater-based effective reproduction number estimates against a case-based standard.

## Figures & Tables
* The [figures](https://github.com/cdphmodeling/wastewaterReff/tree/main/figures) and [tables](https://github.com/cdphmodeling/wastewaterReff/tree/main/tables) folders host all recreated manuscript figures and tables produced when executing [01_ww_re_script.Rmd](https://github.com/cdphmodeling/wastewaterReff/blob/main/01_ww_re_script.Rmd). Prior to running this script, both folders should be empty. 
