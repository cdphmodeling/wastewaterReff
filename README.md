# Estimating effective reproduction numbers using wastewater data from multiple sewersheds for SARS-CoV-2 in California counties

This repository contains the code for the manuscript "Estimating effective reproduction numbers using wastewater data from multiple sewersheds for SARS-CoV-2 in California counties" by Ravuri et al. (2024). This manuscript details a pipeline to estimate county-aggregated, sewershed-restricted wastewater effective reproduction numbers. This repository can inform other jurisdictions’ and public health agencies’ construction of wastewater-based effective reproduction number estimation pipelines at the sewershed or county level.

## Code & Data
_ATTN: This repository is designed to run with Rstudio._ 

* The  pipeline to estimate county-aggregated, sewershed-restricted wastewater effective reproduction numbers is stored within the file labeled 'ww_re_script.Rmd' (R Markdown). All nine steps of the pipeline are described, including data import, calculation of cohort and instantaneous effective reproduction numbers at the sewershed level, county-aggregation of sewershed-level effective reproduction numbers, and comparison of wastewater-based estimates against [CalCAT](calcat.cdph.ca.gov) (the California Communicable diseases Assessment Tool, a publicly available, county-level ensemble of case-based effective reproduction number estimates). 

* General functions necessary to run 'ww_re_script.Rmd' are located in 'ww_re_script_functions.R' as a numbered list with function descriptions. 

* Functions necessary to execute Step 8 of the estimation pipeline (calculation of county-level effective reproduction number confidence intervals using a Gaussian process) are included in the file 'gaussian_process_functions.R'

* The only data file included in this repository is 'CalCATEnsemble_timeseries.csv'; this file is a historical time series of CalCAT Ensemble data, necessary for comparison of wastewater-based effective reproduction number estimates against a case-based standard. 

_Note: While this code is broadly designed to illustrate the pipeline detailed within the manuscript, it has been modified by the authors to include publicly available wastewater and cases data as much as possible. Given this, analyses relying on sewershed-restricted case data (not publicly available) are excluded. Please contact the corresponding author for any queries._

## Updates
_12/17/2024:_ The repository will be updated in the future to include code for producing manuscript figures and tables. 