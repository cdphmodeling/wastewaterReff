# General Functions (Numbered List)

# 1.1. Import publicly available wastewater data

read_in_ww_data <- function(county_list) {
  
  wastewater_url <- "https://data.chhs.ca.gov/dataset/a6ca879a-6014-4b72-9ea6-07ef8b87ae83/resource/2742b824-3736-4292-90a9-7fad98e94c06/download/wastewatersurveillancecalifornia.csv"
  
  ww_df_all <-
    read_csv(wastewater_url)%>% 
    mutate(date = lubridate::mdy(sample_collect_date))
  
  # Reshape & filter wastewater data
  
  ww_df <- ww_df_all %>% dplyr::select(data_source,
                                       label_name,
                                       wwtp_name,
                                       sample_date = sample_collect_date,
                                       population_served,
                                       pcr_target, 
                                       pcr_gene_target, 
                                       raw_concentration = pcr_target_avg_conc, 
                                       County_address = county_treatmentplant) %>% 
    mutate(date = lubridate::parse_date_time(sample_date, c("%d/%m/%Y","%m/%d/%Y", "%Y-%m-%d"))) %>%
    filter(pcr_gene_target %in% c("n"),
           data_source %in% c("Sewage Coronavirus Alert Network (SCAN)", "WastewaterSCAN", "Healthy Central Valley Together (HCVT)"),
           wwtp_name %notin% "CODIGA_All", #Codiga Resource Recovery Center on the Stanford campus is excluded from analysis
           County_address %in% county_list) %>% 
    mutate(data_source_short = ifelse(data_source %in% c("Sewage Coronavirus Alert Network (SCAN)", "WastewaterSCAN"), "scan", "hcvt")) %>%
    arrange(group_by=sample_date)

  return(ww_df)
}

# 1.2. Import publicly available sewershed-restricted case data 

read_in_case_data <- function(sewershed_list) {
  
  # Sewershed-restricted case data does not map sewersheds to corresponding counties served; creating a sewershed to county crosskey is a necessary first step
  
  wastewater_url <- "https://data.chhs.ca.gov/dataset/a6ca879a-6014-4b72-9ea6-07ef8b87ae83/resource/2742b824-3736-4292-90a9-7fad98e94c06/download/wastewatersurveillancecalifornia.csv"
  ww_df_all <-
    read_csv(wastewater_url)%>% 
    mutate(date = lubridate::mdy(sample_collect_date))
  sewershed_county_key <- ww_df_all %>% group_by(wwtp_name) %>% summarise(County_address = unique(county_treatmentplant), label_name = unique(label_name))
  
  # Import sewershed-restricted case data
  case_url <- "https://data.chhs.ca.gov/dataset/a6aadb84-9c1a-4ddc-9d09-6df9e476768a/resource/e148f208-3069-4dc6-a98b-ef6d29dfbc16/download/outputfile.csv"
  
  # Reshape & filter cases data
  case_df_all <-
    read_csv(case_url) 
  
  case_df <- case_df_all %>% 
    rename(wwtp_name = Area, sample_date = date, raw_cases = Cases, label_name = Label_Name) %>%
    filter(wwtp_name %in% sewershed_list) %>%
    full_join(sewershed_county_key, by = c("wwtp_name", "label_name")) %>% 
    dplyr::ungroup() %>%
    dplyr::group_by(sample_date) %>%
    distinct(County_address, wwtp_name, sample_date, raw_cases, .keep_all = T) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(wwtp_name) %>% 
    dplyr::select(sample_date, County_address, wwtp_name, raw_cases) 
  
  return(case_df)
}

# 2. Apply smoothing spline, single user-defined knot value

# Sample input values to test function 
# df <- ww_df %>% filter(wwtp_name == "Turlock_RWQCF") %>% filter(pcr_gene_target == "n")
# y_var <- "raw_concentration"
# knot_spacing <- 7

ww_smoothing <- function(df=ww_df, knot_spacing, y_var = "raw_concentration") {
  
  # Set knots to be located every n days in the dataset (n is defined by the function input variable, knot_spacing). Last knot should be located a few days before the last date of the analysis time period
  knot_locations <- seq(from = min(df$sample_date), to = max(df$sample_date)-3, by = knot_spacing) 
  
  # Add a variable that fills on all non-sampling knot dates with the most recent data point
  df$concentration_var <- df[[y_var]]
  
  df <- df %>%
    mutate(concentration_fill = concentration_var) %>%
    fill(concentration_fill, .direction = 'down') %>%
    mutate(knot_date = ifelse(sample_date %in% knot_locations, TRUE, FALSE)) %>% # Indicator variables for sample dates with knots
    mutate(raw_spline_input = ifelse(is.na(concentration_var) & knot_date, concentration_fill, concentration_var)) # If no sample on a knot date, set to most recent sample
  
  # Filter any NAs in raw_spline_input
  df_for_model <- df %>% filter(!is.na(raw_spline_input))
  
  # Make a date vector for the downstream predict() function
  all_dates <- seq(from = min(df$sample_date), to = max(df$sample_date), by = "1 day")
  
  # Fit smoothing spline to dataset with no NA's, using OCV (ordinary cross validation) to set smoothing parameters
  spline_ocv <- ss(x= df_for_model$sample_date, y = df_for_model$raw_spline_input, method = "OCV", knots = knot_locations)
  
  # Use predict function to get fitted spline values for all dates in the original dataset
  fitted_ocv <- predict(spline_ocv, x = all_dates, se.fit = F)
  
  # Create new output data frame ("out") with predicted spline values
  out <- data.frame(sample_date = fitted_ocv$x,
                    fit_ocv = fitted_ocv$y,
                    data_source_short = unique(df$data_source_short),
                    County_address = unique(df$County_address),
                    wwtp_name = unique(df$wwtp_name),
                    pcr_gene_target = unique(df$pcr_gene_target))
  
  # Identify the column with fitted value
  fit_cols <- colnames(out)[grepl("fit_",colnames(out))]
  
  # If the spline predicted a negative number, set the value to zero for corresponding sample date
  out <- out %>% mutate(across(all_of(fit_cols), function(fit_col) ifelse(fit_col<0, 0, fit_col)))
  
  # Rename column
  colnames(out)[grepl("fit_", colnames(out))] <- paste0(colnames(out)[grepl("fit_", colnames(out))], '_', knot_spacing, "day_knots")

  return(out)
}

# 3. Apply smoothing spline, multiple user-defined knot values

no_knots_vector <- c(7, 10)

ww_data_with_smoothed_concentrations <- function(df=ww_df, no_knots_vector, normalization_var) {
  
  for(i in 1:length(no_knots_vector)) {
    # Create data frame with spline-smoothed wastewater concentrations for each knot value
    no <- no_knots_vector[i]
    out_temp <- ww_smoothing(df, no, y_var = normalization_var)
    # Merge output data frame for each knot value into a single master data frame for multiple knot values
    if(i == 1){out_df = out_temp}else{out_df <- out_df %>%
      left_join(out_temp, by = c("data_source_short", "wwtp_name", "pcr_gene_target", "sample_date"))}
  }
  
  return(out_df)
}

# 4. Determine gamma parameters based on the mean/sd values of a specific distribution
# Based on following source: https://ehp.niehs.nih.gov/doi/10.1289/EHP10050

getGammaParams <- function(meanParam, sdParam){
  shapeParam <- meanParam^2 / (sdParam^2)
  scaleParam <- (sdParam^2) / meanParam
  return(list(shape = shapeParam, scale = scaleParam))
}

# 5. Define infection to shedding delay distrubtion (outputs a list object)

define_infection_to_shedding <- function(mean, sd) {
  #' SARS-CoV-2 Delay between infection to shedding into wastewater
  sars_cov_2_distribution_infection_to_shedding <- list(
    name = "gamma",
    shape = getGammaParams(mean, sd)$shape,
    scale = getGammaParams(mean, sd)$scale)
  
  return(sars_cov_2_distribution_infection_to_shedding)
  
}

# 6. Determine log normal parameters based on the mean/sd values of a specific distribution
getLogNormalParams <- function(meanParam, sdParam){
  logmeanParam <- log(meanParam^2 / sqrt(sdParam^2 + meanParam^2))
  logSDParam <- sqrt(log(1 + (sdParam^2 / meanParam^2)))
  return(list(logMean = logmeanParam, logSD = logSDParam))
}

# 7. Define incubation period and symptom onset to reporting delay distrubtions (outputs a list object) for estimation of case-based effective reproduction numbers

define_incubation <- function(mean, sd) {
  #' SARS-CoV-2  Delay between infection and onset of symptoms (incubation period) in days
  sars_cov_2_distribution_incubation <- list(
    name = "lnorm",
    meanlog = getLogNormalParams(mean, sd)$logMean,
    sdlog = getLogNormalParams(mean,sd)$logSD)
  
  return(sars_cov_2_distribution_incubation)
  
}

define_onset_to_reporting <- function(mean, sd) {
  sars_cov_2_distribution_onset_to_reporting <- list(
    name = "gamma",
    shape =  getGammaParams(mean, sd)$shape,
    scale =  getGammaParams(mean, sd)$scale)
  
  return(sars_cov_2_distribution_onset_to_reporting)
  
}

# 8. Estimation of effective reproduction numbers using a single input variable

# Sample input values to test function 
# df <- ww_data_analysis_period %>% filter(wwtp_name == "PA_all")
# input_col <- "fit_ocv_7day_knots"

reff_single_input_col <- function(df, input_col) {
  if(input_col == "raw_cases"){
  # Input column may be cases or wastewater. If the user posesses cases data of interest, incubation period and symptom onset to reporting delay distributions must be defined. LOESS smoothing is subsequently applied to the time series of case data. 
      delay = list(sars_cov_2_distribution_incubation, sars_cov_2_distribution_onset_to_reporting) # Refer to Functions 6 and 7; necessary for deconvolution step in estimateR
      deconvolution_input_df = df %>% ungroup() %>% dplyr::select(sample_date, all_of(input_col)) %>% arrange(sample_date) %>% na.omit() 
      deconvolution_input_pre_smoothing = deconvolution_input_df %>% pull(input_col)
      
      # Perform LOESS smoothing on deconvoluted time series (unlike estimateR (instantaneous), R0 (cohort) does not execute in-built LOESS smoothing during effective reproduction number estimation)
      deconvolution_input <- smooth_incidence(
        incidence_data = deconvolution_input_pre_smoothing,
        smoothing_method = "LOESS")
      minimum_cumul_incidence=12 # Necessary input variable for 'get_block_bootstrapped_estimate' function below

  } else{
    # The below applies to wastewater input data. 
    delay = list(sars_cov_2_distribution_infection_to_shedding)
    deconvolution_input_df = df %>% ungroup() %>% dplyr::select(sample_date, all_of(input_col)) %>% arrange(sample_date) 
    deconvolution_input = deconvolution_input_df %>% pull(input_col)
    minimum_cumul_incidence=12
  }
  
  # Run estimateR to calculate instantaneous effective reproduction numbers (in-built deconvolution)
  # Refer to https://covid-19-re.github.io/estimateR/ for further details on notation
  estimateR_input = list(
    values = deconvolution_input,
    index_offset = 0)
  
  estimateR_output <- get_block_bootstrapped_estimate(
    estimateR_input$values,
    N_bootstrap_replicates = 1000,
    smoothing_method = "none",
    deconvolution_method = "Richardson-Lucy delay distribution",
    estimation_method = "EpiEstim sliding window",
    estimation_window = 3, #EpiEstim default
    uncertainty_summary_method = "original estimate - CI from bootstrap estimates",
    combine_bootstrap_and_estimation_uncertainties = TRUE,
    delay = delay,
    minimum_cumul_incidence = minimum_cumul_incidence,
    mean_serial_interval = sars_cov_2_mean_serial_interval_days,
    std_serial_interval = sars_cov_2_std_serial_interval_days,
    ref_date = 
      as.Date(first((deconvolution_input_df$sample_date))),
    time_step = "day",
    output_Re_only = F) %>%
    mutate(Re_model = "EpiEstim", input = input_col, .after = date) %>%
    rename(sample_date = date) %>%
    dplyr::select(sample_date, 
                  Re_model, 
                  input, 
                  deconvolved_incidence, 
                  Re_estimate, CI_down_Re_estimate, 
                  CI_up_Re_estimate, CI_down_deconvolved_incidence, CI_up_deconvolved_incidence) 
  
  # Run R0 to calculate cohort effective reproduction numbers (no in-built deconvolution)
  # Extract deconvoluted input time series from estimateR in-built deconvolution
  WallingaTeunis_input <- estimateR_output %>% dplyr::select(sample_date, deconvolved_incidence, CI_down_deconvolved_incidence, CI_up_deconvolved_incidence) %>% filter(!is.na(deconvolved_incidence))
  deconvolution_output_as_named_vector <- setNames(WallingaTeunis_input$deconvolved_incidence, WallingaTeunis_input$sample_date)
  
  # Calculate effective reproduction number 
  WallingaTeunis_output <- R0::est.R0.TD(deconvolution_output_as_named_vector, mGT2, begin = as.integer(1), end = length(deconvolution_output_as_named_vector), nsim = 1000, correct = TRUE)
  
  # Replace last estimated reproduction number with forecasted value + this value's confidence intervals with prediction intervals (smooths abrupt tailing estimate outlier, artifact produced by Wallinga&Teunis)
  fix_tailing_Re_estimates = forecast(WallingaTeunis_output$R[-c((length(WallingaTeunis_output$R) - 1):length(WallingaTeunis_output$R))])
  WallingaTeunis_output$R[c((length(WallingaTeunis_output$R) - 1):length(WallingaTeunis_output$R))] = fix_tailing_Re_estimates$mean[1:2]
  
  WallingaTeunis_output$conf.int[,1][c((length(WallingaTeunis_output$conf.int[,1]) - 1):length(WallingaTeunis_output$conf.int[,1]))] = fix_tailing_Re_estimates$lower[,2][1:2]
  WallingaTeunis_output$conf.int[,2][c((length(WallingaTeunis_output$conf.int[,2]) - 1):length(WallingaTeunis_output$conf.int[,2]))] = fix_tailing_Re_estimates$upper[,2][1:2]
  
  # Create data frame from R0 output
  WallingaTeunis_output <- data.frame(sample_date = as.Date(names(WallingaTeunis_output$R)),
                                      Re_model = "Wallinga&Teunis",
                                      input = input_col,
                                      Re_estimate = WallingaTeunis_output$R,
                                      CI_down_Re_estimate = WallingaTeunis_output$conf.int[,1],
                                      CI_up_Re_estimate = WallingaTeunis_output$conf.int[,2]
  )
  
  WallingaTeunis_output <- left_join(WallingaTeunis_output, WallingaTeunis_input)
  
  # Join estimateR & R0 outputs into a single data frame
  reff <- full_join(WallingaTeunis_output, estimateR_output) %>% mutate(wwtp_name = unique(df$wwtp_name), County_address = unique(df$County_address), .after = sample_date)
  
  return(reff)
  
}


# 9. Estimate effective reproduction numbers for multiple input variables of interest 
reff_all_input_cols <- function(df, input_cols_as_vector) {
  # Estime effective reproduction numbers for each input variable of interest independently
  input_cols_as_dataframe <- data.frame(inputs = input_cols_as_vector)
  for(i in input_cols_as_dataframe$inputs) {
    if(i == input_cols_as_dataframe$inputs[1]) {
      out_df <- reff_single_input_col(df, i)
    } else{
      # Merge data frames for each input variable of interest into a single output data frame
      out_temp <- reff_single_input_col(df, i)
      out_df <- full_join(out_df, out_temp)
    }
  }

  return(out_df)
}

# 10. Functions to calculate county-level populations (relevant for Step 8 in calculation pipeline)
# Sourced from: https://github.com/StateOfCalifornia/CalCAT/tree/master/R

get_counties <- function(State = state_name){

counties <- data.table::fread("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv") %>% 
  .[state == State, .(fips = unique(fips)), by = .(county, state)] %>% 
  na.omit() %>%
  .[, .(county, state, fips, fips_char = ifelse(fips < 9999, paste0("0",fips),as.character(fips)))] %>%
  merge(data.table::fread("https://raw.githubusercontent.com/CDPHrusers/List-of-US-States/master/states.csv"), by.x = "state", by.y = "State")

return(counties)

}

get_state_fips <- function(State = state_name, type = "integer"){
  states <- data.table(state = c("Alabama",
                                 "Alaska",
                                 "Arizona",
                                 "Arkansas",
                                 "California",
                                 "Colorado",
                                 "Connecticut",
                                 "Delaware",
                                 "District of Columbia",
                                 "Florida",
                                 "Georgia",
                                 "Hawaii",
                                 "Idaho",
                                 "Illinois",
                                 "Indiana",
                                 "Iowa",
                                 "Kansas",
                                 "Kentucky",
                                 "Louisiana",
                                 "Maine",
                                 "Maryland",
                                 "Massachusetts",
                                 "Michigan",
                                 "Minnesota",
                                 "Mississippi",
                                 "Missouri",
                                 "Montana",
                                 "Nebraska",
                                 "Nevada",
                                 "New Hampshire", 
                                 "New Jersey",
                                 "New Mexico",
                                 "New York",
                                 "North Carolina",
                                 "North Dakota",
                                 "Ohio",
                                 "Oklahoma", 
                                 "Oregon",
                                 "Pennsylvania",
                                 "Rhode Island",
                                 "South Carolina",
                                 "South Dakota",
                                 "Tennessee",
                                 "Texas",
                                 "Utah",
                                 "Vermont",
                                 "Virginia",
                                 "Washington", 
                                 "West Virginia",
                                 "Wisconsin",
                                 "Wyoming"), 
                       fips = c(
                         1,
                         2,
                         4,
                         5,
                         6,
                         8,
                         9,
                         10,
                         11,
                         12,
                         13,
                         15,
                         16,
                         17,
                         18,
                         19,
                         20,
                         21,
                         22,
                         23,
                         24,
                         25,
                         26,
                         27,
                         28,
                         29,
                         30,
                         31,
                         32,
                         33,
                         34,
                         35,
                         36,
                         37,
                         38,
                         39,
                         40,
                         41,
                         42,
                         44,
                         45,
                         46,
                         47,
                         48,
                         49,
                         50,
                         51,
                         53,
                         54,
                         55,
                         56
                       ))
  
  
  num <-  states$fips[which(states$state == State)]
  char <- ifelse(num<10,as.character(paste0("0",num)),as.character(num))
  
  return({
    if(type == "integer"){
      as.integer(num)
    } else if(type == "character"){
      as.character(char)
    } else {errorCondition("Please select type = 'integer' or 'character'")}
    
  })
  
}

get_county_populations <- function(State = state_name){
  
  goo <- get_counties(State) %>% .[, .(county, fips = as.integer(str_sub(fips, -3)))]
  
  foo <- fread("https://www2.census.gov/programs-surveys/popest/datasets/2010-2019/counties/totals/co-est2019-alldata.csv") %>%
    .[STNAME == State, .(fips = ifelse(COUNTY<10,paste0(get_state_fips(State = State, type = "character"),"00",COUNTY),
                                       ifelse(COUNTY<100,paste0(get_state_fips(State = State, type = "character"),"0",COUNTY),
                                       ifelse(COUNTY>100,paste0(get_state_fips(State = State, type = "character"), COUNTY), COUNTY))), 
                         pop2020 = as.double(POPESTIMATE2019))] %>% 
    .[fips == 0, county := State] %>% select(fips, pop2020) %>% tibble()
  
  return(foo)  
  
}

CA_county_populations <- function(counties_of_interest){
counties <- get_counties("California") %>% 
  filter(county %in% counties_of_interest) %>% 
  dplyr::select(-c(fips)) %>% 
  rename(fips = fips_char)
county_pop = right_join(get_county_populations("California"), counties, by = c("fips"))}

# 11. Calculate cross-correlation coefficients between case-based and wastewater-based effective reproduction number estimates for a single location

# Sample input values to test function 
# df = all_re_df %>% filter(location=="Santa Clara")
# models = c("EpiEstim", "Wallinga&Teunis")
# calcat_comparison = T
# raw_cases_comparison = F
# inputs_vec = c("root_transform")

ccf_analysis <- function(df, 
                         models = c("EpiEstim", "Wallinga&Teunis"), 
                         calcat_comparison = T, # Does the user want to compute cross-correlation coefficients between wastewater-based effective reproduction numbers and the CalCAT Ensemble?
                         raw_cases_comparison = F, # Does the user want to compute cross-correlation coefficients between wastewater- and sewershed-restricted, case-based effective reproduction numbers?
                         inputs_vec # User-defined vector of either: exclusively the wastewater-based effective reproduction number time series, or two time series variables for which to compute cross-correlation coefficients (e.g., wastewater- and case-based effective reproduction numbers). In the case of inputting exclusively the wastewater-based effective reproduction number time series, the user is exclusively interested in a comparison against the CalCAT ensemble.
                         ) {
  # Initialize output data frame
  lag_analysis = data.frame(location=character(),
                            location_level=character(),
                            stationary_variable = character(),
                            sliding_variable = character(),
                            Re_model=character(),
                            lag=numeric(), 
                            acf=numeric()) 
  
  # Use the below 'for' loop to calculate metrics comparing two time series of interest (in the manuscript, this loop was utilized to compare wastewater-based and sewershed-restricted, case-based effective reproduction number estimates)
  for(i in models) {
    if(raw_cases_comparison == T){
    for(j in 1:ncol(combn(inputs_vec,2))){
      inputs = combn(inputs_vec,2)[,j] %>% sort()
      input_1 = inputs[1]
      input_2 = inputs[2]
      
      INPUT1 <- df %>% filter(Re_model == i, input == input_1) %>% pull(Re_estimate)
      INPUT2 <- df %>% filter(Re_model == i, input == input_2) %>% pull(Re_estimate)
      if( length(INPUT1)==0 | length(INPUT2)==0){next}
      
      # Rank effective reproduction number estimates prior to computing Pearson cross-correlation coefficients 
      INPUT1 <- frank(INPUT1)
      INPUT2 <- frank(INPUT2)
      
      # Compute cross-correlation coefficients 
      corr<-ccf(INPUT2, INPUT1, lag.max=20, 
                main=paste0(unique(df$location),
                            ": ", i))
      # Store results of cross-correlation in temporary data frame
      lag_analysis_temp <- data.frame(location=unique(df$location),
                                      location_level=unique(df$location_level),
                                      stationary_variable = input_1,
                                      sliding_variable = input_2,
                                      Re_model= i,
                                      lag=corr$lag, 
                                      acf=corr$acf) 
      # Merge temporary data frame with initialized output data frame
      lag_analysis = lag_analysis %>% bind_rows(lag_analysis_temp)
    }}
    
    # Use the below 'if' statement to compute cross-correlation coefficients between wastewater-based effective reproduction numbers and the CalCAT ensemble
    # The below 'if' statement will also calculate cross-correlation coefficients between the sewershed-restricted case- and CalCAT-based effective reproduction numbers as a secondary product
    if(calcat_comparison == T){  # Only T if user is conducting county-level analyses; for sewershed-level, this should be F
      for(k in 1:length(inputs_vec)){
        input_1 = inputs_vec[k]
        input_2 = "CalCAT"
        # Remaining steps are identical to prior 'for' loop for two user-defined input time series
        INPUT1 <- df %>% filter(Re_model == i, input == input_1) %>% pull(Re_estimate)
        INPUT2 <- df %>% filter(Re_model == "CalCAT ensemble", input == input_2) %>% pull(Re_estimate)
        if(length(INPUT1)==0 | length(INPUT2)==0){next}
        
        INPUT1 <- frank(INPUT1)
        INPUT2 <- frank(INPUT2)
        corr<-ccf(INPUT1, INPUT2, lag.max=20, 
                  main=paste0(unique(df$location),
                              ": ", i))
        
        lag_analysis_temp_calcat <- data.frame(location=unique(df$location),
                                               location_level=unique(df$location_level),
                                               stationary_variable = input_2,
                                               sliding_variable = input_1,
                                               Re_model= i,
                                               lag=corr$lag, 
                                               acf=corr$acf) 
        
        lag_analysis = lag_analysis %>% bind_rows(lag_analysis_temp_calcat)
        
      }
    }
  }
  
  return(lag_analysis)
}

# 12. Function to aid in formatting/aesthetic fine-tuning of final Table 6 of manuscript

max_lag_table = function(df, county, stationary, sliding){
  
  data = df %>% filter(County_address == county, stationary_variable == stationary, 
                       sliding_variable == sliding) 

  data_processed = data %>% 
    mutate(label_name = ifelse(label_name == "COUNTY", paste0(County_address), label_name)) %>%
    tidyr::pivot_wider(id_cols = c(label_name, County_address, stationary_variable, sliding_variable), 
                       names_from = Re_model, 
                       values_from = c(lag_range, maximum_acf_range))  %>% dplyr::select(-c(County_address))
  
  return(data_processed)
}

# 13. Create a 'geom_tile' plot for cross-correlation analyses (as in Figures S5 & S6 of the manuscript)

ccf_plot_all_counties = function(df = final_lag_county_only, county = NULL, stationary, sliding){
  
  if(stationary == "raw_cases") {
    subtitle_hjust = 0.5
    subtitle=expression(paste(R[ww]," vs. Sewershed-restricted ", R[cc]))
    title = county} else if(stationary == "CalCAT") {
      subtitle_hjust = 0.5
      subtitle=expression(paste(R[ww]," vs. CalCAT ensemble"))
      title = county
    }
  
  if(is.null(county)) {
    data = df %>% 
      filter(stationary_variable == stationary, sliding_variable == sliding) %>% mutate(location_level = as.factor(location_level)) %>%
      arrange(location_level)
    
  } else{
    data = df %>% 
      filter(County_address == county, stationary_variable == stationary, sliding_variable == sliding) %>% mutate(location_level = as.factor(location_level)) %>%
      arrange(location_level)}
  
  p = data %>% ggplot(aes(label_name, lag, fill=acf)) + 
    geom_tile(color = "black") +
    geom_text(aes(label = round(acf, 2)), color = "black", size = 4.5) +
    #scale_fill_gradient(low=swatch()[3], high=swatch()[4]) +
    scale_fill_gradientn(colours = colorspace::divergingx_hcl(palette = "Spectral", n = 10, rev = T)) +   
    facet_grid(Re_model ~ location) +
    #coord_fixed(ratio = .2) + 
    #facet_wrap(~Re_model)+
    #labs(title = paste0(county, ": ", sliding, " vs. ", stationary))+
    labs(title = title, subtitle = subtitle)+
    ylab("Lag")+
    theme(
      legend.title = element_text(size=13),
      legend.key.size = unit(0.3, 'cm'),
      legend.text = element_text(size=13),
      plot.title = element_text(size = 19, hjust = 0.41),
      plot.subtitle = element_text(size = 17, hjust = subtitle_hjust),
      axis.text.x = element_blank(), #element_text(size = 15, angle = 45, vjust = 1, hjust=1),
      axis.text.y = element_text(size=13), 
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=13),
      strip.text.x = element_text(size = 13, hjust = 0.5),  #vjust = 1
      strip.text.y = element_text(size = 13),
      legend.position = "right",
      legend.box = "vertical")
  
  p <- p + guides(fill=guide_legend(title="CC"))
  
  
  return(p)
}

# 14. Produce data frame outputting confusion matrix metrics of interest (for Tables 4 & 5 of the manuscript)

# Sample input variables for function testing
# df = all_re_df_confusionmatrix
# loc_level = "county" 
# loc = "Stanislaus"
# actuals_model = "CalCAT ensemble"
# preds_model_vector = c("EpiEstim","Wallinga&Teunis")
# preds_model_input = "root_transform"
# j = 1
# j = 2

confusionmatrix_metrics = function(df = all_re_df_confusionmatrix, 
                                   loc_level = "county", 
                                   loc, 
                                   actuals_model = "CalCAT ensemble", 
                                   preds_model_vector = c("EpiEstim","Wallinga&Teunis"), 
                                   preds_model_input = "root_transform"){
  
  test = df %>% filter(location_level == loc_level, location == loc)
  actual = test %>% filter(Re_model %in% actuals_model) %>% rename(actual = category)
  
  for(j in 1:length(preds_model_vector)){
    i = preds_model_vector[j]
    pred = test %>% filter(Re_model == i, input == preds_model_input) %>% rename(pred = category)
    test2_pre = full_join(actual, pred, by = c("sample_date", "location")) %>% dplyr::select(sample_date, location, actual, pred)
    levels = c("<0.7", "0.7-0.9", "0.9-1.1", "1.1-1.3",">1.3")
    test2 = test2_pre %>% mutate(actual = factor(actual, levels = levels), pred = factor(pred, levels = levels))
    cm = caret::confusionMatrix(test2$pred, test2$actual)
    cmdf = data.frame(cm$byClass) %>% dplyr::select(-c(5:11)) %>% mutate(Class = str_remove(rownames(cm$byClass), "Class: "), .before = Sensitivity) %>%
      filter(Class %notin% c("<0.7", ">1.3")) %>% 
      mutate('Overall Accuracy' = cm$overall['Accuracy']) %>%
      mutate_if(is.numeric,
                round,
                digits = 2)
    
    rownames(cmdf) = NULL
    cmdf$location = unique(test2_pre$location)
    cmdf$Re_model = unique(pred$Re_model)
    cmdf = cmdf %>% rowwise() %>%
      mutate(`Sensitivity, Specificity, Pos. Pred. Value, Neg.Pred. Value`=paste(Sensitivity,Specificity,Pos.Pred.Value,Neg.Pred.Value, sep = ", "))
    cmdf = cmdf %>%
      tidyr::pivot_wider(id_cols = c(location, Class), 
                         names_from = Re_model, 
                         names_sep = ": ",
                         values_from =c(`Sensitivity, Specificity, Pos. Pred. Value, Neg.Pred. Value`, `Overall Accuracy`))
    
    
    if(j == 1) {cm_final = cmdf} else {cm_final = full_join(cm_final, cmdf)}
  }
  return(cm_final)
}

# 15. Produce data frame outputting confusion matrix-based frequency of agreement between wastewater- and CalCAT-based effective reproduction numbers for the counties of interest
confusionmatrix_as_df = function(df = all_re_df_confusionmatrix, 
                                 loc_level = "county", 
                                 loc, 
                                 actuals_model = "CalCAT ensemble", 
                                 preds_model_vector = c("EpiEstim","Wallinga&Teunis"), 
                                 preds_model_input = "root_transform"){
  
  test = df %>% filter(location_level == loc_level, location == loc)
  actual = test %>% filter(Re_model %in% actuals_model) %>% rename(actual = category)
  
  for(j in 1:length(preds_model_vector)){
    i = preds_model_vector[j]
    pred = test %>% filter(Re_model == i, input == preds_model_input) %>% rename(pred = category)
    test2_pre = full_join(actual, pred, by = c("sample_date", "location")) %>% dplyr::select(sample_date, location, actual, pred)
    levels = c("<0.7", "0.7-0.9", "0.9-1.1", "1.1-1.3",">1.3")
    test2 = test2_pre %>% mutate(actual = factor(actual, levels = levels), pred = factor(pred, levels = levels))
    cm = caret::confusionMatrix(test2$pred, test2$actual)
    cmtable = cm$table
    cmdf = as.data.frame(cmtable)
    cmdf$location = unique(test2_pre$location)
    cmdf$Re_model = unique(pred$Re_model)
    if(j == 1) {cm_final = cmdf} else {cm_final = rbind(cm_final, cmdf)}
  }
  return(cm_final)
}

# 16. Calculate mean absolute error and Spearman's rank correlations to compare wastewater- and case-based effective reproduction numbers

# Sample input parameters to test function 
# df = all_re_df %>% filter(County_address=="Yolo") # county-level example for CalCAT section
# models = c("EpiEstim", "Wallinga&Teunis")
# raw_cases_comparison = F
# calcat_comparison = T
# inputs_vec = c("root_transform")
# i = "Wallinga&Teunis"
# j = 1
# k = 2

compare_Re <- function(df, 
                       models = c("EpiEstim", "Wallinga&Teunis"), 
                       raw_cases_comparison = F, # Does the user want to compute cross-correlation coefficients between wastewater- and sewershed-restricted, case-based effective reproduction numbers?
                       calcat_comparison = T, # Only T if user is conducting county-level analyses; for sewershed-level, this should be F
                       inputs_vec) {
  
  comparison_analysis = data.frame(location=character(),
                                   location_level=character(),
                                   Re_model=character(),
                                   actual_x_variable = character(),
                                   predicted_y_variable = character(),
                                   mae = numeric(),
                                   spearman_rho=numeric(), 
                                   spearman_pvalue=character()) 
  
  for(i in models) { # Instantaneous or cohort
    if(raw_cases_comparison == T){
      # Use the below 'for' loop to calculate metrics comparing two time series of interest (in the manuscript, this loop was utilized to compare wastewater-based and sewershed-restricted, case-based effective reproduction number estimates)
      for(j in 1:ncol(combn(inputs_vec,2))){
        inputs = combn(inputs_vec,2)[,j] %>% sort()
        input_1 = inputs[1]
        input_2 = inputs[2]
        INPUT1 <- df %>% filter(Re_model == i, input == input_1) %>% pull(Re_estimate)
        INPUT2 <- df %>% filter(Re_model == i, input == input_2) %>% pull(Re_estimate)
        
        if( length(INPUT1)==0 | length(INPUT2)==0){next}
        if( length(INPUT1)!=length(INPUT2)){
          
          INPUT1 <- df %>% filter(Re_model == i, input == input_1) 
          RANGE1 = range(INPUT1$sample_date)
          
          INPUT2 <- df %>% filter(Re_model == i, input == input_2) 
          RANGE2 = range(INPUT2$sample_date)
          
          threshold_date = max(c(min(RANGE1), min(RANGE2)))
          
          INPUT1 <- INPUT1 %>% filter(sample_date >= threshold_date) %>% pull(Re_estimate)
          INPUT2 <- INPUT2 %>% filter(sample_date >= threshold_date) %>% pull(Re_estimate)
        }
        
        mae <- Metrics::mae(INPUT1, INPUT2) # Mean absolute error
        spearman <- cor.test(INPUT1, INPUT2, method = 'spearman', exact=FALSE) # Spearman's correlation
        comparison_analysis_temp <- data.frame(location=unique(df$location),
                                               location_level=unique(df$location_level),
                                               Re_model=i,
                                               actual_x_variable = input_1,
                                               predicted_y_variable = input_2,
                                               mae = mae, 
                                               spearman_rho=spearman$estimate, 
                                               spearman_pvalue=spearman$p.value)
        
        # Format to improve data frame aesthetics
        rownames(comparison_analysis_temp) <- NULL
        comparison_analysis_temp$spearman_pvalue <- formatC(comparison_analysis_temp$spearman_pvalue, format = "e", digits = 1)
        
        comparison_analysis_temp <- comparison_analysis_temp %>%
          mutate('spearman_pvalue' = as.character(comparison_analysis_temp$'spearman_pvalue')) %>%
          mutate(across(where(is.numeric), function(x) round(x, digits = 3)))
        
        comparison_analysis = comparison_analysis %>% bind_rows(comparison_analysis_temp)
      }}
    # Repeat above process, but now comparing wastewater-based estimates to CalCAT-based estimates
    # The below 'if' statement will also calculate cross-correlation coefficients between the sewershed-restricted case- and CalCAT-based effective reproduction numbers as a secondary product
    if(calcat_comparison == T){  # Only T if user is conducting county-level analyses; for sewershed-level, this should be F
      for(k in 1:length(inputs_vec)){
        input_1 = "CalCAT"
        input_2 = inputs_vec[k]
        
        INPUT1 <- df %>% filter(Re_model == "CalCAT ensemble", input == input_1) %>% pull(Re_estimate)
        INPUT2 <- df %>% filter(Re_model == i, input == input_2) %>% pull(Re_estimate)
        
        if(length(INPUT1)==0 | length(INPUT2)==0){next}
        
        mae <- Metrics::mae(INPUT1, INPUT2)
        spearman <- cor.test(INPUT1, INPUT2, method = 'spearman', exact=FALSE)
        
        comparison_analysis_temp_calcat <- data.frame(location=unique(df$location),
                                                      location_level=unique(df$location_level),
                                                      Re_model=i,
                                                      actual_x_variable = input_1,
                                                      predicted_y_variable = input_2,
                                                      mae = mae, 
                                                      spearman_rho=spearman$estimate, 
                                                      spearman_pvalue=spearman$p.value)
        
        rownames(comparison_analysis_temp_calcat) <- NULL
        
        comparison_analysis_temp_calcat$spearman_pvalue <- formatC(comparison_analysis_temp_calcat$spearman_pvalue, format = "e", digits = 1)
        
        comparison_analysis_temp_calcat <- comparison_analysis_temp_calcat %>% mutate('spearman_pvalue' = as.character(comparison_analysis_temp_calcat$'spearman_pvalue')) %>%
          mutate(across(where(is.numeric), function(x) round(x, digits = 3)))
        
        comparison_analysis = comparison_analysis %>% bind_rows(comparison_analysis_temp_calcat)
      }
    }
  }
  return(comparison_analysis)
}

# 17. Function to aid in recreation of Table 3 within manuscript

# Sample input variables to test function 
# df = final_comparison_restricted
# county = "Santa Clara"
# actual_x = "CalCAT"
# predicted_y = "root_transform"

comparison_table = function(df, county, actual_x, predicted_y){
  
  data = df %>% filter(County_address == county, actual_x_variable == actual_x, predicted_y_variable == predicted_y) 
  
  data_processed = data %>% mutate(label_name = ifelse(label_name == "COUNTY", paste0(County_address), label_name)) %>%
    tidyr::pivot_wider(id_cols = c(label_name, County_address, actual_x_variable, predicted_y_variable), 
                       names_from = Re_model, 
                       names_sep = ", ",
                       values_from = c(mae, spearman_rho, spearman_pvalue))  %>% dplyr::select(-c(County_address))
  
  return(data_processed)
}

