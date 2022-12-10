/*****************************************************************************/
/*** DESCRIPTION: REPLICATION CODE FOR "DIVERSITY AND ECONOMIC FREEDOM"            ***/
/*** AUTHORS:     ***/
/*****************************************************************************/

/***********************************************************************************************************************************************************/
/* TABLE: Population Diversity and Economic Freedom across Countries - Robustness to Accounting for Spatial Autocorrelation in Errors */
/***********************************************************************************************************************************************************/

# delimit 
clear all 

set more off 


/**********************************************************/
/* Define wrapper programs for running Conley regressions */
/**********************************************************/

/* Define a wrapper program around "x_ols.ado" that implements Conley's (1999) OLS estimation with spatial standard errors */
capture program drop conleyreg 
program define conleyreg, eclass 
  version 14.2 
  syntax varlist [if] 

  tokenize `varlist' 
  local lhs_var `1' 
  local rhs_var `2' 
  macro shift    2  
  local ctl_lst `*' 

  local rhs_lst `"`rhs_var' `ctl_lst' _cons"' 
  local reg_num : list sizeof rhs_lst 

  quietly { 
    preserve 
    keep `if' 
    x_ols conley_coord1 conley_coord2 conley_cutoff1 conley_cutoff2 `lhs_var' `rhs_var' `ctl_lst' conley_const, xreg(`reg_num') coord(2) 
    restore 
  } 

  tempname b V 
  matrix `b' = e(b) 
  matrix `V' = cov_dep 

  matrix drop cov_dep 
  scalar drop fix 

  matrix colnames `b' = `rhs_lst' 
  matrix rownames `b' = `lhs_var' 
  matrix colnames `V' = `rhs_lst' 
  matrix rownames `V' = `rhs_lst' 

  local N    = e(N) 
  local r2   = e(r2) 
  local r2_a = e(r2_a) 

  ereturn post `b' `V' 

  ereturn local cmd     `"conleyreg"' 
  ereturn local cmdline `"conleyreg `varlist'"' 
  ereturn local coords  `"conley_coord1 conley_coord2"' 
  ereturn local cutoffs `"conley_cutoff1 conley_cutoff2"' 
  ereturn local depvar  `"`lhs_var'"' 
  ereturn local model   `"conley ols"' 
  ereturn local title   `"Conley OLS regression"' 
  ereturn local vcetype `"Conley"' 

  ereturn scalar N    = `N' 
  ereturn scalar r2   = `r2' 
  ereturn scalar r2_a = `r2_a' 

  ereturn display 
end 

/* Define a wrapper program around "x_gmm.ado" that implements Conley's (1999) spatial GMM estimation */
capture program drop conleygmm 
program define conleygmm, eclass 
  version 14.2 
  syntax varlist [if], iv(varlist) 

  tokenize `varlist' 
  local lhs_var `1' 
  local rhs_var `2' 
  macro shift    2  
  local ctl_lst `*' 

  local rhs_lst `"`rhs_var' `ctl_lst' _cons"' 
  local reg_num : list sizeof rhs_lst 
  local ins_lst `"`iv' `ctl_lst' _cons"' 
  local ins_num : list sizeof ins_lst 

  quietly { 
    preserve 
    keep `if' 
    x_gmm conley_coord1 conley_coord2 conley_cutoff1 conley_cutoff2 `lhs_var' `rhs_var' `ctl_lst' conley_const `iv' `ctl_lst' conley_const, xreg(`reg_num') inst(`ins_num') coord(2) 
    restore 
  } 

  tempname b V 
  matrix `b' = betagmm' 
  matrix `V' = cov_dep 

  matrix drop J cov_dep betagmm bgmm1 cov2SLS b2SLS 

  matrix colnames `b' = `rhs_lst' 
  matrix rownames `b' = `lhs_var' 
  matrix colnames `V' = `rhs_lst' 
  matrix rownames `V' = `rhs_lst' 

  qui count `if' 
  local N = r(N) 

  ereturn post `b' `V' 

  ereturn local cmd       `"conleygmm"' 
  ereturn local cmdline   `"conleygmm `varlist', iv(`iv')"' 
  ereturn local coords    `"conley_coord1 conley_coord2"' 
  ereturn local cutoffs   `"conley_cutoff1 conley_cutoff2"' 
  ereturn local depvar    `"`lhs_var'"' 
  ereturn local instd     `"`rhs_var'"' 
  ereturn local insts     `"`ctl_lst' `iv'"' 
  ereturn local title     `"Conley GMM regression"' 
  ereturn local vcetype   `"Conley"' 
  ereturn local estimator `"conley gmm"' 
  ereturn local exogr     `"`ctl_lst'"' 

  ereturn scalar N = `N' 

  ereturn display 
end 

/*******************/
/* Open a log file */
/*******************/

capture log close

/*****************************************/
/* Load the AAGK cross-sectional dataset */
/*****************************************/


/****************************************************************************/
/* Assign appropriate labels to the RHS variables for the regression output */
/****************************************************************************/

/*Economic Freedom Variables*/     
label variable a_overall_index "Economic Freedom Fraser Index"
label variable a_SizeofGovernment "Size of Government"
label variable a_LegalSystemandPropertyRights "Legal System and Property Rights"
label variable a_SoundMoney "Sound Money"
label variable a_FreedomtoTradeInternationally "Freedom to Trade Internationally"
label variable a_Regulation "Regulation"

/* Diversity variables */
label variable avg_pdiv_aa   "Population diversity (ancestry adjusted)" 
label variable pdiv_aa_w "Within-group population diversity" 
label variable pdiv_aa_b "Between-group population diversity" 
label variable avg_efrac     "Ethnic fractionalization" 
label variable avg_des_pol15 "Ethnolinguistic polarization" 
label variable avg_pdiv_aa_sq "Population diversity (ancestry adjusted) Squared"

       
/* Geographical variables */
label variable avg_abslat      "Absolute latitude" 
label variable avg_ruggavg     "Ruggedness" 
label variable avg_elevavg     "Mean elevation" 
label variable avg_elevrng     "Range of elevation" 
label variable avg_suitavg     "Mean land suitability" 
label variable avg_suitrng     "Range of land suitability" 
label variable avg_distcr      "Distance to nearest waterway" 
label variable avg_island      "Island nation dummy" 
label variable avg_mdist_addis "Migratory distance from East Africa (in 10,000 km)" 
label variable avg_mdist_addis_sq "Migratory distance from East Africa (in 10,000 km) Squared"
   
/* Deep institutional variables */
label variable avg_iscol_gbr "Ever a U.K. colony dummy" 
label variable avg_iscol_fra "Ever a French colony dummy" 
label variable avg_iscol_oth "Ever a non-U.K./non-French colony dummy" 
label variable avg_legor_uk    "British legal origin dummy" 
label variable avg_legor_fr    "French legal origin dummy" 
  
/* Contemporary institutional variables */
label variable avg_xconst_p4_v17 "Executive constraints, 1960--2017 average" 
label variable avg_ddemoc_p4_v17  "Fraction of years under democracy, 1960--2017" 
label variable avg_dautoc_p4_v17  "Fraction of years under autocracy, 1960--2017" 

/*****************************************/
/* Oil, population, and income variables */
/*****************************************/ 

label variable avg_anypetroleum_pet   "Oil or gas reserve discovery" 
label variable avg_lnpop_wdi          "Log population, 1960--2017 average" 
label variable avg_lngdppc_cu_usd_wdi "Log GDP per capita, 1960--2017 average" 

/****************************************/
/* Corruption Variables                 */
/****************************************/    
label variable avg_v2exbribe "Executive bribery and corrupt exchanges"
label variable avg_v2excrptps "Public sector corrupt exchanges"
label variable avg_v2lgcrrpt "Legislature corrupt activities"
label variable avg_v2jucorrdc "Judicial corruption decision"
label variable avg_v2x_corr "Political corruption index"

/****************************************/
/*    Continent Dummies                 */
/****************************************/  
label variable avg_africa "Africa Dummy"
label variable avg_asia "Asia Dummy"
label variable avg_namerica "North America Dummy" 
label variable avg_samerica "South America Dummy"
label variable avg_oceania "Oceania Dummy"
label variable avg_europe "Europe Dummy"

/********************************************************/
/* Create global macro lists of the relevant covariates */
/********************************************************/

global ethfrac "avg_efrac " 
global ethpolr "avg_des_pol15" 
global ylist "a_overall_index"
global xlist "avg_pdiv_aa"
global pd_vars "pdiv_aa_w pdiv_aa_b"
global ivlist "avg_mdist_addis"
global geor " avg_abslat avg_distcr avg_ruggavg avg_suitavg avg_suitrng avg_elevavg avg_elevrng avg_island "
global colhist " avg_iscol_gbr avg_iscol_fra avg_iscol_oth" 
global legalor " avg_legor_uk avg_legor_fr"  
global popsize " avg_lnpop_wdi" 
global ypercap " avg_lngdppc_cu_usd_wdi "  
global contall " avg_africa avg_asia avg_namerica avg_samerica avg_oceania" 
global contold "avg_africa avg_asia"
global inst " avg_xconst_p4_v17 avg_ddemoc_p4_v17 avg_dautoc_p4_v17"
global oilprod " avg_anypetroleum_pet" 
global climate "avg_ecofrac avg_ecopol avg_tmp_cru_v401 avg_pre_cru_v401 avg_tmp_vol_cru_v401 avg_pre_vol_cru_v401"
global deepdet "avg_ln_yst_aa avg_ln_statehist avg_ln_origtime avg_ln_frontdist1500" 
global ethineq "avg_grg_a_lum00pc avg_v1_a_lum00pc"

/******************************************************/
/* Generate the sample indicators for the regressions */
/******************************************************/

qui egen smpl_flag1 = rowmiss($ylist $xlist ${geog} $climate ${ethfrac} ${ethpolr}) 
qui replace smpl_flag1 = (smpl_flag1 == 0) 

qui egen smpl_flag2 = rowmiss($ylist $xlist ${geog} $climate ${ethfrac} ${ethpolr} ${colhist} ${legalor} ${inst} ${oilprod} ${popsize} ${ypercap}) 
qui replace smpl_flag2 = (smpl_flag2 == 0) 


/***********************/
/* Run the regressions */
/***********************/

/* -------- */
/* COLUMN 1 */
/* -------- */
eststo col1: conleyreg $ylist $xlist if smpl_flag1 == 1 

/* Obtain the adjusted R2 */
estadd scalar adjr2 = e(r2_a) : col1 

/* -------- */
/* COLUMN 2 */
/* -------- */
eststo col2: conleyreg $ylist $xlist ${geor} if smpl_flag1 == 1 

/* Obtain the adjusted R2 */
estadd scalar adjr2 = e(r2_a) : col2 

/* -------- */
/* COLUMN 3 */
/* -------- */
eststo col3: conleyreg $ylist $xlist ${geor} ${climate} if smpl_flag1 == 1 

/* Obtain the adjusted R2 */
estadd scalar adjr2 = e(r2_a) : col3 

/* -------- */
/* COLUMN 4 */
/* -------- */
eststo col4: conleyreg $ylist $xlist ${geor} ${climate} ${ethfrac} ${ethpolr}  if smpl_flag1 == 1 

/* Obtain the adjusted R2 */
estadd scalar adjr2 = e(r2_a) : col4 

/* -------- */
/* COLUMN 5 */
/* -------- */
eststo col5: conleyreg  $ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod if smpl_flag2 == 1 

/* Obtain the adjusted R2 */
estadd scalar adjr2 = e(r2_a) : col5 

/* -------- */
/* COLUMN 6 */
/* -------- */
eststo col6: conleyreg $ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst if smpl_flag2 == 1 

/* Obtain the adjusted R2 */
estadd scalar adjr2 = e(r2_a) : col6 

/* -------- */
/* COLUMN 7 */
/* -------- */
eststo col7: conleyreg $ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $legalor $contall if smpl_flag2 == 1 

/* Obtain the adjusted R2 */
estadd scalar adjr2 = e(r2_a) : col7 

/* -------- */
/* COLUMN 8 */
/* -------- */
eststo col8: conleygmm $ylist $xlist $geor $climate if smpl_flag1 == 1, iv($ivlist) 

/* -------- */
/* COLUMN 9 */
/* -------- */
eststo col9: conleygmm $ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $legalor $contall if smpl_flag2 == 1, iv($ivlist)


/**************************************/
/* Print the results to a LaTeX table */
/**************************************/

estout col1 col2 col3 col4 col5 col6 col7 col8 col9 using "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\CROSS\Results\stab_conley.tex", style(tex) replace cells(b(star fmt(a2)) se(fmt(a2) par)) keep($xlist, relax)   indicate("Continent dummies=${contall}"            "Controls for geography=${geor}"      "Controls for climate=${climate}"      "Controls for ethnic diversity=${ethfrac} ${ethpolr}"            "Controls for institutions=  ${inst}"     "Controls for Legal Origins= ${legalor}"       "Controls for oil, population, and income=${oilprod} ${popsize} ${ypercap}", labels("\$\times\$" " "))   stats(N adjr2, fmt(%9.0f a2) labels("Observations" "Adjusted \$R^2\$") layout(@ @))   varwidth(40) msign("\$-\$") label   prehead("\begin{tabular*}{610pt}{@{\extracolsep{\fill}}lccccccccccc}"           "\toprule"           "Cross-country sample:&\multicolumn{5}{c}{Global}&\multicolumn{2}{c}{Old World}&\multicolumn{2}{c}{Global}\\"           "\cmidrule(r){2-6}\cmidrule(lr){7-8}\cmidrule(l){9-10}")   numbers mlabels("OLS" "OLS" "OLS" "OLS" "OLS" "OLS" "OLS" "GMM" "GMM") collabels(none)   posthead("\midrule"            "&\multicolumn{9}{c}{Log number of new PRIO25 civil conflict onsets per year, 1960--2017}\\"            "\cmidrule{2-10}")   varlabels(, elist(pdiv_aa \addlinespace))   prefoot("\midrule") postfoot("\bottomrule\addlinespace" "\end{tabular*}") 

/***********************************************************/
/* Clean-up stored estimates, close the log file, and exit */
/***********************************************************/

est clear 
log close 

exit 


