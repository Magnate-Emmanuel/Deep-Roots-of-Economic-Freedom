/*****************************************************************************/
/*** DESCRIPTION: REPLICATION CODE FOR "DIVERSITY AND ECONOMIC FREEDOM"            ***/
/*** AUTHORS:      ***/
/*****************************************************************************/

/***********************************************************************************************************************************************************************/
/* TABLE: Population Diversity and Economic Freedom across Countries - Robustness to Accounting for Population Diversity as a Generated Regressor */
/***********************************************************************************************************************************************************************/

# delimit 
clear all 

set more off 

/*******************************************************************************/
/* Define programs for the two-step bootstrapped estimation of standard errors */
/*******************************************************************************/

capture program drop getbeta 
program define getbeta, rclass 
  version 14.2 
  syntax varlist if [, iv(varname numeric)] 

  tokenize `varlist' 
  local confvar `1' 
  local pdivvar `2' 
  macro shift 2 
  local ctlvars `*' 

  marksample touse 
  if "`iv'" == "" { 
    qui regress `confvar' `pdivvar' `ctlvars' if `touse' 
  } 
  else { 
    qui ivregress 2sls `confvar' `ctlvars' (`pdivvar' = `iv') if `touse' 
  } 
  qui sum `confvar' if e(sample) == 1 
  local avgconf = r(mean) 
  _pctile `pdivvar' if e(sample) == 1, percentiles(10 90) 
  local delta = r(r2) - r(r1) 
  qui lincom ((1 + `avgconf') * `delta' * `pdivvar') 
  matrix B = (e(b), r(estimate)) 

  return matrix beta = B 
  return scalar r2_a = e(r2_a) 
  return scalar Nobs = e(N) 
end 

capture program drop tseboot 
program define tseboot, eclass 
  version 14.2 
  syntax varlist if [, iv(varname numeric)] 

  tokenize `varlist' 
  local confvar `1' 
  local pdivvar `2' 
  macro shift 2 
  local ctlvars `*' 

  use "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Project\Diversity and conflict\Replication files\Replication\country\data\hgdp.dta", clear 
  bsample 1, strata(country)
  regress adiv mdist 

  use "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\CROSS\Data\fraser_1970.dta", clear 
  marksample touse 
  bsample if `touse' 
  replace avg_pdiv_aa = _b[_cons] + _b[mdist] * mdist_addis_aa 
  if "`iv'" == "" { 
    regress `confvar' `pdivvar' `ctlvars'  
  } 
  else { 
    ivregress 2sls `confvar' `ctlvars' (`pdivvar' = `iv') 
  } 
  sum `confvar' 
  local avgconf = r(mean) 
  _pctile `pdivvar', percentiles(10 90) 
  local delta = r(r2) - r(r1) 
  lincom ((1 + `avgconf') * `delta' * `pdivvar') 

  ereturn scalar meff = r(estimate) 
end 

/*******************/
/* Open a log file */
/*******************/

capture log close 

/*****************************************/
/* Load the AAGK cross-sectional dataset */
/*****************************************/

use "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\CROSS\Data\fraser_1970.dta" , clear 




/********************************************************/
/* Create global macro lists of the relevant covariates */
/********************************************************/

global ylist "a_overall_index"
global xlist "avg_pdiv_aa"
global pd_vars "pdiv_aa_w pdiv_aa_b"
global ivlist "avg_mdist_addis"
global geor " avg_abslat avg_distcr avg_ruggavg avg_suitavg avg_suitrng avg_elevavg avg_elevrng avg_island "
global colhist " avg_iscol_gbr avg_iscol_fra" 
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
qui egen smpl_flag1 = rowmiss($ylist $xlist $geor $climate $ethfrac $ethpolr) 
qui replace smpl_flag1 = (smpl_flag1 == 0) 

qui egen smpl_flag2 = rowmiss($ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $inst $legalor $colhist $oilprod) 
qui replace smpl_flag2 = (smpl_flag2 == 0)

/**************************************************/
/* Save the data changes to a temporary data file */
/**************************************************/

save "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\CROSS\Data\fraser_1970_tmp.dta", replace 

/***********************/
/* Run the regressions */
/***********************/

/* -------- */
/* COLUMN 1 */
/* -------- */
use "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\CROSS\Data\fraser_1970_tmp.dta", clear 

/* Obtain the relevant point estimates */
getbeta $ylist $xlist if smpl_flag1 == 1 
matrix beta = r(beta) 

/* Store the adjusted R2 and the number of observations for subsequent use */
scalar r2_a = r(r2_a) 
scalar Nobs = r(Nobs) 

/* Obtain the bootstrapped standard error estimates */
simulate _b meff = e(meff), reps(1000) seed(12345): tseboot $ylist $xlist if smpl_flag1 == 1 
eststo col1: bstat, stat(beta) 

/* Obtain the adjusted R2 and the number of observations */
estadd scalar adjr2 = r2_a : col1 
estadd scalar N_obs = Nobs : col1 

/* -------- */
/* COLUMN 2 */
/* -------- */
use "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\CROSS\Data\fraser_1970_tmp.dta", clear 

/* Obtain the relevant point estimates */
getbeta $ylist $xlist $geor if smpl_flag1 == 1
matrix beta = r(beta) 

/* Store the adjusted R2 and the number of observations for subsequent use */
scalar r2_a = r(r2_a) 
scalar Nobs = r(Nobs) 

/* Obtain the bootstrapped standard error estimates */
simulate _b meff = e(meff), reps(1000) seed(12345): tseboot $ylist $xlist $geor if smpl_flag1 == 1
eststo col2: bstat, stat(beta) 

/* Obtain the adjusted R2 and the number of observations */
estadd scalar adjr2 = r2_a : col2 
estadd scalar N_obs = Nobs : col2 

/* -------- */
/* COLUMN 3 */
/* -------- */
use "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\CROSS\Data\fraser_1970_tmp.dta", clear 

/* Obtain the relevant point estimates */
getbeta $ylist $xlist $geor $climate if smpl_flag1 == 1
matrix beta = r(beta) 

/* Store the adjusted R2 and the number of observations for subsequent use */
scalar r2_a = r(r2_a) 
scalar Nobs = r(Nobs) 

/* Obtain the bootstrapped standard error estimates */
simulate _b meff = e(meff), reps(1000) seed(12345): tseboot $ylist $xlist $geor $climate if smpl_flag1 == 1
eststo col3: bstat, stat(beta) 

/* Obtain the adjusted R2 and the number of observations */
estadd scalar adjr2 = r2_a : col3 
estadd scalar N_obs = Nobs : col3 

/* -------- */
/* COLUMN 4 */
/* -------- */
use "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\CROSS\Data\fraser_1970_tmp.dta", clear 

/* Obtain the relevant point estimates */
getbeta $ylist $xlist $geor $climate $ethfrac $ethpolr if smpl_flag1 == 1
matrix beta = r(beta) 

/* Store the adjusted R2 and the number of observations for subsequent use */
scalar r2_a = r(r2_a) 
scalar Nobs = r(Nobs) 

/* Obtain the bootstrapped standard error estimates */
simulate _b meff = e(meff), reps(1000) seed(12345): tseboot $ylist $xlist $geor $climate $ethfrac $ethpolr if smpl_flag1 == 1
eststo col4: bstat, stat(beta) 

/* Obtain the adjusted R2 and the number of observations */
estadd scalar adjr2 = r2_a : col4 
estadd scalar N_obs = Nobs : col4 

/* -------- */
/* COLUMN 5 */
/* -------- */
use "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\CROSS\Data\fraser_1970_tmp.dta", clear 

/* Obtain the relevant point estimates */
getbeta $ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod if smpl_flag1 == 1
matrix beta = r(beta) 

/* Store the adjusted R2 and the number of observations for subsequent use */
scalar r2_a = r(r2_a) 
scalar Nobs = r(Nobs) 

/* Obtain the bootstrapped standard error estimates */
simulate _b meff = e(meff), reps(1000) seed(12345): tseboot $ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod if smpl_flag1 == 1
eststo col5: bstat, stat(beta) 

/* Obtain the adjusted R2 and the number of observations */
estadd scalar adjr2 = r2_a : col5 
estadd scalar N_obs = Nobs : col5 

/* -------- */
/* COLUMN 6 */
/* -------- */
use "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\CROSS\Data\fraser_1970_tmp.dta", clear 

/* Obtain the relevant point estimates */
getbeta $ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst if smpl_flag1 == 1
matrix beta = r(beta) 

/* Store the adjusted R2 and the number of observations for subsequent use */
scalar r2_a = r(r2_a) 
scalar Nobs = r(Nobs) 

/* Obtain the bootstrapped standard error estimates */
simulate _b meff = e(meff), reps(1000) seed(12345): tseboot $ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst if smpl_flag1 == 1
eststo col6: bstat, stat(beta) 

/* Obtain the adjusted R2 and the number of observations */
estadd scalar adjr2 = r2_a : col6 
estadd scalar N_obs = Nobs : col6 

/* -------- */
/* COLUMN 7 */
/* -------- */
use "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\CROSS\Data\fraser_1970_tmp.dta", clear 

/* Obtain the relevant point estimates */
getbeta $ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $legalor $colhist $contall if smpl_flag1 == 1
matrix beta = r(beta) 

/* Store the adjusted R2 and the number of observations for subsequent use */
scalar r2_a = r(r2_a) 
scalar Nobs = r(Nobs) 

/* Obtain the bootstrapped standard error estimates */
simulate _b meff = e(meff), reps(1000) seed(12345): tseboot $ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $legalor $colhist $contall if smpl_flag1 == 1
eststo col7: bstat, stat(beta) 

/* Obtain the adjusted R2 and the number of observations */
estadd scalar adjr2 = r2_a : col7 
estadd scalar N_obs = Nobs : col7 


/* -------- */
/* COLUMN 8 */
/* -------- */
use "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\CROSS\Data\fraser_1970_tmp.dta", clear 

/* Obtain the relevant point estimates */
getbeta $ylist $xlist $geor $climate if smpl_flag1 == 1, iv($ivlist) 
matrix beta = r(beta) 

/* Store the number of observations for subsequent use */
scalar Nobs = r(Nobs) 

/* Obtain the bootstrapped standard error estimates */
simulate _b meff = e(meff), reps(1000) seed(12345): tseboot $ylist $xlist $geor $climate  if smpl_flag1 == 1, iv($ivlist) 
eststo col8: bstat, stat(beta) 

/* Obtain the number of observations */
estadd scalar N_obs = Nobs : col8 

/* -------- */
/* COLUMN 9 */
/* -------- */
use "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\CROSS\Data\fraser_1970_tmp.dta", clear 

/* Obtain the relevant point estimates */
getbeta $ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $legalor $colhist $contall if smpl_flag1 == 1, iv($ivlist) 
matrix beta = r(beta) 

/* Store the number of observations for subsequent use */
scalar Nobs = r(Nobs) 

/* Obtain the bootstrapped standard error estimates */
simulate _b meff = e(meff), reps(1000) seed(12345): tseboot $ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $legalor $colhist $contall  if smpl_flag2 == 1, iv($ivlist) 
eststo col9: bstat, stat(beta) 

/* Obtain the number of observations */
estadd scalar N_obs = Nobs : col9 


/**************************************/
/* Print the results to a LaTeX table */
/**************************************/
           
estout col1 col2 col3 col4 col5 col6 col7 col8 col9 using "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\CROSS\Results\xtab_boots.tex", style(tex) replace   cells(b(star fmt(a2)) se(fmt(a2) par)) keep(_b_avg_pdiv_aa, relax)   indicate("Continent dummies=_b_avg_africa _b_avg_asia _b_avg_namerica _b_avg_samerica _b_avg_oceania"            "Controls for geography=_b_avg_abslat _b_avg_distcr _b_avg_ruggavg _b_avg_suitavg _b_avg_suitrng _b_avg_elevavg _b_avg_elevrng _b_avg_island"  "Controls for climate= _b_avg_ecofrac _b_avg_ecopol _b_avg_tmp_cru_v401 _b_avg_pre_cru_v401 _b_avg_tmp_vol_cru_v401 _b_avg_pre_vol_cru_v401" "Controls for ethnic diversity=_b_avg_efrac _b_avg_des_pol15"  "Controls for Colonial history= _b_avg_iscol_gbr " "Controls for Legal Origin=_b_avg_legor_uk _b_avg_legor_fr" "Controls for institutions= _b_avg_xconst_p4_v17 _b_avg_ddemoc_p4_v17 _b_avg_dautoc_p4_v17"            "Controls for oil, population, and income=_b_avg_anypetroleum_pet _b_avg_lnpop_wdi _b_avg_lngdppc_cu_usd_wdi", labels("\$\times\$" " "))   stats(N_obs adjr2, fmt(%9.0f a2) labels("Observations" "Adjusted \$R^2\$") layout(@ @))   varwidth(44) msign("\$-\$") nolabel   prehead("\begin{tabular*}{645pt}{@{\extracolsep{\fill}}lccccccccc}"           "\toprule"           "Cross-country sample:&\multicolumn{5}{c}{Global}&\multicolumn{2}{c}{Old World}&\multicolumn{2}{c}{Global}\\"           "\cmidrule(r){2-6}\cmidrule(lr){7-8}\cmidrule(l){9-10}")   numbers mlabels("OLS" "OLS" "OLS" "OLS" "OLS" "OLS" "OLS" "2SLS" "2SLS") collabels(none)   posthead("\midrule"            "&\multicolumn{9}{c}{Log number of new PRIO25 civil conflict onsets per year, 1960--2017}\\"            "\cmidrule{2-10}")   varlabels(_b_pdiv_aa "Population diversity (ancestry adjusted)", elist(_b_pdiv_aa \addlinespace))   prefoot("\midrule") postfoot("\addlinespace")


 estout col1 col2 col3 col4 col5 col6 col7 col8 col9 using "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\CROSS\Results\xtab_boots.tex", style(tex) append   cells(b(nostar fmt(a2)) se(fmt(a2) par)) keep(_eq2_meff, relax)   varwidth(44) msign("\$-\$") nolabel   mlabels(none) collabels(none)   varlabels(_eq2_meff "Effect of 10th--90th \%ile move in diversity")   postfoot("\bottomrule\addlinespace" "\end{tabular*}") 
 
 
/****************************************************************************************************************/
/* Erase the temporary data file from disk, clean-up stored estimates from memory, close the log file, and exit */
/****************************************************************************************************************/

erase "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\CROSS\Data\fraser_1970_tmp.dta"
est clear 
log close 

exit 
