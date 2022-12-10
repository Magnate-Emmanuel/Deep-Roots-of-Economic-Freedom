/*****************************************************************************/
/*** DESCRIPTION: REPLICATION CODE FOR "DIVERSITY AND ECONOMIC FREEDOM"            ***/
/*** AUTHORS:      ***/
/*****************************************************************************/

/********************************************************************************************************************************************/
/* TABLE: Population Diversity and Economic Freedom - Robustness to Accounting for Spatial Dependence */
/********************************************************************************************************************************************/

# delimit 
clear all 

set more off 

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

qui egen smpl_flag1 = rowmiss(a_overall_index avg_pdiv_aa $geor $climate $ethfrac $ethpolr) 
qui replace smpl_flag1 = (smpl_flag1 == 0) 

qui egen smpl_flag2 = rowmiss(a_overall_index avg_pdiv_aa $geor $climate $ethfrac $ethpolr $ypercap $popsize $inst $legalor $colhist $oilprod) 
qui replace smpl_flag2 = (smpl_flag2 == 0) 

/*****************************************************************************************************/
/* Generate the Mata inverse-distance spatial-weighting matrices for the relevant estimation samples */
/*****************************************************************************************************/

clear mata 

drop if cenlat == . | cenlong == . 
drop if avg_ddemoc_p4_v17 == .
drop if avg_xconst_p4_v17 == .
drop if avg_dautoc_p4_v17 == .

spmat idistance dobj1_global cenlong cenlat if smpl_flag1 == 1, id(cid) dfunction(dhaversine) normalize(minmax) dtruncate(10 10) replace 
spmat idistance dobj1_oworld cenlong cenlat if smpl_flag1 == 1 & oldw == 1, id(cid) dfunction(dhaversine) normalize(minmax) dtruncate(10 10) replace 

spmat idistance dobj2_global cenlong cenlat if smpl_flag2 == 1, id(cid) dfunction(dhaversine) normalize(minmax) dtruncate(10 10) replace 
spmat idistance dobj2_oworld cenlong cenlat if smpl_flag2 == 1 & oldw == 1, id(cid) dfunction(dhaversine) normalize(minmax) dtruncate(10 10) replace 


/***********************/
/* Run the regressions */
/***********************/

/* -------- */
/* COLUMN 1 */
/* -------- */
eststo col1: spreg gs2sls $ylist $xlist if smpl_flag1 == 1, id(cid) dlmat(dobj1_global) elmat(dobj1_global) heteroskedastic

/* Obtain the effect of a 10th to 90th percentile move in diversity */
qui sum $ylist if e(sample) == 1 
local avgconf = r(mean) 
_pctile $xlist if e(sample) == 1, percentiles(10 90) 
eststo col1_me: nlcom (me_1090: _b[$xlist] * (1 + `avgconf') * (`r(r2)' - `r(r1)')), post 

/* -------- */
/* COLUMN 2 */
/* -------- */
eststo col2: spreg gs2sls $ylist $xlist $geor if smpl_flag1 == 1, id(cid) dlmat(dobj1_global) elmat(dobj1_global) heteroskedastic

/* Obtain the effect of a 10th to 90th percentile move in diversity */
qui sum $ylist if e(sample) == 1 
local avgconf = r(mean) 
_pctile $xlist if e(sample) == 1, percentiles(10 90) 
eststo col2_me: nlcom (me_1090: _b[$xlist] * (1 + `avgconf') * (`r(r2)' - `r(r1)')), post

/* -------- */
/* COLUMN 3 */
/* -------- */
eststo col3: spreg gs2sls $ylist $xlist $geor $climate if smpl_flag1 == 1, id(cid) dlmat(dobj1_global) elmat(dobj1_global) heteroskedastic

/* Obtain the effect of a 10th to 90th percentile move in diversity */
qui sum $ylist if e(sample) == 1 
local avgconf = r(mean) 
_pctile $xlist if e(sample) == 1, percentiles(10 90) 
eststo col3_me: nlcom (me_1090: _b[$xlist] * (1 + `avgconf') * (`r(r2)' - `r(r1)')), post


/* -------- */
/* COLUMN 4 */
/* -------- */
eststo col4: spreg gs2sls $ylist $xlist $geor $climate $ethfrac $ethpolr if smpl_flag1 == 1, id(cid) dlmat(dobj1_global) elmat(dobj1_global) heteroskedastic

/* Obtain the effect of a 10th to 90th percentile move in diversity */
qui sum $ylist if e(sample) == 1 
local avgconf = r(mean) 
_pctile $xlist if e(sample) == 1, percentiles(10 90) 
eststo col4_me: nlcom (me_1090: _b[$xlist] * (1 + `avgconf') * (`r(r2)' - `r(r1)')), post

/* -------- */
/* COLUMN 5 */
/* -------- */
eststo col5: spreg gs2sls $ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod if smpl_flag1 == 1, id(cid) dlmat(dobj1_global) elmat(dobj1_global) heteroskedastic

/* Obtain the effect of a 10th to 90th percentile move in diversity */
qui sum $ylist if e(sample) == 1 
local avgconf = r(mean) 
_pctile $xlist if e(sample) == 1, percentiles(10 90) 
eststo col5_me: nlcom (me_1090: _b[$xlist] * (1 + `avgconf') * (`r(r2)' - `r(r1)')), post

/* -------- */
/* COLUMN 6 */
/* -------- */
eststo col6: spreg gs2sls $ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $legalor $colhist if smpl_flag1 == 1, id(cid) dlmat(dobj1_global) elmat(dobj1_global) heteroskedastic

/* Obtain the effect of a 10th to 90th percentile move in diversity */
qui sum $ylist if e(sample) == 1 
local avgconf = r(mean) 
_pctile $xlist if e(sample) == 1, percentiles(10 90) 
eststo col6_me: nlcom (me_1090: _b[$xlist] * (1 + `avgconf') * (`r(r2)' - `r(r1)')), post

/* -------- */
/* COLUMN 7 */
/* -------- */
eststo col7: spreg gs2sls $ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $legalor $colhist $contall if smpl_flag1 == 1, id(cid) dlmat(dobj1_global) elmat(dobj1_global) heteroskedastic

/* Obtain the effect of a 10th to 90th percentile move in diversity */
qui sum $ylist if e(sample) == 1 
local avgconf = r(mean) 
_pctile $xlist if e(sample) == 1, percentiles(10 90) 
eststo col7_me: nlcom (me_1090: _b[$xlist] * (1 + `avgconf') * (`r(r2)' - `r(r1)')), post

/* -------- */
/* COLUMN 8 */
/* -------- */
eststo col8: spivreg $ylist $geor $climate ($xlist = $ivlist) if smpl_flag1 == 1, id(cid) dlmat(dobj1_global) elmat(dobj1_global) heteroskedastic 

/* Obtain the effect of a 10th to 90th percentile move in diversity */
qui sum $ylist if e(sample) == 1 
local avgconf = r(mean) 
_pctile $xlist if e(sample) == 1, percentiles(10 90) 
eststo col8_me: nlcom (me_1090: _b[$xlist] * (1 + `avgconf') * (`r(r2)' - `r(r1)')), post

/* -------- */
/* COLUMN 9 */
/* -------- */
eststo col9: spivreg $ylist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $legalor $colhist $contall ($xlist = $ivlist) if smpl_flag2 == 1, id(cid) dlmat(dobj2_global) elmat(dobj2_global) heteroskedastic 

/* Obtain the effect of a 10th to 90th percentile move in diversity */
qui sum $ylist if e(sample) == 1 
local avgconf = r(mean) 
_pctile $xlist if e(sample) == 1, percentiles(10 90) 
eststo col9_me: nlcom (me_1090: _b[$xlist] * (1 + `avgconf') * (`r(r2)' - `r(r1)')), post 

/**************************************/
/* Print the results to a LaTeX table */
/**************************************/

estout col1 col2 col3 col4 col5 col6 col7 col8 col9 using "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\CROSS\Results\xtab_sarar.tex", style(tex) replace cells(b(star fmt(a2)) se(fmt(a2) par)) keep($xlist lambda:_cons rho:_cons, relax)order($xlist lambda:_cons rho:_cons) indicate("Continent dummies=${contall}"  "Controls for geography=${geor}" "Controls for Climate=${climate}" "Controls for ethnic diversity=${ethfrac} ${ethpolr}" "Controls for institutions=${inst} ${legalor}" "Controls for Colonial History=${colhist}" "Controls for oil, population, and income=${oilprod} ${popsize} ${ypercap}", labels("\$\times\$" " ")) stats(N, fmt(%9.0f) labels("Observations") layout(@)) varwidth(44) msign("\$-\$") label prehead("\begin{tabular*}{665pt}{@{\extracolsep{\fill}}lccccccccc}" "\toprule" "Cross-country sample:&\multicolumn{9}{c}{Global}\\" "\cmidrule(r){2-6}\cmidrule(lr){7-8}\cmidrule(l){9-10}") numbers mlabels("OLS" "OLS" "OLS" "OLS" "OLS" "OLS" "OLS" "2SLS" "2SLS") collabels(none) posthead("\midrule"  "&\multicolumn{9}{c}{Economic Freedom per year, 1970--2017}\\" "\cmidrule{2-10}") eqlabels("" "Spatial lag AR(1) of Economic Freedom ($\lambda$)" "Spatial lag AR(1) of error ($\rho$)", none) varlabels(, elist($xlist \addlinespace rho:_cons \addlinespace)) prefoot("\midrule") postfoot("\addlinespace")


estout col1_me col2_me col3_me col4_me col5_me col6_me col7_me col8_me col9_me using "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\CROSS\Results\xtab_sarar.tex", style(tex) append cells(b(nostar fmt(a2)) se(fmt(a2) par)) varwidth(44) msign("\$-\$") mlabels(none) collabels(none) varlabels(me_1090 "Effect of 10th--90th \%ile move in diversity") postfoot("\bottomrule\addlinespace" "\end{tabular*}")
  
  /***********************************************************/
/* Clean-up stored estimates, close the log file, and exit */
/***********************************************************/

est clear 
log close 

exit