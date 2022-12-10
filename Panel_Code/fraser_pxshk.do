/*****************************************************************************/
/*** DESCRIPTION: REPLICATION CODE FOR "DIVERSITY AND ECONOMIC FREEDOM"            ***/
/*** AUTHORS:     ***/
/*****************************************************************************/

/***********************************************************************************************************************************************************/
/* TABLE: Population Diversity and the Onset of Civil Conflict in Repeated Cross-Country Data - Robustness to Accounting for Commodity Export Price Shocks */
/***********************************************************************************************************************************************************/


# delimit 
clear all 

set more off 

/*******************/
/* Open a log file */
/*******************/

capture log close 

/************************************/
/* Load the AAGK 1-yr panel dataset */
/************************************/

/*************************/
/* Generate time dummies */
/*************************/

xi, noomit i.year

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
label variable iscol_gbr "Ever a U.K. colony dummy" 
label variable iscol_fra "Ever a French colony dummy" 
label variable iscol_oth "Ever a non-U.K./non-French colony dummy" 
label variable legor_uk    "British legal origin dummy" 
label variable legor_fr    "French legal origin dummy" 
  
/* Contemporary institutional variables */
label variable xconst_p4_v17 "Executive constraints, 1960--2017 average" 
label variable ddemoc_p4_v17  "Fraction of years under democracy, 1960--2017" 
label variable dautoc_p4_v17  "Fraction of years under autocracy, 1960--2017" 

/* Oil, population, and income variables */
label variable anypetroleum_pet                 "Oil or gas reserve discovery" 
label variable lnpop_wdi          "Log population, 1960--2017 average" 
label variable lngdppc_cu_usd_wdi "Log GDP per capita, 1960--2017 average" 

/****************************************/
/* Corruption Variables                 */
/****************************************/    
label variable v2exbribe "Executive bribery and corrupt exchanges"
label variable v2excrptps "Public sector corrupt exchanges"
label variable v2lgcrrpt "Legislature corrupt activities"
label variable v2jucorrdc "Judicial corruption decision"
label variable v2x_corr "Political corruption index"

/****************************************/
/*    Continent Dummies                 */
/****************************************/  
label variable africa "Africa Dummy"
label variable asia "Asia Dummy"
label variable namerica "North America Dummy" 
label variable samerica "South America Dummy"
label variable oceania "Oceania Dummy"
label variable europe "Europe Dummy"

/********************************************************/
/* Create global macro lists of the relevant covariates */
/********************************************************/

global timedum "_Iyear_1970-_Iyear_2007" 
global ethfrac "efrac " 
global ethpolr "des_pol15" 
global ylist "overall_index"
global xlist "pdiv_aa"
global pd_vars "pdiv_aa_w pdiv_aa_b"
global ivlist "mdist_addis"
global geor " abslat distcr ruggavg suitavg suitrng elevavg elevrng island "
global colhist " iscol_gbr iscol_fra iscol_oth" 
global legalor " legor_uk legor_fr"  
global popsize " lnpop_wdi" 
global ypercap " lngdppc_cu_usd_wdi "  
global contall " africa asia namerica samerica oceania" 
global contold "africa asia"
global inst " xconst_p4_v17 ddemoc_p4_v17 dautoc_p4_v17"
global oilprod " anypetroleum_pet" 
global climate "ecofrac ecopol tmp_cru_v401 pre_cru_v401 tmp_vol_cru_v401 pre_vol_cru_v401"
global deepdet "ln_yst_aa ln_statehist ln_origtime ln_frontdist1500" 
global ethineq "grg_a_lum00pc v1_a_lum00pc"
global ethalt1 "ef_fearon des_pol15" 
global ethalt2 "lfrac des_pol15" 
global ethalt3 "rfrac des_pol15" 
global ethalt4 "frac_fear_emr er_fear_delta005_emr gini_fear_delta005_PERCAPTA_emr" 
global moregeo "lmtnest_fl03 ncontig_fl03 disease"
global shocks1 "pshock_npi_p pshock_npi_p1 pshock_npi_p2" 
global shocks2 "pshock_annual_npi_p pshock_annual_npi_p1 pshock_annual_npi_p2 pshock_perennial_npi_p pshock_perennial_npi_p1 pshock_perennial_npi_p2 pshock_extract_npi_p pshock_extract_npi_p1 pshock_extract_npi_p2"

/******************************************************/
/* Generate the sample indicators for the regressions */
/******************************************************/

qui egen smpl_flag1 = rowmiss($ylist $xlist ${geor} $shocks1 $shocks2 $ethineq $climate ${ethfrac} ${ethpolr}) 
qui replace smpl_flag1 = (smpl_flag1 == 0) 

qui egen smpl_flag2 = rowmiss($ylist $xlist ${geor} $deepdet $ethineq $climate $shocks1 $shocks2 ${ethfrac} ${ethpolr} ${colhist} ${legalor} ${inst}) 
qui replace smpl_flag2 = (smpl_flag2 == 0) 

/* -------- */
/* COLUMN 1 */
/* -------- */
eststo col1: reg $ylist $xlist $shocks1 $geor $climate $timedum if smpl_flag1 == 1, vce(cluster cid) noomitted vsquish
*outreg2 using econfree.rtf, append keep(avg_pdiv_aa $geor) ctitle(OLS) 

/* Obtain the adjusted R2 */
estadd scalar adjr2 = e(r2_a) : col1 

/* Obtain the effect of a 10th to 90th percentile move in diversity */
qui sum $ylist if e(sample) == 1 
local avgconf = r(mean) 
_pctile $xlist if e(sample) == 1, percentiles(10 90) 
eststo col1_me: nlcom (me_1090: _b[$xlist] * (1 + `avgconf') * (`r(r2)' - `r(r1)')), post 

/* Obtain the partial R2 of diversity */
pcorr $ylist $xlist $shocks1 $geor $climate if e(sample) == 1 
matrix pc_mtrx = r(p_corr) 
estadd scalar pc_pdiv = pc_mtrx[rownumb(pc_mtrx, "$xlist"), 1]^2 : col1 


/* -------- */
/* COLUMN 2 */
/* -------- */
eststo col2: reg $ylist $xlist $shocks1 $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $legalor $colhist $contall $timedum if smpl_flag1 == 1, vce(cluster cid) noomitted vsquish

/* Obtain the adjusted R2 */
estadd scalar adjr2 = e(r2_a) : col2 

/* Obtain the effect of a 10th to 90th percentile move in diversity */
qui sum $ylist if e(sample) == 1 
local avgconf = r(mean) 
_pctile $xlist if e(sample) == 1, percentiles(10 90) 
eststo col2_me: nlcom (me_1090: _b[$xlist] * (1 + `avgconf') * (`r(r2)' - `r(r1)')), post 

/* Obtain the partial R2 of diversity */
pcorr $ylist $xlist $geor $shocks1 $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $legalor $colhist $contall if e(sample) == 1 
matrix pc_mtrx = r(p_corr) 
estadd scalar pc_pdiv = pc_mtrx[rownumb(pc_mtrx, "$xlist"), 1]^2 : col2


/* -------- */
/* COLUMN 3 */
/* -------- */
eststo col3: reg $ylist $xlist $shocks2 $geor $climate $timedum if smpl_flag1 == 1, vce(cluster cid) noomitted vsquish
*outreg2 using econfree.rtf, append keep(avg_pdiv_aa $geor) ctitle(OLS) 

/* Obtain the adjusted R2 */
estadd scalar adjr2 = e(r2_a) : col3 

/* Obtain the effect of a 10th to 90th percentile move in diversity */
qui sum $ylist if e(sample) == 1 
local avgconf = r(mean) 
_pctile $xlist if e(sample) == 1, percentiles(10 90) 
eststo col3_me: nlcom (me_1090: _b[$xlist] * (1 + `avgconf') * (`r(r2)' - `r(r1)')), post 

/* Obtain the partial R2 of diversity */
pcorr $ylist $xlist $shocks2 $geor $climate if e(sample) == 1 
matrix pc_mtrx = r(p_corr) 
estadd scalar pc_pdiv = pc_mtrx[rownumb(pc_mtrx, "$xlist"), 1]^2 : col3 


/* -------- */
/* COLUMN 4*/
/* -------- */
eststo col4: reg $ylist $xlist $shocks2 $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $legalor $colhist $contall $timedum if smpl_flag1 == 1, vce(cluster cid) noomitted vsquish

/* Obtain the adjusted R2 */
estadd scalar adjr2 = e(r2_a) : col4 

/* Obtain the effect of a 10th to 90th percentile move in diversity */
qui sum $ylist if e(sample) == 1 
local avgconf = r(mean) 
_pctile $xlist if e(sample) == 1, percentiles(10 90) 
eststo col4_me: nlcom (me_1090: _b[$xlist] * (1 + `avgconf') * (`r(r2)' - `r(r1)')), post 

/* Obtain the partial R2 of diversity */
pcorr $ylist $xlist $shocks2 $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $legalor $colhist $contall if e(sample) == 1 
matrix pc_mtrx = r(p_corr) 
estadd scalar pc_pdiv = pc_mtrx[rownumb(pc_mtrx, "$xlist"), 1]^2 : col4



/* --------- */
/* COLUMN 5 */
/* --------- */
eststo col5: ivregress 2sls $ylist $geor $shocks1 $climate $timedum ($xlist = $ivlist) if smpl_flag1 == 1, vce(cluster cid) first noomitted vsquish 

/* Store the first-stage statistics for subsequent use */
estat firststage 
matrix fs_mtrx = r(singleresults) 

/* Obtain the effect of a 10th to 90th percentile move in diversity */
qui sum $ylist if e(sample) == 1 
local avgconf = r(mean) 
_pctile $xlist if e(sample) == 1, percentiles(10 90) 
eststo col5_me: nlcom (me_1090: _b[$xlist] * (1 + `avgconf') * (`r(r2)' - `r(r1)')), post 

/* Obtain the first-stage regression coefficients and the first-stage F-statistic */
eststo col5_fs: reg $xlist $ivlist $shocks1 $geor $climate if smpl_flag1 == 1, vce(cluster cid) noomitted vsquish
estadd scalar fs_fsta = fs_mtrx[1,4] : col5_fs 


/* --------- */
/* COLUMN 6 */
/* --------- */
eststo col6: ivregress 2sls $ylist $shocks1 $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $legalor $colhist $contall $timedum ($xlist = $ivlist) if smpl_flag2 == 1, vce(cluster cid) first noomitted vsquish

/* Store the first-stage statistics for subsequent use */
estat firststage 
matrix fs_mtrx = r(singleresults) 

/* Obtain the effect of a 10th to 90th percentile move in diversity */
qui sum $ylist if e(sample) == 1 
local avgconf = r(mean) 
_pctile $xlist if e(sample) == 1, percentiles(10 90) 
eststo col6_me: nlcom (me_1090: _b[$xlist] * (1 + `avgconf') * (`r(r2)' - `r(r1)')), post 

/* Obtain the first-stage regression coefficients and the first-stage F-statistic */
eststo col6_fs: reg $xlist $ivlist $shocks1 $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $legalor $colhist $contall if smpl_flag2 == 1, vce(cluster cid) noomitted vsquish
estadd scalar fs_fsta = fs_mtrx[1,4] : col6_fs 

/* --------- */
/* COLUMN 7 */
/* --------- */
eststo col7: ivregress 2sls $ylist $shocks2 $geor $climate $timedum ($xlist = $ivlist) if smpl_flag1 == 1, vce(cluster cid) first noomitted vsquish 

/* Store the first-stage statistics for subsequent use */
estat firststage 
matrix fs_mtrx = r(singleresults) 

/* Obtain the effect of a 10th to 90th percentile move in diversity */
qui sum $ylist if e(sample) == 1 
local avgconf = r(mean) 
_pctile $xlist if e(sample) == 1, percentiles(10 90) 
eststo col7_me: nlcom (me_1090: _b[$xlist] * (1 + `avgconf') * (`r(r2)' - `r(r1)')), post 

/* Obtain the first-stage regression coefficients and the first-stage F-statistic */
eststo col7_fs: reg $xlist $ivlist $shocks2 $geor $climate if smpl_flag1 == 1, vce(cluster cid) noomitted vsquish
estadd scalar fs_fsta = fs_mtrx[1,4] : col7_fs 


/* --------- */
/* COLUMN 8 */
/* --------- */
eststo col8: ivregress 2sls $ylist $shocks2 $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $legalor $colhist $contall $timedum ($xlist = $ivlist) if smpl_flag2 == 1, vce(cluster cid) first noomitted vsquish

/* Store the first-stage statistics for subsequent use */
estat firststage 
matrix fs_mtrx = r(singleresults) 

/* Obtain the effect of a 10th to 90th percentile move in diversity */
qui sum $ylist if e(sample) == 1 
local avgconf = r(mean) 
_pctile $xlist if e(sample) == 1, percentiles(10 90) 
eststo col8_me: nlcom (me_1090: _b[$xlist] * (1 + `avgconf') * (`r(r2)' - `r(r1)')), post 

/* Obtain the first-stage regression coefficients and the first-stage F-statistic */
eststo col8_fs: reg $xlist $ivlist $shocks2 $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $legalor $colhist $contall if smpl_flag2 == 1, vce(cluster cid) noomitted vsquish
estadd scalar fs_fsta = fs_mtrx[1,4] : col8_fs 


/**************************************/
/* Print the results to a LaTeX table */
/**************************************/

estout col1 col2 col3 col4 col5 col6 col7 col8 using "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\PANEL\Results\mtab_panel_pxshk.tex", style(tex) replace cells(b(star fmt(a2)) se(fmt(a2) par)) keep($xlist $shocks1 $shocks2, relax)  indicate( "Continent dummies=${contall}" "Controls for geography=$geor"  "Controls of ethnic diversity=$ethfrac $ethpolr" "Controls for oil, population, and income=$ypercap $popsize $oilprod" "Controls for Institutions=$inst"  "Colonial History=$colhist" "Legal origin dummies=${legalor}" , labels("\$\times\$" " ")) stats(N pc_pdiv  adjr2, fmt(%9.0f  a2 a2) labels("Observations"  "Partial \$R^2\$ of population diversity"   "Adjusted \$R^2\$") layout(@ @ @ @ @)) varwidth(50) msign("\$-\$") label prehead("\begin{tabular*}{895pt}{@{\extracolsep{\fill}}lcccccccccc}" "\toprule" "Cross-country sample:&\multicolumn{7}{c}{Global}&\multicolumn{1}{c}{Old World}&\multicolumn{2}{c}{Global}\\" "\cmidrule(r){2-7}\cmidrule(lr){9-10}\cmidrule(l){10-11}") numbers mlabels("OLS" "OLS" "OLS" "OLS" "OLS" "OLS" "OLS" "OLS" "2SLS" "2SLS") collabels(none) posthead("\midrule"  "&\multicolumn{12}{c}{Economic Freedom per year, 1970--2007}\\"  "\cmidrule{2-13}") varlabels(, elist(${inst} \addlinespace)) prefoot("\midrule") postfoot("\addlinespace") 

estout col1_me col2_me col3_me col4_me col5_me col6_me col7_me col8_me  using "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\PANEL\Results\mtab_panel_pxshk.tex", style(tex) append cells(b(nostar fmt(a2)) se(fmt(a2) par)) varwidth(50) msign("\$-\$") extracols(8)  mlabels(none) collabels(none) varlabels(me_1090 "Effect of 10th--90th \%ile move in diversity") 

estout col5_fs col6_fs col7_fs col8_fs using "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\PANEL\Results\mtab_panel_pxshk.tex", style(tex) append cells(b(nostar fmt(a2)) se(fmt(a2) par)) keep(mdist_addis, relax) stats(fs_fsta, fmt(a2) labels("First-stage \$F\$ statistic") layout(@)) varwidth(50) msign("\$-\$") extracols(-9/0) label mlabels(none) collabels(none) posthead("FIRST STAGE&\multicolumn{10}{c}{}&\multicolumn{2}{c}{Population diversity}\\"  "&\multicolumn{10}{c}{}&\multicolumn{2}{c}{(ancestry adjusted)}\\"  "\cmidrule{12-13}") varlabels(,elist(mdist_addis \addlinespace)) postfoot("\bottomrule\addlinespace" "\end{tabular*}") 