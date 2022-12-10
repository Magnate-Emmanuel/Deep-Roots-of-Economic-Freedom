/*****************************************************************************/
/*** DESCRIPTION: REPLICATION CODE FOR "DIVERSITY AND ECONOMIC FREEDOM"            ***/
/*** AUTHORS:      ***/
/*****************************************************************************/

/****************************************************************************************************************************************************************************************************************/
/* TABLE: Population Diversity and the Incidence or Onset of Civil Conflict in Repeated Cross-Country Data - Robustness to Accounting for Spatiotemporal Dependence using Two-Way Clustering of Standard Errors */
/****************************************************************************************************************************************************************************************************************/

# delimit 
clear all 

set more off 

/*******************/
/* Open a log file */
/*******************/

capture log close 

/************************************/
/* Load the AAGK 5-yr panel dataset */
/************************************/


/*************************/
/* Generate time dummies */
/*************************/

xi, noomit i.period


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

global timedum "_Iyear_1970-_Iyear_2017" 
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

/******************************************************/
/* Generate the sample indicators for the regressions */
/******************************************************/

qui egen smpl_flag1 = rowmiss($ylist $xlist ${geor} $deepdet $ethineq $climate ${ethfrac} ${ethpolr}) 
qui replace smpl_flag1 = (smpl_flag1 == 0) 

qui egen smpl_flag2 = rowmiss($ylist $xlist ${geor} $deepdet $ethineq $climate ${ethfrac} ${ethpolr} ${colhist} ${legalor} ${inst} ${oilprod} ${popsize} ${ypercap}) 
qui replace smpl_flag2 = (smpl_flag2 == 0) 


/***********************/
/* Run the regressions */
/***********************/

/* -------- */
/* COLUMN 1 */
/* -------- */
eststo col1: cluster2 $ylist $xlist $timedum if smpl_flag1 == 1, fcluster(cid) tcluster(year)
*outreg2 using econfree.rtf, replace keep(avg_pdiv_aa) ctitle(OLS) 
estadd scalar adjr2 = e(r2_a) : col1 
 


/* -------- */
/* COLUMN 2 */
/* -------- */
eststo col2: cluster2 $ylist $xlist $geor $timedum if smpl_flag1 == 1, fcluster(cid) tcluster(year)
*outreg2 using econfree.rtf, append keep(avg_pdiv_aa $geor) ctitle(OLS) 
/* Obtain the adjusted R2 */
estadd scalar adjr2 = e(r2_a) : col2 



/* -------- */
/* COLUMN 3 */
/* -------- */
eststo col3: cluster2 $ylist $xlist $geor $climate $timedum if smpl_flag1 == 1, fcluster(cid) tcluster(year)
*outreg2 using econfree.rtf, append keep(avg_pdiv_aa $geor) ctitle(OLS) 
estadd scalar adjr2 = e(r2_a) : col3 

/* -------- */
/* COLUMN 4 */
/* -------- */
eststo col4: cluster2 $ylist $xlist $geor $climate $ethfrac $ethpolr $timedum if smpl_flag1 == 1, fcluster(cid) tcluster(year)
estadd scalar adjr2 = e(r2_a) : col4 

/* -------- */
/* COLUMN 5 */
/* -------- */
eststo col5: cluster2 $ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $timedum if smpl_flag1 == 1, fcluster(cid) tcluster(year)
estadd scalar adjr2 = e(r2_a) : col5

/* -------- */
/* COLUMN 6 */
/* -------- */
eststo col6: cluster2 $ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $timedum if smpl_flag1 == 1, fcluster(cid) tcluster(year)
estadd scalar adjr2 = e(r2_a) : col6 

/* -------- */
/* COLUMN 7 */
/* -------- */
eststo col7: cluster2 $ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $legalor $colhist $timedum if smpl_flag1 == 1, fcluster(cid) tcluster(year)
estadd scalar adjr2 = e(r2_a) : col7

/* -------- */
/* COLUMN 9 */
/* -------- */
eststo col9: cluster2 $ylist $xlist $geor $climate $ethfrac $ethpolr $ypercap $popsize $oilprod $inst $legalor $colhist $contall $timedum if smpl_flag1 == 1, fcluster(cid) tcluster(year)
estadd scalar adjr2 = e(r2_a) : col9


/**************************************/
/* Print the results to a LaTeX table */
/**************************************/

estout col1 col2 col3 col4 col5 col6 col7 col9 using "C:\Users\sarfo\Desktop\GRADUATE SCHOOL\SCHOOLS\HSE\COURSES\RESEARCH\Economic Freedom Predictors\Codes\Fraser Institute\PANEL\Results\mtab_panel_tclust.tex", style(tex) replace cells(b(star fmt(a2)) se(fmt(a2) par)) keep($xlist)   indicate(  "Continent dummies=${contall}" "Controls for geography=$geor"  "Controls for climate=$climate" "Controls of ethnic diversity=$ethfrac $ethpolr" "Controls for oil, population, and income=$ypercap $popsize $oilprod" "Controls for Institutions=$inst"  "Colonial History=$colhist" "Legal origin dummies=${legalor}", labels("\$\times\$" " "))   stats(N N_fclust adjr2, fmt(%9.0fc %9.0f a2) labels("Observations" "Countries" "Adjusted \$R^2\$") layout(@ @ @))   varwidth(40) msign("\$-\$") label   prehead("\midrule"  "&\multicolumn{12}{c}{Economic Freedom per year, 1970--2017}\\"  "\cmidrule{2-13}")   numbers mlabels("OLS" "OLS" "OLS" "OLS" "OLS" "OLS" "OLS" "OLS") collabels(none)   posthead("\midrule"  "&\multicolumn{12}{c}{Economic Freedom per year, 1970--2017}\\"  "\cmidrule{2-13}") varlabels(, elist(${inst} \addlinespace)) prefoot("\midrule") postfoot("\addlinespace" "\end{tabular*}")


  
  
/***********************************************************/
/* Clean-up stored estimates, close the log file, and exit */
/***********************************************************/

est clear 
log close 

exit 

