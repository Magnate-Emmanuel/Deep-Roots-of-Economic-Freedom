global ethfrac "avg_efrac " 
global ethpolr "avg_des_pol15" 
global geor " avg_abslat avg_distcr avg_ruggavg avg_suitavg avg_suitrng avg_elevavg avg_elevrng avg_island "
global colhist " avg_iscol_gbr avg_iscol_fra avg_iscol_oth" 
global legalor " avg_legor_uk avg_legor_fr"  
global popsize " avg_lnpop_wdi" 
global ypercap " avg_lngdppc_cu_usd_wdi "  
global contall " avg_africa avg_asia avg_namerica avg_samerica avg_oceania avg_europe " 
global oilprod " avg_anypetroleum_pet" 
global contold "avg_africa avg_asia" 
global climate "avg_ecofrac avg_ecopol avg_tmp_cru_v401 avg_pre_cru_v401 avg_tmp_vol_cru_v401 avg_pre_vol_cru_v401"


reg a_overall_index avg_pdiv_aa $timedum, vce(robust) noomitted vsquish
*outreg2 using econfree.rtf, replace keep(avg_pdiv_aa) ctitle(OLS) 

reg a_overall_index avg_pdiv_aa $climate $geor $timedum, vce(robust) noomitted vsquish
*outreg2 using econfree.rtf, append keep(avg_pdiv_aa $geor) ctitle(OLS) 

reg a_overall_index avg_pdiv_aa $climate $geor $ypercap $popsize $timedum, vce(robust) noomitted vsquish
*outreg2 using econfree.rtf, append keep(avg_pdiv_aa $geor $ypercap $popsize) ctitle(OLS) 

reg a_overall_index avg_pdiv_aa $climate $geor $ypercap $popsize $inst $timedum, vce(robust) noomitted vsquish
*outreg2 using econfree.rtf, append keep(avg_pdiv_aa $geor $ypercap $popsize $inst) ctitle(OLS) 

reg a_overall_index avg_pdiv_aa $climate $geor $ypercap $popsize $inst $legalor $timedum, vce(robust) noomitted vsquish
*outreg2 using econfree.rtf, append keep(avg_pdiv_aa $geor $ypercap $popsize $inst) ctitle(OLS) addtext(Legal Origin, YES)

reg a_overall_index avg_pdiv_aa $climate $geor $ypercap $popsize $inst $legalor $colhist $timedum, vce(robust) noomitted vsquish
*outreg2 using econfree.rtf, append keep(avg_pdiv_aa $geor $ypercap $popsize $inst) ctitle(OLS) addtext(Colonial History, YES, Legal Origin, YES)

reg a_overall_index avg_pdiv_aa $climate $geor $ypercap $popsize $inst $legalor $colhist $contall $timedum, vce(robust) noomitted vsquish
*outreg2 using econfree.rtf, append keep(avg_pdiv_aa $geor $ypercap $popsize $inst) ctitle(OLS) addtext(Colonial History, YES, Legal Origin, YES, Continental Dummies, YES)
