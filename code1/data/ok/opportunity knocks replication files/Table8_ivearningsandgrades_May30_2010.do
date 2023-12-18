*Tyler Williams
*5/30/2010
*This file uses datasets created by OK_gradesupdater_Feb5_2010.do
*It runs regressions to estimate treatment effects controlling for strata (gender, year, hs grade quartile) and covariates
*The regressions are IV with treatment assignment instrumenting for those that showed awareness of being in the program

*Set stata options
clear
set more off
set mem 200m
capture log close
cd "C:\Users\twill0k0\Downloads"

/* LOAD THE INDIVIDUAL LEVEL DATA */

use OKgradesUpdate_Feb5_2010, clear

/* SET THE STRATA CONTROLS LIST */

local stratacontrols ""
tab s_group_quart, gen(s_group_quart)
forvalues i=2(1)16 {
	local stratacontrols "`stratacontrols' s_group_quart`i'"
}

/* ADD IN ALL OTHER CONTROLS TO GET FULL CONTROLS LIST */

local fullcontrols "s_hsgrade3 s_mtongue_english s_mothergraddegree s_test1correct s_test2correct s_motherhsdegree s_mothercolldegree s_mothergraddegree s_mothereducmiss s_fatherhsdegree s_fathercolldegree s_fathergraddegree s_fathereducmiss `stratacontrols'" 

/* TABLE 8: ACADEMIC OUTCOME IV REGS FOR FULL YEAR, SPRING, AND FALL VARS WITH VARIOUS CONTROLS */

local tabnum=8

*Erase old tables
capture erase nbertable`tabnum'.csv
estimates clear

*Loop over subgroups
foreach restriction in "" "& s_test2correct==1" {

*FIRST STAGE REGRESSIONS
*Loop over subgroups: ^F is all women, ^M is all men, 1$ is all 1st years, 0$ is all 2nd years, . is everyone
foreach strata in F_1 F_0 ^F M_1 M_0 ^M 1$ 0$ . {
		
	foreach controlset in "`fullcontrols'" {
		qui eststo: reg anycontact T `controlset' if regexm(s_group,"`strata'") & earnings2008!=. `restriction', r
	}
}

*Output all effects in a table
qui esttab using nbertable`tabnum'.csv, cells(b(fmt(%9.3f)) se(fmt(%9.3f) star)) append nonumber keep(T) ///
varlabels(T "1st Stage Treatment Effect") mlabels(none) collabels(none) starlevels(* .1 ** .05 *** .01) stardetach
	
qui estimates clear

*Loop over academic variables
foreach depvar in earnings avggrade gpa earned ptsover70 {

	*SECOND STAGE REGRESSIONS
	*Loop over subgroups: ^F is all women, ^M is all men, 1$ is all 1st years, 0$ is all 2nd years, . is everyone
	foreach strata in F_1 F_0 ^F M_1 M_0 ^M 1$ 0$ . {
		
		foreach controlset in "`fullcontrols'" {
			qui eststo: ivregress 2sls `depvar'2008 `controlset' (anycontact = T) if regexm(s_group,"`strata'") `restriction', vce(robust)
		}
	}

	*Output all effects in a table
	qui esttab using nbertable`tabnum'.csv, cells(b(fmt(%9.3f)) se(fmt(%9.3f) star)) append nonumber keep(anycontact) ///
	varlabels(anycontact "2nd Stage") mlabels(none) collabels(none) starlevels(* .1 ** .05 *** .01) stardetach
	
	qui estimates clear

*Close loop over academic variables
}

*Close loop over subgroups
}

/* DISPLAY CONTROLS LISTS */

disp "`stratacontrols'"
disp "`fullcontrols'"
