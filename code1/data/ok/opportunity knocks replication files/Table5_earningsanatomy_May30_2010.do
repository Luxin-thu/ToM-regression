*Tyler Williams
*5/30/2010
*This file uses datasets created by OK_gradesupdater_Feb5_2010.do
*It runs OLS regressions to estimate treatment effects on academic outcomes controlling for strata (gender, year, hs grade 
*quartile) and covariates

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

/* TABLE 5: CONTROL MEANS AND ACADEMIC OUTCOME REGRESSIONS FOR FULL YEAR, SPRING, AND FALL VARS WITH ALL CONTROLS */

*Loop over academic variables
local tabnum = 5
*Erase old tables
capture erase nbertable`tabnum'.csv
estimates clear
foreach depvar in earned ptsover70 {

*RESULTS

*Loop over full year
foreach length in "" {

	*CONTROL MEANS
	*Loop over subgroups: ^F is all women, ^M is all men, 1$ is all 1st years, 0$ is all 2nd years, . is everyone
	foreach strata in F_1 F_0 ^F M_1 M_0 ^M 1$ 0$ . {

		qui eststo: reg `depvar'`length'2008 if regexm(s_group,"`strata'"), r
		qui sum `depvar'`length'2008 if regexm(s_group,"`strata'") & T==0
		qui estadd r(mean)
		qui estadd r(sd)
	}

	*Output all effects in a table
	qui esttab using nbertable`tabnum'.csv, cells(none) append nonumber ///
	stats(mean sd, fmt(%9.3f %9.3f) labels("Control Mean" "SD")) mlabels(none) collabels(none)
	
	qui estimates clear

	*TREATMENT EFFECTS
	*Loop over subgroups: ^F is all women, ^M is all men, 1$ is all 1st years, 0$ is all 2nd years, . is everyone
	local labnum = 1
	foreach strata in F_1 F_0 ^F M_1 M_0 ^M 1$ 0$ . {

		foreach controlset in "`fullcontrols'" {
			qui eststo: reg `depvar'`length'2008 T `controlset' if regexm(s_group,"`strata'"), r
		}
	}

	*Output all effects in a table
	qui esttab using nbertable`tabnum'.csv, cells(b(fmt(%9.3f)) se(fmt(%9.3f) star)) append nonumber keep(T) ///
	varlabels(T "Treatment Effect") mlabels(none) collabels(none) starlevels(* .1 ** .05 *** .01) stardetach stats(N, fmt(0))
	
	qui estimates clear
}

*Close loop over academic variables
}

/* DISPLAY CONTROLS LISTS */

disp "`stratacontrols'"
disp "`fullcontrols'"
