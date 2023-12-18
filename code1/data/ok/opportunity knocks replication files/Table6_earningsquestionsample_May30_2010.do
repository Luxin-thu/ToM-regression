*Tyler Williams
*5/30/2010
*This file uses datasets created by OK_gradesupdater_Feb5_2010.do
*It runs regressions to estimate academic outcome treatment effects for various subgroups of the data, controlling for
*strata (gender, year, hs grade quartile) and covariates

*Set stata options
clear
set more off
set mem 200m
capture log close
cd "C:\Documents and Settings\Julie Gross\My Documents\RA Work\Pay for Grades\Experiment Analysis\Dta"

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

/* TABLE 6: ACADEMIC OUTCOME REGRESSIONS FOR FULL YEAR, SPRING, AND FALL VARS WITH VARIOUS
 CONTROLS, RESTRICTED TO SUBGROUPS */

local tabnum=6
*Erase old tables
capture erase nbertable`tabnum'.csv
estimates clear

*Loop over subgroup variables
foreach subgroupvar in s_test2correct {

*Loop over academic variables
foreach depvar in earnings avggrade gpa earned ptsover70 {

*RESULTS

*Loop over full year
foreach length in "" {

	*TREATMENT EFFECTS
	*Loop over subgroups: ^F is all women, ^M is all men, 1$ is all 1st years, 0$ is all 2nd years, . is everyone
	local labnum = 1
	foreach strata in F_1 F_0 ^F M_1 M_0 ^M 1$ 0$ . {

		foreach controlset in "`fullcontrols'" {
			qui eststo: reg `depvar'`length'2008 T `controlset' if regexm(s_group,"`strata'") & `subgroupvar'==1, r
		}
	}

	*Output all effects in a table
	qui esttab using nbertable`tabnum'.csv, cells(b(fmt(%9.3f)) se(fmt(%9.3f) star)) append nonumber keep(T) ///
	varlabels(T "Treatment Effect") mlabels(none) collabels(none) starlevels(* .1 ** .05 *** .01) stardetach stats(N, fmt(0))
	
	qui estimates clear
}

*Close loop over academic variables
}

*Close loop over subgroup variables
}

/* DISPLAY CONTROLS LISTS */

disp "`stratacontrols'"
disp "`fullcontrols'"
