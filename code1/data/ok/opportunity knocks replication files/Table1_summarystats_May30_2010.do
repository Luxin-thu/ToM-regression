*Tyler Williams
*5/30/2010
*This file uses datasets created by OK_gradesupdater_Feb5_2010.do
*It outputs summary statistics by treatment status, gender, and year

*Set stata options
clear
set more off
set mem 200m
capture log close
cd "C:\Users\twill0k0\Downloads"

/* LOAD THE INDIVIDUAL LEVEL GRADE AND ADMINISTRATIVE DATA */

use OKgradesUpdate_Feb5_2010, clear

/* CREATE A FEW VARIABLES */

gen C = 1-T
gen s_second_year = 1-s_first_year
gen s_female = 1-s_male

/* CHANGE COLLEGE GRAD/HIGH SCHOOL GRAD VARIABLES TO INCLUDE THOSE WITH HIGHER DEGREES */

replace s_motherhsdegree = 1 if s_mothercolldegree==1 | s_mothergraddegree==1
replace s_fatherhsdegree = 1 if s_fathercolldegree==1 | s_fathergraddegree==1
replace s_mothercolldegree = 1 if s_mothergraddegree==1
replace s_fathercolldegree = 1 if s_fathergraddegree==1

/* GENERATE CONTROLS HYPOTHETICAL EARNINGS VARIABLES */

gen controlswhoearned = gradeover702008 if T==0
gen controlsearnings = earnings2008 if T==0

/* SET THE STRATA CONTROLS LIST */

local stratacontrols ""
tab s_group_quart, gen(s_group_quart)
forvalues i=2(1)16 {
	local stratacontrols "`stratacontrols' s_group_quart`i'"
}

/* DESCRIPTIVE STATISTICS ON DEMOGRAPHIC VARIABLES BY STRATIFICATION GROUP AND TREATMENT */

*Make variable name labels
local var1 "Age"
local var2 "High school grade average"
local var3 "1st language is English"
local var4 "Mother finished college"
local var5 "Father finished college"
local var6 "Answered earnings test question correctly"
local var7 "Controls who would have been paid"
local var8 "Mean hypothetical earnings for controls"

*Erase the old table
capture erase nbertable1.csv
estimates clear

*TREATMENT DIFFERENCES CONTROLLING FOR STRATA
*Loop over variables to be summarized
local i = 1
foreach sumvar in s_age s_hsgrade3 s_mtongue_english s_mothercolldegree s_fathercolldegree s_test2correct {

	*Loop over subgroups: . is everyone
	foreach strata in F_1 F_0 M_1 M_0 . {
		qui eststo: reg `sumvar' T `stratacontrols' if regexm(s_group,"`strata'"), r
	}

	*Output all effects in a table
	qui esttab using nbertable1.csv, cells(b(fmt(3)) se(fmt(3) star)) append nonumber keep(T) noobs ///
	varlabels(T "`var`i''") mlabels(none) collabels(none) starlevels(* .1 ** .05 *** .01) stardetach
	
	qui estimates clear
	local ++i
}

*TREATMENT SAMPLE SIZES
*Loop over subgroups: . is everyone
foreach strata in F_1 F_0 M_1 M_0 . {

	qui eststo: reg s_age T `stratacontrols' if regexm(s_group,"`strata'"), r
	qui sum s_age if regexm(s_group,"`strata'") & T==1
	qui estadd scalar obs = r(N)
}

*Output all effects in a table
qui esttab using nbertable1.csv, cells(none) append nonumber ///
stats(obs, fmt(0) label("Observations")) mlabels(none) collabels(none) starlevels(* .1 ** .05 *** .01)

qui estimates clear

*CONTROL MEANS
*Loop over variables to be summarized
local i = 1
foreach sumvar in s_age s_hsgrade3 s_mtongue_english s_mothercolldegree s_fathercolldegree s_test2correct ///
controlswhoearned controlsearnings {

	*Loop over subgroups: . is everyone
	foreach strata in F_1 F_0 M_1 M_0 . {
		
		qui eststo: reg `sumvar' if regexm(s_group,"`strata'") & T==0, r
		qui sum `sumvar' if regexm(s_group,"`strata'") & T==0
		qui estadd r(mean)
		qui estadd r(sd)
	}

	*Output all effects in a table
	qui esttab using nbertable1.csv, cells(none) append nonumber ///
	stats(mean sd, fmt(3 3) labels("`var`i''" " ")) mlabels(none) collabels(none)
	
	qui estimates clear
	local ++i
}

*CONTROL SAMPLE SIZES
*Loop over subgroups: . is everyone
foreach strata in F_1 F_0 M_1 M_0 . {

	qui eststo: reg s_age T `stratacontrols' if regexm(s_group,"`strata'"), r
	qui sum s_age if regexm(s_group,"`strata'") & T==0
	qui estadd scalar obs = r(N)
}

*Output all effects in a table
qui esttab using nbertable1.csv, cells(none) append nonumber ///
stats(obs, fmt(0) label("Observations")) mlabels(none) collabels(none) starlevels(* .1 ** .05 *** .01)

qui estimates clear

*F TESTS OF JOINT SIGNIFICANCE (estimate each regression, save results, use "suest" to combine, then test)
*Loop over variables to be summarized
foreach sumvar in s_age s_hsgrade3 s_mtongue_english s_mothercolldegree s_fathercolldegree s_test2correct {

	*Loop over subgroups: . is everyone
	local stratanum = 1
	foreach strata in F_1 F_0 M_1 M_0 . {
		qui eststo `sumvar'`stratanum': reg `sumvar' T `stratacontrols' if regexm(s_group,"`strata'")
		local ++stratanum
	}
}

*Loop over subgroups
forvalues stratanum=1(1)5 {

	qui eststo model`stratanum': suest s_age`stratanum' s_hsgrade3`stratanum' s_mtongue_english`stratanum' s_mothercolldegree`stratanum' ///
	s_fathercolldegree`stratanum' s_test2correct`stratanum'
	qui test [s_age`stratanum'_mean]T [s_hsgrade3`stratanum'_mean]T [s_mtongue_english`stratanum'_mean]T ///
	[s_mothercolldegree`stratanum'_mean]T [s_fathercolldegree`stratanum'_mean]T [s_test2correct`stratanum'_mean]T
	qui local F = r(chi2)/r(df)
	qui estadd scalar F = `F'
	qui estadd r(p)
}

*Output all effects in a table
qui esttab model1 model2 model3 model4 model5 using nbertable1.csv, cells(none) append nonumber ///
stats(F p, fmt(3 3) labels("F test for joint significance" " ")) mlabels(none) collabels(none)
