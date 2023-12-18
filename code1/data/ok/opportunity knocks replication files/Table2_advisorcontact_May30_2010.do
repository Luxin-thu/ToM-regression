*Tyler Williams
*5/30/2010
*This file uses the dataset created by OK_gradesupdater_Feb5_2010.do
*It outputs summary statistics on advisor contacts by gender and year

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

/* DESCRIPTIVE STATISTICS ON ADVISOR CONTACT VARIABLES BY STRATIFICATION GROUP */

*Make variable name labels
local var1 "At least one email to advisor (fall)"
local var2 "At least one email to advisor (spring)"
local var3 "At least one email to advisor (fall or spring)"
local var4 "Checked scholarship earnings online"
local var5 "Other email to program website"
local var6 "Any contact"

*Erase the old table
capture erase nbertable2.csv
estimates clear

*MEAN AND SD
*Loop over variables to be summarized
local i = 1
foreach sumvar in springadvisorcontact falladvisorcontact advisorcontact checkedscholarship ///
websitecontact anycontact {

	*Loop over subgroups: . is everyone
	foreach strata in F_1 F_0 ^F M_1 M_0 ^M 1$ 0$ . {
		qui eststo: reg `sumvar' if regexm(s_group,"`strata'") & T==1, r
		qui sum `sumvar' if regexm(s_group,"`strata'") & T==1
		qui estadd r(mean)
		qui estadd r(sd)
	}

	*Output all effects in a table
	qui esttab using nbertable2.csv, cells(none) append nonumber mlabels(none) collabels(none) noobs ///
	stats(mean sd, fmt(%9.3f %9.3f) labels("`var`i''" " "))
	
	qui estimates clear
	local ++i
}

*ADD OBSERVATIONS TO THE BOTTOM OF THE TABLE

*Add sample sizes to the tables
qui estimates clear

*Loop over subgroups: . is everyone
foreach strata in F_1 F_0 ^F M_1 M_0 ^M 1$ 0$ . {
	qui eststo: reg anycontact if regexm(s_group,"`strata'") & T==1, r
}

qui esttab using nbertable2.csv, cells(none) append nonumber ///
stats(N, fmt(0) labels("Observations")) mlabels(none) collabels(none)
