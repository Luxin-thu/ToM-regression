*Tyler Williams
*7/1/2010
*This file uses datasets created by OK_gradesupdater_Feb5_2010.do
*It runs OLS regressions to estimate treatment effects on academic outcomes controlling for strata (gender, year, hs grade 
*quartile) and covariates, and plots coefficients along with the control grades kernel density
*Uses command "renvars" - to install type "findit renvars" in the command line and click through the links

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

/* MAKE MATRIX TO HOLD THE THRESHOLDS, COEFFICIENTS, AND CONFIDENCE INTERVALS */

matrix thresholdeffects = J(84,6,.)

*PUT COEFFICIENTS AND CONFIDENCE INTERVALS IN A MATRIX

*Loop over male/female
local row = 1
foreach male in 1 0 {

	*Loop over first year or second year
	foreach firstyear in 1 0 {

		*Loop over grade thresholds
		forvalues grade=60(1)80 {

			*Put the thresholds in the matrix
			qui matrix thresholdeffects[`row',1] = `grade'

			*Mark whether the row is for males or females
			qui matrix thresholdeffects[`row',5] = `male'

			*Mark whether the row is for first years or second years
			qui matrix thresholdeffects[`row',6] = `firstyear'

			foreach controlset in "`fullcontrols'" {
				qui reg gradeover`grade'2008 T `fullcontrols' if s_male==`male' & s_first_year==`firstyear', r
				qui matrix thresholdeffects[`row',2] = _b[T]
				qui matrix thresholdeffects[`row',3] = _b[T]-1.96*_se[T]
				qui matrix thresholdeffects[`row',4] = _b[T]+1.96*_se[T]
			}
		local ++row
		}
	}
}

/* LOAD THE COURSE-LEVEL PANEL DATA */

use OKgradesUpdate_Jan18_2010_panel, clear

/* GET KDENSITY POINTS FOR COURSE GRADE DISTRIBUTION FOR CONTROL AND OUTSHEET ALL THE DATA POINTS */

*Add the treatment effects on grade thresholds to the dataset
qui svmat thresholdeffects, names(var)
renvars var1-var6 \ threshold betaY lowerY upperY effectmale effectfirstyear

*Loop over male/female
foreach male in 0 1 {

	*Loop over first year/second year
	foreach firstyear in 1 0 {

		kdensity grade if grade>=60 & grade<=80 & T==0 & s_male==`male' & s_first_year==`firstyear', kernel(gaussian) bwidth(2) gen(Tgrade`male'`firstyear' Tdens`male'`firstyear') nograph
	}
}

outsheet threshold - effectfirstyear Tdens01 - Tgrade10 using nberfigure2points.csv, comma replace
