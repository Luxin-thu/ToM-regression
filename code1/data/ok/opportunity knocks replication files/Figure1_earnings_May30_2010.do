*Tyler Williams
*5/30/2010
*This file uses datasets created by OK_gradesupdater_Feb5_2010.do
*It creates treatment/control kernel density plots by strata (gender, year, hs grade quartile)

*Set stata options
clear
set more off
set mem 200m
capture log close
cd "C:\Users\twill0k0\Downloads"

/* LOAD THE INDIVIDUAL LEVEL DATA */

use OKgradesUpdate_Feb5_2010, clear

/* PLOT GRADES AND EARNINGS DISTRIBUTIONS FOR TREATMENT AND CONTROL BY STRATA GROUPS */

*Loop over grades variables
foreach depvar in earnings {

	*Loop over time period
	local i = 1
	foreach length in 2008 {

		*Loop over strata
		local k = 1
		foreach strata in F_1 F_0 M_1 M_0 {
			qui ksmirnov `depvar'`length' if s_group=="`strata'", by(T)
			local p = round(r(p_cor),.001)
			local p = substr("`p'",1,4)
			if `p'==0 {
				local p = ".000"
			}
			kdensity `depvar'`length' if `depvar'`length'>=0 & `depvar'`length'<5000 & s_group=="`strata'" & T==1, kernel(gaussian) generate(Tearnings`strata' Tdens`strata') nograph
			kdensity `depvar'`length' if `depvar'`length'>=0 & `depvar'`length'<5000 & s_group=="`strata'" & T==0, kernel(gaussian) generate(Cearnings`strata' Cdens`strata') nograph
			disp "KS test p-value for `strata': 0`p'"
		}
	}
}

outsheet TdensF_1 - CearningsM_0 using nberfigure1points.csv, comma replace
