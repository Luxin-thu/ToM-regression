/********************************** 
/ Replication code  
/ "Social Trust in Polarized Times: 
   How Perceptions of Political Polarization 
   Affect Americans' Trust in Each Other" 
/ Amber Hye-Yon Lee 
/ March 3 2022 
**********************************/ 

use "PerPolarization_SocialTrust_Experiment.dta", clear


/* Figure 1 Manipulation Check */ 
mean perpolindx, over(condition)

reg perpolindx i.condition 
contrast {condition 1 0 -1}, effects //control vs less
contrast {condition 0 1 -1}, effects //control vs more
contrast {condition -1 1 0}, effects //less vs more  
margins i.condition, post
coefplot (., keep(2.condition) fc(gs2)) (., keep(1.condition) fc(gs8)) (., keep(3.condition) fc(white)), recast(bar) barwidth(0.5) lc(black) vert ylab(1(0.5)4, angle(0) format(%5.1fc) nogrid) leg(off) nooff ciop(recast(rcap) lc(black) lw(vthin)) citop plotr(m(sides)) xlab(1 `" "More" "Polarization" "' 2 `" "Less" "Polarization" "' 3 "Control") aspectratio(1.1) xtitle("Condition", m(top)) 


/* Figure 2 Social Trust Outcomes */ 

global covariates "i.female i.nonwhite i.nocollege i.divsep c.income c.qage" 

reg gstrustindx i.condition $covariates
reg ameritrustindx i.condition $covariates
reg strangertrindx i.condition $covariates

foreach var in gstrustindx ameritrustindx strangertrindx {
quietly reg `var' i.condition $covariates
eststo f2`var': margins i.condition, at((means) female nonwhite nocollege divsep income qage) post 
}
coefplot (f2gstrustindx, keep(2.condition) fc(gs2)) (f2gstrustindx, keep(1.condition) fc(gs8)) (f2gstrustindx, keep(3.condition) fc(white)), bylabel(Generalized Social Trust) || (f2ameritrustindx, keep(2.condition) fc(gs2)) (f2ameritrustindx, keep(1.condition) fc(gs8)) (f2ameritrustindx, keep(3.condition) fc(white)), bylabel(Trust in the American Citizenry) ||  (f2strangertrindx, keep(2.condition) fc(gs2)) (f2strangertrindx, keep(1.condition) fc(gs8)) (f2strangertrindx, keep(3.condition) fc(white)), bylabel(Willingess to Trust Strangers) ||, recast(bar) barwidth(0.5) lc(black) vert ylab(0.2(0.1)0.7, angle(0) format(%5.1fc) nogrid labs(small)) nooff ciop(recast(rcap) lc(black) lw(vthin)) citop plotr(m(top)) xlab(1 `" "Perceive" "More" "Polarization" "' 2 `" "Perceive" "Less" "Polarization" "' 3 "Control", labs(small) labgap(medium)) byopts(cols(3) leg(off)) subt(, s(medsmall) bc(gs12)) aspectratio(1.5)


/* Figure 3 Perceived Polarization x Donation */

reg charitypct i.condition i.donation $covariates
testparm i.condition
testparm i.donation

reg charitypct i.condition##i.donation $covariates
testparm i.condition#i.donation
margins i.condition#i.donation, at((means) female nonwhite nocollege divsep income qage) post
coefplot (., keep(2.condition#1.donation) fc(gs7) offset(-0.3)) (., keep(2.condition#2.donation) fc(white) offset(-0.4)) (., keep(1.condition#1.donation) fc(gs7) offset(0.1)) (., keep(1.condition#2.donation) fc(white)) (., keep(3.condition#1.donation) fc(gs7) offset(0.4)) (., keep(3.condition#2.donation) fc(white) offset(0.3)), recast(bar) barwidth(0.7) lc(black) vert ylab(0(10)60, angle(h) nogrid) ciopts(recast(rcap) lc(black) lw(vthin)) citop plotr(m(sides)) xlab(1.1 `" "Perceive" "More" "Polarization" "' 3.5 `""Perceive" "Less" "Polarization" "' 5.8 "Control") yti("Share of Bonus" "Donated to Charity (%)", orientation(h) yoffset(30)) leg(cols(1) pos(1) symxsize(3) yoffset(-13) order(1 3) lab(1 "Baseline") lab(3 "Donation Matching")) 


/* Figure 4 Americans Share Similar Values/Beliefs */

reg mostshareint i.condition

reg mostshareval i.condition
margins i.condition, post
coefplot (., keep(2.condition) fc(gs2)) (., keep(1.condition) fc(gs8)) (., keep(3.condition) fc(white)), recast(bar) barwidth(0.5) lc(black) vert ylabel(0 "0" 0.1 "10" 0.2 "20" 0.3 "30" 0.4 "40" 0.5 "50" 0.6 "60", angle(0) nogrid) leg(off) nooff ciop(recast(rcap) lc(black) lw(vthin)) citop plotr(m(sides)) xlab(1 `" "More" "Polarization" "' 2 `" "Less" "Polarization" "' 3 "Control") aspectratio(1.1) xtitle("Condition", m(top))




**************************
/* Coding of Variables	*/
**************************

* perceived polarization index (higher = perceive more polarization) 
alpha perpol_1-perpol_5
egen perpolindx=rowmean(perpol_1-perpol_5) 

* generalized social trust index (higher = more trusting) 
recode gstrust2 (1=2) (2=1), gen(gstrust2r)
recode trust2 (1=0) (2=1), gen(soctrust2)

egen gstrustindx=anycount(gstrust1 gstrust2r gstrust3), value(1)

* trust in the American people index (higher = more trusting) 
alpha ameritrust_*
egen ameritrustindx = rowmean(ameritrust_*)

* willingness to trust strangers index (higher = more trusting) 
alpha stranger_*
egen strangertrindx = rowmean(stranger_*)

foreach x in ameritrustindx gstrustindx strangertrindx {
	quietly summarize `x'
	replace `x'=(`x'-r(min))/(r(max)-r(min))
	}

* donation to charity 
egen tocharity = rowtotal(unitedway cancersociety goodwill redcross)
gen charitypct = tocharity / 2 

* most, almost all, all americans share R's values or beleifs 
recode sharevalues (1 2 3=1) (4/7=0), gen(mostshareval)
recode shareinterests (1 2 3=1) (4/7=0), gen(mostshareint)


* covariates 
recode qgender (2=1) (1=0), gen(female)
gen race=.
replace race=1 if qwhite==1 & qhispanic==0 //white, nonhispanic
replace race=2 if qblack==1 & qhispanic==0 //black, nonhispanic
replace race=3 if qhispanic==1 //hispanic 
replace race=4 if qhispanic==0 & qwhite!=1 & qblack!=1  //other 
recode race (1=0) (2 3 4 =1), gen(nonwhite) 
recode qincome (.=6), gen(income)
recode qeducation (1/5=1) (6/9=0), gen(nocollege)
recode marital (1/3=0) (4=1), gen(divsep)






