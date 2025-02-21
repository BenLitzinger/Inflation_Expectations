version 10.1
drop _all
set memory 512m
set more off
set matsize 800
set type double, permanently


//local dr "C:\Users\benli\Documents\Econometrics Chair\Michigan_Inflation_Rates"
//local outdir "Users\benli\OneDrive - WHU\Econometrics Chair\Michigan_Inflation_ Rates\output"
//cd "`dr'/Dropbox/temp/Inflexp"

global qtr = 1 

capture program drop ipolmedian
program define ipolmedian
  version 8.2
  syntax varlist [, by(varlist)] 
  tempvar med m ntot cf cl ch f l u medval loval hival
  sort `by'
  by `by': egen `m' = mean(`varlist')
  by `by': egen `med' = median(`varlist')
  by `by': egen `ntot' = count(`varlist')   // count only uses non-missing values
  quietly gen `medval' = `varlist' if `varlist' == `med'
  quietly gen `loval' = `varlist' if `varlist' < `med'
  quietly gen `hival' = `varlist' if `varlist' > `med'
  by `by': egen `cl' = count(`loval')
  by `by': egen `ch' = count(`hival')
  by `by': egen `cf' = count(`medval')
  by `by': egen `l' = max(`loval')
  by `by': egen `u' = min(`hival')
  local medname "`varlist'_med" 
  quietly gen `medname' = .
  //quietly replace `medname' = `med' + 0.5*(`ch'/`cf')*(`u'-`med') + 0.5*(`cl'/`cf')*(`l'-`med') if (`cf' > 0 & `ntot' > 1)   // collapses to simple formula for uniform grid
  //quietly replace `medname' = `l' + ((`ntot'/2-`cl')/`cf')*(`u'-`l') if (`cf' > 0 & `ntot' > 1)   // collapses to simple formula for uniform grid
  quietly replace `medname' = `med' + max((`ntot'/2-`cl'-`cf'/2)/(0.5*`cf'),0)*0.5*(`u'-`med') + min((`ntot'/2-`cl'-`cf'/2)/(0.5*`cf'),0)*0.5*(`med'-`l')   if (`cf' > 0 & `ntot' > 1)
  quietly replace `medname' = `med' if (`cf' == 0 | `ch' < 2 | `cl' < 2 | `ntot' == 1 ) 
end


capture program drop wtipolmedian
program define wtipolmedian
  version 8.2
  syntax varlist [, by(varlist)] 
  tempvar med m ntot cf cl ch f l u medval wtx wts lowts hiwts medwts runsum grmed loval hival
  sort `by' `varlist'
  qui gen `wtx' = wt*`varlist'
  qui gen `wts' = wt if `wtx' ~= .
  by `by': egen `ntot' = total(`wts')
  by `by': egen `m' = total(`wtx')
  qui replace `m' = `m'/`ntot'
  
  by `by': gen `runsum' = sum(`wts')
  qui replace `runsum' = `runsum'/`ntot'
  qui gen `grmed' = 999
  qui replace `grmed' = `varlist' if `runsum' >= 0.5
  by `by': egen `med' = min(`grmed') 
  
  quietly gen `medwts' = `wts' if `varlist' == `med'
  quietly gen `lowts' = `wts' if `varlist' < `med'
  quietly gen `hiwts' = `wts' if (`varlist' > `med' & `varlist' < .)
  by `by': egen `cl' = total(`lowts') 
  by `by': egen `ch' = total(`hiwts') 
  by `by': egen `cf' = total(`medwts') 
  quietly gen `loval' = `varlist' if `varlist' < `med'
  quietly gen `hival' = `varlist' if (`varlist' > `med' & `varlist' < .)
  by `by': egen `l' = max(`loval') 
  by `by': egen `u' = min(`hival') 
  //qui replace chk1 = `cl' 
  //qui replace chk2 = `ch'
  //qui replace chk3 = `cf' 
  //qui replace chk4 = `med'
  //qui replace chk5 = `ntot' 
 
  local medname "`varlist'_med" 
  quietly gen `medname' = .
  //quietly replace `medname' = `l' + ((`ntot'/2-`cl')/`cf')*(`u'-`l') if (`cf' > 0 & `ntot' > 1)   // collapses to simple formula for uniform grid
  quietly replace `medname' = `med' + max((`ntot'/2-`cl'-`cf'/2)/(0.5*`cf'),0)*0.5*(`u'-`med') + min((`ntot'/2-`cl'-`cf'/2)/(0.5*`cf'),0)*0.5*(`med'-`l')   if (`cf' > 0 & `ntot' > 1)
  quietly replace `medname' = `med' if (`cf' == 0 | `ch' < 2 | `cl' < 2 | `ntot' == 1 )  // cf = 0 can happen when Stata sets median with tie break to ~.5
end


// simple regression with two-way clustered s.e. 

capture program drop clust2reg
program clust2reg, eclass
   version 10
   syntax varlist [if] [, noconstant]
   marksample touse
   gettoken depvar explvar: varlist
   tempname Vc Vt V eV eb
   reg `depvar' `explvar' if `touse', `constant' vce(cluster cohort) 
   matrix `Vc' = e(V)
   matrix `eb' = e(b)
   qui reg `depvar' `explvar' if `touse', `constant' vce(cluster yyyyq) 
   matrix `Vt' = e(V)
   qui reg `depvar' `explvar' if `touse', `constant' vce(robust) 
   matrix `V' = e(V)
   matrix `eV' = `Vc' + `Vt' - `V'
   
   scalar e_N=e(N)
   scalar e_df_m = e(df_m)
   scalar e_df_r = e(df_r)
   scalar e_F = e(F)
   scalar e_r2 = e(r2)
   scalar e_rmse = e(rmse)
   scalar e_mss = e(mss)
   scalar e_rss = e(rss)
   scalar e_r2_a = e(r2_a)
   scalar e_ll = e(ll)
   scalar e_ll_0 = e(ll_0)
   
   ereturn post `eb' `eV'
   //ereturn clear
   ereturn scalar N = e_N
   //ereturn scalar df_m = e_df_m
   //ereturn scalar df_r = e_df_r
   //ereturn scalar F= e_F
   //ereturn scalar r2= e_r2
   ereturn scalar rmse = e_rmse
   //ereturn scalar mss = e_mss
   //ereturn scalar rss = e_rss
   ereturn scalar r2_a = e_r2_a
   //ereturn scalar ll = e_ll
   //ereturn scalar ll_0 = e_ll_0
   
end   
    

	
	

// --------------------------------------------------------------------------- //
//																			   //
// Form cohorts and impute missing percentage expectations                     //
//   from categorical responses                                                //
//                                                                             //
// --------------------------------------------------------------------------- //


// UNCOMMENT THIS SECTION TO COMPUTE COHORT LEVEL DATA FROM MICHIGAN SURVEY 
// MICRO DATA




/******************* load michigan survey data****************************************/


//insheet using "michsurvey_stata.txt", clear
import delimited using "michsurvey_nec.csv", clear varnames(1)
foreach var in px1q1 px1q2 px1 px5q1 px5q2 px5 {
    destring `var', replace
}

******Inflation Percentage******
gen infl1pct = px1q2 / 100
gen infl1cat = .
replace infl1cat = 0 if px1q1 == 5
replace infl1cat = 1 if px1q1 == 3
replace infl1cat = 2 if (px1q1 == 1 | px1q1 == 2)

gen infl5pct = px5q2/100
gen infl5cat = .
replace infl5cat = 0 if px5q1 == 5
replace infl5cat = 1 if px5q1 == 3
replace infl5cat = 2 if (px5q1 == 1 | px5q1 == 2)

gen infl1cat_0 = 0
replace infl1cat_0 =1 if infl1cat == 0 
replace infl1cat_0 =. if infl1cat == .
gen infl1cat_1 = 0
replace infl1cat_1 =1 if infl1cat == 1 
replace infl1cat_1 =. if infl1cat == .
gen infl1cat_2 = 0
replace infl1cat_2 =1 if infl1cat == 2 
replace infl1cat_2 =. if infl1cat == .

gen infl5cat_0 = 0
replace infl5cat_0 =1 if infl5cat == 0 
replace infl5cat_0 =. if infl5cat == .
gen infl5cat_1 = 0
replace infl5cat_1 =1 if infl5cat == 1 
replace infl5cat_1 =. if infl5cat == .
gen infl5cat_2 = 0
replace infl5cat_2 =1 if infl5cat == 2 
replace infl5cat_2 =. if infl5cat == .


********************* sample requirements *****************************
destring age, replace
keep if (age > 24 & age < 75) //nur Personen im erwerbsfähigen Alter
keep if yyyy < 2024 //allen neuen Daten mitnehmen


gen mm = yyyymm-yyyy*100 //Monat wird berechnet
gen q = 1 //Quartal wird eingesetzt
replace q = 2 if (mm > 3 & mm < 7)   
replace q = 3 if (mm > 6 & mm < 10)
replace q = 4 if mm > 9
//gen yyyyq = yyyy*10+q
gen cohort = yyyy - age //Gleich wie var birthy
replace wt = 1 if wt == .   // no missing values in time periods in which sample weights are available

save basefile, replace



// Note that finfl and finfl5 defined starting from the beginning of the month 
// following the survey month to 12 or 60 months out. 
// By forming cohort data by quarter, we are averaging these monthly 
// numbers within a quarter

	
//overall median
wtipolmedian infl1pct, by(yyyymm)
wtipolmedian infl5pct, by(yyyymm) 
gen fullmedinfl1pct = infl1pct_med
gen fullmedinfl5pct = infl5pct_med
drop infl1pct_med
drop infl5pct_med


collapse (mean) fullmedinfl1pct fullmedinfl5pct infl1pct infl5pct  ///
             infl1cat_0 infl5cat_0 infl1cat_1 infl5cat_1 infl1cat_2 infl5cat_2 /// 
            (count) n5 = infl5pct n1 = infl1pct [aweight=wt], by(yyyymm yyyyq yyyy q age)


gen cohort = yyyy-age

xi: clust2reg infl5pct i.yyyymm infl5cat_0 infl5cat_2 if (infl5pct~=. & yyyy > 1981)  
xi: reg infl5pct i.yyyymm infl5cat_0 infl5cat_2  if (infl5pct~=. & yyyy > 1981)   
predict infl5pctimp

xi: clust2reg infl1pct i.yyyymm infl1cat_0 infl1cat_2  if (infl1pct~=.) // & yyyy > 1981)
xi: reg infl1pct i.yyyymm infl1cat_0 infl1cat_2  if (infl1pct~=.) // & yyyy > 1981) 
predict infl1pctimp


collapse (mean) yyyyq cohort fullmedinfl1pct fullmedinfl5pct infl1pct infl5pct ///
             infl1cat_0 infl5cat_0 infl1cat_1 infl5cat_1 infl1cat_2 infl5cat_2 infl1pctimp infl5pctimp /// 
            (sum) n5 = n5 n1 = n1, by(yyyy q age)
 
egen chkinfl5 = mean(infl5pct), by(yyyy q)    // only impute when entire year missing. Sometimes 1 to 3 age groups missing. Don't impute those. 
egen chkinfl1 = mean(infl1pct), by(yyyy q)    

gen imp5flag = 0
replace imp5flag = 1 if (chkinfl5 == . & infl5pctimp ~= .)
replace infl5pct = infl5pctimp if chkinfl5 ==.

gen imp1flag = 0
replace imp1flag = 1 if (chkinfl1 == . & infl1pctimp ~= .)
replace infl1pct = infl1pctimp if chkinfl1  ==.

//gen err5 = infl5pct-finfl5
//gen err1 = infl1pct-finfl

ipolmedian infl1pctimp, by(yyyy q)
ipolmedian infl5pctimp, by(yyyy q) 
gen medinfl1pctimp = infl1pctimp_med
gen medinfl5pctimp = infl5pctimp_med
drop infl1pctimp_med 
drop infl5pctimp_med

ipolmedian infl1pct, by(yyyy q)
ipolmedian infl5pct, by(yyyy q) 
gen medinfl1pct = infl1pct_med
gen medinfl5pct = infl5pct_med
drop infl1pct_med 
drop infl5pct_med

bysort yyyy q: egen meaninfl1pctimp = mean(infl1pctimp)
bysort yyyy q: egen meaninfl5pctimp = mean(infl5pctimp)

bysort yyyy q: egen meaninfl1pct = mean(infl1pct)
bysort yyyy q: egen meaninfl5pct = mean(infl5pct)


save cohfile, replace

use macrofile, clear
gen stryyyyq = string(yyyy)+"q"+string(q)
gen yqdate = quarterly(stryyyyq, "yq")
format yqdate %tq
list yqdate in 1/30
tsset yqdate
tssmooth ma mainfl = yrinfl, window(3 1) 
keep if yyyy > 1876
tsline mainfl, clcolor(blue) clpattern(solid) clwidth(medthick) ///
        xtitle("Quarter") ytitle("Inflation") ///
             graphregion(color(white))
graph export "`outdir'Infl_hist.eps", replace

keep if yyyy > 1950
tsline yrinfl, clcolor(blue) clpattern(solid) clwidth(medthick) ///
        xtitle("Quarter") ytitle("Inflation") ///
             graphregion(color(white))
graph export "`outdir'Infl_hist_short.eps", replace


use cohfile, clear

gen infl5cat_2_young = .
replace infl5cat_2_young = infl5cat_2 if age < 40
gen infl5cat_2_old = .
replace infl5cat_2_old = infl5cat_2 if age > 60

gen infl5cat_1_young = .
replace infl5cat_1_young = infl5cat_1 if age < 40
gen infl5cat_1_old = .
replace infl5cat_1_old = infl5cat_1 if age > 60

gen infl5cat_0_young = .
replace infl5cat_0_young = infl5cat_0 if age < 40
gen infl5cat_0_old = .
replace infl5cat_0_old = infl5cat_0 if age > 60

gen infl5pct_young = .
replace infl5pct_young = infl5pct if age < 40
gen infl5pct_old = .
replace infl5pct_old = infl5pct if age > 60

gen infl5pctimp_young = .
replace infl5pctimp_young = infl5pctimp if age < 40
gen infl5pctimp_old = .
replace infl5pctimp_old = infl5pctimp if age > 60

//gen err5_young = .
//replace err5_young = err5 if age < 40
//gen err5_old = .
//replace err5_old = err5 if age > 60

//gen sqerr5_young = .
//replace sqerr5_young = err5^2 if age < 40
//gen sqerr5_old = .
//replace sqerr5_old = err5^2 if age > 60


gen infl1cat_2_young = .
replace infl1cat_2_young = infl1cat_2 if age < 40
gen infl1cat_2_old = .
replace infl1cat_2_old = infl1cat_2 if age > 60

gen infl1cat_1_young = .
replace infl1cat_1_young = infl1cat_1 if age < 40
gen infl1cat_1_old = .
replace infl1cat_1_old = infl1cat_1 if age > 60

gen infl1cat_0_young = .
replace infl1cat_0_young = infl1cat_0 if age < 40
gen infl1cat_0_old = .
replace infl1cat_0_old = infl1cat_0 if age > 60

gen infl1pct_young = .
replace infl1pct_young = infl1pct if age < 40
gen infl1pct_old = .
replace infl1pct_old = infl1pct if age > 60
gen infl1pct_mid = .
replace infl1pct_mid = infl1pct if (age < 61 & age > 39)

gen infl1pctimp_young = .
replace infl1pctimp_young = infl1pctimp if age < 40
gen infl1pctimp_old = .
replace infl1pctimp_old = infl1pctimp if age > 60
gen infl1pctimp_mid = .
replace infl1pctimp_mid = infl1pctimp if (age < 61 & age > 39)

//gen err1_young = .
//replace err1_young = err1 if age < 40
//gen err1_old = .
//replace err1_old = err1 if age > 60

//gen sqerr1_young = .
//replace sqerr1_young = err1^2 if age < 40
//gen sqerr1_old = .
//replace sqerr1_old = err1^2 if age > 60

//qui gen spfloc = (yyyy-1872)*4+q 
//qui gen spf = spf[spfloc,1]/100
//qui gen spf10 = spf10[spfloc,1]/100


/*
collapse (mean) infl1pct_mid infl1pctimp_mid infl5pct_young infl5pct_old  ///
                err1_young sqerr1_young err1_old sqerr1_old finfl5 imp5flag meaninfl5pct infl5pctimp_old ///
                infl5cat_2_young infl5cat_2_old infl5cat_1_young infl5cat_1_old ///
                infl5cat_0_young infl5cat_0_old infl5pctimp_young meaninfl5pctimp ///
                infl1pct_young infl1pct_old finfl imp1flag meaninfl1pct infl1pctimp_old ///
                infl1cat_2_young infl1cat_2_old infl1cat_1_young infl1cat_1_old ///
                infl1cat_0_young infl1cat_0_old infl1pctimp_young meaninfl1pctimp spf infl1pct infl5pct ///
                fullmedinfl1pct fullmedinfl5pct, by (yyyy q)

gen rmse1_young = sqrt(sqerr1_young)
gen rmse5_young = sqrt(sqerr5_young)
gen rmse1_old = sqrt(sqerr1_old)
gen rmse5_old = sqrt(sqerr5_old)
*/
gen dinfl1pct_young = infl1pct_young-meaninfl1pct
gen dinfl1pctimp_young = infl1pctimp_young-meaninfl1pctimp
gen dinfl1pct_old = infl1pct_old-meaninfl1pct
gen dinfl1pctimp_old = infl1pctimp_old-meaninfl1pctimp
gen dinfl1pct_mid = infl1pct_mid-meaninfl1pct
gen dinfl1pctimp_mid = infl1pctimp_mid-meaninfl1pctimp

gen dinfl5pct_young = infl5pct_young-meaninfl5pct
gen dinfl5pctimp_young = infl5pctimp_young-meaninfl5pctimp
gen dinfl5pct_old = infl5pct_old-meaninfl5pct
gen dinfl5pctimp_old = infl5pctimp_old-meaninfl5pctimp


gen stryyyyq = string(yyyy)+"q"+string(q)
gen yqdate = quarterly(stryyyyq, "yq")
format yqdate %tq
sort yqdate
quietly by yqdate: gen dup = cond(_N==1,0,_n)
drop if dup>1
list yqdate in 1/10
tsset yqdate


//construct 4-qtr moving average of cross-sectional deviations for young and old
tssmooth ma dinfl5pct_young=dinfl5pct_young, window(3 1 0) replace
tssmooth ma dinfl5pct_old=dinfl5pct_old, window(3 1 0) replace
tssmooth ma dinfl1pct_young=dinfl1pct_young, window(3 1 0) replace
tssmooth ma dinfl1pct_old=dinfl1pct_old, window(3 1 0) replace
tssmooth ma dinfl1pct_mid=dinfl1pct_mid, window(3 1 0) replace
tssmooth ma dinfl1pctimp_young=dinfl1pctimp_young, window(3 1 0) replace
tssmooth ma dinfl1pctimp_old=dinfl1pctimp_old, window(3 1 0) replace
tssmooth ma dinfl1pctimp_mid=dinfl1pctimp_mid, window(3 1 0) replace

twoway tsline dinfl1pct_young, clcolor(black) clpattern(solid) clwidth(medthick) ms(dh) mc(black) mlwidth(thin) recast(scatter) || ///
       tsline dinfl1pct_mid, clcolor(red) clpattern(solid) clwidth(medthick) ms(x) mc(red) mlwidth(thin) recast(scatter) || ///
       tsline dinfl1pct_old, clcolor(blue) clpattern(solid) clwidth(medthick) ms(o) mc(blue) mlwidth(thin) recast(scatter) ||, ///
       xtitle("Quarter") ytitle("Expected Inflation") ylabel(-.02(.01).02) ///
       legend(label(1 "Age < 40") label(2 "Age 40 to 60") label(3 "Age > 60") cols(3))  ///
       graphregion(color(white))
graph export "`outdir'Inflexp1_ages_igntro.eps", replace

 
