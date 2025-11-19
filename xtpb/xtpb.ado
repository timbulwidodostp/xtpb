// This program implements Pooled Bewley estimator by 
// A. Chudik, M. H. Pesaran and R. Smith (2023), "Pooled Bewley Estimator of Long-Run Relationships in Dynamic Heterogenous Panels", revised November 2023
// available at https://doi.org/10.1016/j.ecosta.2023.11.001.
// Please report any errors to alexander.chudik@gmail.com
// Code version: 1.3

// Pooled Bewley estimator is suitable for estimation of cointegrating relationship
// in large T panels where the cross-section dimension (N) could be moderate or large
// See Chudik, Pesaran and Smith (2023) for details.

capture program drop xtpb 
capture program drop simulation 
capture mata mata drop PBuncorrected()
capture mata mata drop compute_omega()
capture mata mata drop PBsimul()
capture mata mata drop SR_Param_EST()
capture mata mata drop GEN_Simulated_Data()
capture mata mata drop PBjack()
capture mata mata drop selhs
capture mata mata drop crm
capture mata mata drop pbestim()


program define xtpb, eclass sortpreserve
	version 15.1 // Stata version where this code was tested. It is possible it may run on older versions as well.
	
    syntax varlist(numeric ts) [if] [in] [,  ///
				BIAScorrect(string) ///
				LAGorder(integer 1) ///
				BOOTstrap(string) ///
				FULLdisplay ///
				ERRorcorrect(string) ///
				RESiduals(string) ///
				]
				
				// store the name of the dependent variable in the local macro depvar
				// store names of the regressors in the local macro indepvars. 
				gettoken depvar indepvars : varlist
				_fv_check_depvar `depvar'

				fvexpand `indepvars'
				local cnames `r(varlist)'
				quietly xtset 
				local d_idvar  `r(panelvar)'
				local d_tvar  `r(timevar)'
				tempvar id_t tvar idvar
				// create in program indicators for time and unit variable	
				egen `idvar' = group(`d_idvar')	
				egen `tvar' = group(`d_tvar')
				sort `idvar' `tvar'  
				qui xtset `d_idvar' `d_tvar' 
				
				local biascorrect `biascorrect'
				local lagorder `lagorder'
				local bootstrap `bootstrap'
				local errorcorrect `errorcorrect'
				local residuals `residuals'

				// store bootstrapping suboptions 
				if ("`bootstrap'" != ""){
					simulation `bootstrap'
					
					local csrobust `s(csrobust)'
					local btx `s(btx)'
					local btx_lagorder `s(btx_lagorder)'
					local bcialpha `s(bcialpha)'
					local seed `s(seed)'
					local simulrep `s(simulrep)'
				}			
				

				
/* Clarification on optional inputs:
 	lagorder:   number of lags for the dependent variable and regressors in the level representation
			    - default value is lagorder=1 (one lag in levels representation, which means zero lags for the first-differenced variables in the error-correction representation
			   
	biascorrect: choice of small-T bias-correction
			   - default value is bc="none" (no bias correction is implemented)
			   - bc="jackknife" implements half-panel jackknife bias correction
			   - bc="bootstrap" implements bias correction based on stochastic simulation. Use bootstrap option below to set random seed 
	
	fulldisplay: option to display all regression outputs used for the error correction coefficients. 
	
	errorcorrect: generates a variable that contains the error correction terms
	
	residuals: generates a variable that contains fitted value residuals
			   
	bootstrap(#reps, [csrobust btx btx_lagorder bcialpha seed]): suboptions associated with bootstrapping
				- #reps, required. Defualt to 2000. 
				- csrobust, choice of resampling errors. Allows for arbitrary cross-sectional dependence of errors by resampling column vectors of cross-sectionally stacked residuals 
					(T needs be sufficiently large). Defauly is iid
				- btx,  choice of bootstrapping algorithm for resampling regressors in x
					- default is btx="fixed", which conduct bootstrapping conditional on x, namely regressors x are fixed across the bootstrap replications
					- btx="varx", regressors in x are resampled in bootstrap replications according to a var model in dx
					- btx="varxy",regressors in x are resampled in bootstrap replications according to a var model in dx augmented with lags of dy
				- btx_lagorder, choice of lag order for the marginal model for regressors in x for bootstrapping
					- when not specified, btx_lagorder is set equal to lagorder
				- bcialpha: choice of p-value for 1-p bootstrapped confidence intervals. 
				    - default value if bcialpha=0.95	
					- any value between 0 and 1 can be chosen
				- seed, set for reproducibility. Default is 123456.  
					  
*/
	
marksample touse 
	
//check if moremata package is installed

capture which lmoremata.mlib
if _rc==111 {
	noi {
			di in gr "xtpb command requires moremata package."
			di in gr "Please install moremata package by typing:"
			di in gr ". ssc install moremata"
			}	
}
else {


	// Set seed based on input with default option to ensure replicability
	if ("`bootstrap'" != ""){
		if ("`seed'"!= ""){
			set seed `seed'
			else {
				set seed 123456
			}
		}	
	}
	
	// stata objects needed below to store output
	tempname b obs N min max avg V u temp_matrix CV CI CIss CIjack
	
	local lag = "`lagorder'"
	local RMC = "`simulrep'"
	local btxlg = "`btx_lagorder'"
	
	if ("`btxlg'"=="" | "`btxlg'"=="-1") local lagx = "`lagorder'"
			else local lagx = "`btx_lagorder'"
			

	// define default choices for btu, bci, bc and btx
	if ("`csrobust'" == "") local btu = "iid" // default choice for btu
	if ("`csrobust'" != "") local btu = "csrobust" // alternative choice for btu
	if ("`biascorrect'" == "") local biascorrect = "none" // default choice for bc
    if ("`bci'" == ""){ // default choices for bci (an indicator for producing bootstapped CIs)
		if ("`biascorrect'" == "bootstrap") local bci = "yes"
			else local bci = "no"
		}
	if ("`bootstrap'" != "") local bci = "yes"
	if ("`btx'" == "") local btx = "fixed" // default choice for btx
	if ("`bcialpha'" == "") local bcialpha = 0.95 // default CI level
	if ("`RMC'" == "") local RMC = 2000 // default bootstrap replications
	if ("`residuals'" == "") local residuals = "residuals" // default, any old data not lost if not specified
	if ("`errorcorrect'" == "") local errorcorrect = "error_correct" // default, any old data not lost if not specified
	local bc = "`biascorrect'" //shorthand for biascorrect
	
	
	//Ensure string options contain an allowed string 
	if  ("`biascorrect'" != ""){
		if ("`biascorrect'" != "none" & "`biascorrect'" != "bootstrap" & "`biascorrect'" != "jackknife"){
			noi di in gr "option 'biascorrect()' incorrectly specified"
			error 198 
		}
	}
	if  ("`btx'" != ""){
		if ("`btx'" != "fixed" & "`btx'" != "varx" & "`btx'" != "varxy"){
			noi di in gr "option 'btx()' incorrectly specified"
			error 198 
		}
	}
	
	//Make sure lagorder is at least 1
	if (`lagorder' < 1){
		noi di in gr "option 'lagorder()' incorrectly specified. 'lagorder()' must be 1 or larger."
		error 198 
		
	}

	// conduct estimations	
	if ("`bc'" == "none"){
		qui xtset `d_idvar' `d_tvar' 
		mata: PBuncorrected("`depvar'", "`cnames'","`idvar'", "`touse'", "`constant'", ///
       `lag', "`b'", "`obs'", "`N'", "`min'", "`max'", "`avg'", "`V'")
	   
	noi {
		di _newline
		di in gr "Pooled Bewley Estimation of Long-Run Relationship in Dynamic Heterogenous Panel"
		di in gr "-------------------------------------------------------------------------------"
		di in gr "Group variable: "	as result "`d_idvar'" _column(45) as text "Total number of observations = "`obs'
		di in gr                  _column(45) as text "Number of groups = " `N'
		di in gr                  _column(45) as text "Obs per group: min = " `min'
		di in gr                  _column(60) as text "max = " `max' 
		di in gr                  _column(60) as text "avg = " `avg' 
		di in gr "Original Pooled Bewley (PB) estimator without bias correction."
		di in gr "Inference below conducted based on asymptotic standard errors."
	  }

		matrix colnames `b' = `cnames'
		matrix colnames `V' = `cnames'
		matrix rownames `V' = `cnames'
		ereturn post `b' `V', depname(`depvar') esample(`touse') 
		ereturn local cmd = "xtpb"
		ereturn local panelvar = "`d_idvar'" 
		ereturn local timevar = "`d_tvar'" 
		ereturn scalar N = `N' 
		ereturn scalar Tavg = `avg' 
		ereturn scalar Tmin = `min' 
		ereturn scalar Tmax = `max'
		ereturn display
		matrix beta = e(b)
		matrix variance = e(V)

		if ("`bci'" == "no"){
			noi {
			di in gr "Bootstrap confidence intervals (CIs) were not computed."
			di in gr "To compute bootstrap CIs, use option 'bootstrap()' and associated suboptions."
			di in gr "-----------------------------------------------------------------------------"
			}
		}
		else {
		noi {
			di in gr "Computing bootstrap CIs based on " `RMC' " bootstrap replications..."
			}
			
		marksample touse 
		qui xtset `idvar' `tvar'
		mata: PBsimul("`depvar'", "`cnames'","`idvar'", "`tvar'", "`touse'", "`constant'", ///
       `lag', `RMC', "`btu'", "`b'", "`obs'", "`N'", "`min'", "`max'", "`avg'","`V'", ///
	    "`CI'", "`CIss'", "`CIjack'", `bcialpha', "`bc'" , "`btx'", `lagx' )
		
		noi {
			//di in gr "... bootstrap computations are finished."
			di in gr _newline
			di in gr `bcialpha'*100 " percent bootstrapped confidence intervals:" 
			}
	    matrix rownames `CI' = `cnames'
		matlist `CI', names(rows)  
		ereturn mat BCI = `CI'
		ereturn scalar bcialpha=`bcialpha'
		matrix bci = e(BCI)
		scalar bcialpha=`bcialpha'
		if ("`btu'" == "iid") {
		 noi { 
		   // di in gr _newline
			di in gr "Bootstrapping assumed no cross sectional dependence of errors."
			di in gr "To allow for arbitrary cross sectional dependence, use option 'btu(csrobust)'"
			di in gr "To change confidence interval coverage, use suboption 'bcialpha()' of option 'bootstrap()'"
			di in gr "-----------------------------------------------------------------------------"
			}
		 }
		if ("`btu'" == "csrobust") {
		noi {
		   // di in gr _newline
		    di in gr "Bootstrapping allowed for arbitrary cross sectional dependence of errors."
			di in gr "To change confidence interval coverage, use suboption 'bcialpha()' of option 'bootstrap()'"
			di in gr "-----------------------------------------------------------------------------"
			}	
		 }
		}
	  }
	
	else if ("`bc'" == "jackknife"){
		qui xtset `d_idvar' `d_tvar' 
		mata: PBjack("`depvar'", "`cnames'","`idvar'", "`touse'", "`constant'", ///
       `lag', "`b'", "`obs'", "`N'", "`min'", "`max'", "`avg'", "`V'")
	   
	   noi {
		di _newline
		di in gr "Pooled Bewley Estimation of Long-Run Relationship in Dynamic Heterogenous Panel"
		di in gr "-------------------------------------------------------------------------------"
		di in gr "Group variable: "	as result "`d_idvar'" _column(45) as text "Total number of observations= " `obs'
		di in gr                  _column(45) as text "Number of groups = " `N'
		di in gr                  _column(45) as text "Obs per group: min = " `min'
		di in gr                  _column(60) as text "max = " `max' 
		di in gr                  _column(60) as text "avg = " `avg' _newline
		di in gr "jackknife bias-corrected Pooled Bewley (PB) estimator."
		di in gr "Inference below conducted based on asymptotic standard errors."
	    }
	
	    matrix colnames `b' = `cnames'
		matrix colnames `V' = `cnames'
		matrix rownames `V' = `cnames'
		ereturn post `b' `V', depname(`depvar')
		ereturn local cmd = "xtpb"
		ereturn local panelvar = "`d_idvar'" 
		ereturn local timevar = "`d_tvar'" 
		ereturn scalar N = `N' 
		ereturn scalar Tavg = `avg' 
		ereturn scalar Tmin = `min' 
		ereturn scalar Tmax = `max'
        ereturn display
		matrix beta = e(b)
		matrix variance = e(V)
		matrix reg_tab = r(table)
		
	if ("`bci'" == "no"){
			noi {
			di in gr "Bootstrap confidence intervals (CIs) were not computed."
			di in gr "To compute bootstrap CIs, use option 'bootstrap()' and associated suboptions."
			di in gr "-----------------------------------------------------------------------------"
			}
		}
	else {
		noi {
			di in gr "Computing bootstrap CIs based on " `RMC' " bootstrap replications..."
			}
		qui xtset `d_idvar' `d_tvar' 
		mata: PBsimul("`depvar'", "`cnames'","`idvar'", "`tvar'", "`touse'", "`constant'", ///
       `lag', `RMC', "`btu'", "`b'", "`obs'", "`N'", "`min'", "`max'", "`avg'","`V'", ///
	    "`CI'", "`CIss'", "`CIjack'", `bcialpha', "`bc'", "`btx'", `lagx' )

        noi {
			//di in gr "... bootstrap computations are finished."
			di in gr _newline
			di in gr `bcialpha'*100 " percent bootstrapped confidence intervals:" 
			}
			
		matrix rownames `CIjack' = `cnames'
		matlist `CIjack', names(rows)  
		ereturn mat BCI = `CIjack'
		ereturn scalar bcialpha=`bcialpha'
		matrix bci = e(BCI)
		scalar bcialpha=`bcialpha'
		if ("`btu'" == "iid") {
		 noi { 
		    //di in gr _newline
			di in gr "Bootstrapping assumed no cross sectional dependence of errors."
			di in gr "To allow for arbitrary cross sectional dependence, use suboption 'csrobust' of option 'bootstrap()'"
			di in gr "To change confidence interval coverage, use suboption 'bcialpha()' of option 'bootstrap()'"
			di in gr "-----------------------------------------------------------------------------"
			}
		 }
		if ("`btu'" == "csrobust") {
		noi {
		   // di in gr _newline
		    di in gr "Bootstrapping allowed for arbitrary cross sectional dependence of errors."
			di in gr "To change confidence interval coverage, use suboption 'bcialpha()' of option 'bootstrap()'"
			di in gr "-----------------------------------------------------------------------------"
			}	
		 }
	 }
	}
	
	else if ("`bc'" == "bootstrap"){
		qui xtset `d_idvar' `d_tvar' 
		mata: PBsimul("`depvar'", "`cnames'","`idvar'", "`tvar'", "`touse'", "`constant'", ///
       `lag', `RMC', "`btu'", "`b'", "`obs'", "`N'", "`min'", "`max'", "`avg'","`V'", ///
	    "`CI'","`CIss'","`CIjack'", `bcialpha', "`bc'", "`btx'", `lagx')
	
		
	   noi {
		di _newline
		di in gr "Pooled Bewley Estimation of Long-Run Relationship in Dynamic Heterogenous Panel"
		di in gr "-------------------------------------------------------------------------------"
		di in gr "Group variable: "	as result "`d_idvar'" _column(45) as text "Total number of observations= " `obs'
		di in gr                  _column(45) as text "Number of groups = " `N'
		di in gr                  _column(45) as text "Obs per group: min = " `min'
		di in gr                  _column(60) as text "max = " `max' 
		di in gr                  _column(60) as text "avg = " `avg' _newline
		di in gr "Bias-corrected Pooled Bewley (PB) estimator using stochastic simulations"
		di in gr "based on " `RMC' " replications."
		di in gr "Inference below conducted based on asymptotic standard errors."
	    }
		
		//matrix list `temp_matrix'
		matrix colnames `b' = `cnames'
		matrix colnames `V' = `cnames'
		matrix rownames `V' = `cnames'
		ereturn post `b' `V', depname(`depvar')
		ereturn local cmd = "xtpb"
		ereturn local panelvar = "`d_idvar'" 
		ereturn local timevar = "`d_tvar'" 
		ereturn scalar N = `N' 
		ereturn scalar Tavg = `avg' 
		ereturn scalar Tmin = `min' 
		ereturn scalar Tmax = `max'
		ereturn display
		
		matrix beta = e(b)
		matrix variance = e(V)
		
		noi {
			//di in gr _newline
			di in gr `bcialpha'*100 " percent bootstrapped confidence intervals:" 
			}
		matrix rownames `CIss' = `cnames'
		matlist `CIss', names(rows)  
		ereturn mat BCI = `CIss'
		ereturn scalar bcialpha=`bcialpha'
		scalar bcialpha=`bcialpha'
		matrix bci = e(BCI)
		if ("`btu'" == "iid") {
		 noi { 
		    //di in gr _newline
			di in gr "Bootstrap CIs are based on " `RMC' " bootstrap replications."
			di in gr "Simulations/bootstrapping assumed no cross sectional dependence of errors."
			di in gr "To allow for arbitrary cross-sectional dependence, use suboption 'csrobust' of option 'bootstrap()'"
			di in gr "To change confidence interval coverage, use suboption 'bcialpha()' of option 'bootstrap()'"
			di in gr "-----------------------------------------------------------------------------"
			}
		 }
		if ("`btu'" == "csrobust") {
		noi {
		    //di in gr _newline
			di in gr "Bootstrap CIs are based on " `RMC' " bootstrap replications."
		    di in gr "Simulations/bootstrapping allowed for arbitrary cross sectional dep. of errors."
			di in gr "To change confidence interval coverage, use suboption 'bcialpha()' of option 'bootstrap()'"
			di in gr "-----------------------------------------------------------------------------"
			}	
		 }	
	}
	//Save regression table
	if ("`bc'" ~= "jackknife"){
	matrix reg_tab = r(table)
	}
	preserve 
	// error correction coeficients
	cap drop pred_`depvar' 
	cap drop `errorcorrect'_temp 
	cap drop `errorcorrect'
	cap drop `residuals'
	
	qui predict pred_`depvar'
	qui gen `errorcorrect' = `depvar' - pred_`depvar'
	qui bys `idvar': egen `errorcorrect'_temp=mean(`errorcorrect')
	qui replace `errorcorrect' = `errorcorrect' - `errorcorrect'_temp
	drop `errorcorrect'_temp
	
	qui levelsof `idvar', local(panel_vals)
	local total_levels = r(r)
	
	
	if ("`fulldisplay'" ~= ""){
		matrix ec = J(`total_levels',1,.)
		matrix rsq = J(`total_levels',2,.)
		foreach i in `panel_vals'{
				qui xtset `d_idvar' `d_tvar' 
				disp "Error correction estimation for panel variable group: `i'"
				qui local lagd= `lag' - 1
				if (`lagd' == 0){
					reg D.`depvar' L.`errorcorrect' D.(`indepvars') if `idvar'==`i', nohe
					qui predict `residuals'`i', residual
					qui replace `residuals'`i'=. if `idvar'~=`i'
					matrix ec[`i',1] = e(b)[1,1]
					matrix err_`i' = r(table)
					matrix rsq[`i',1] = e(r2)
					matrix rsq[`i',2] = e(r2_a)
					matrix rownames rsq = `panel_vals'
				}
				if (`lagd' ~=0){
					reg D.`depvar' L.`errorcorrect' L(1/`lagd').D.`depvar' L(0/`lagd').D.(`indepvars') if `idvar'==`i', nohe
					qui predict `residuals'`i', residual
					qui replace `residuals'`i'=. if `idvar'~=`i'
					matrix ec[`i',1] = e(b)[1,1]
					matrix err_`i' = r(table)
					matrix rsq[`i',1] = e(r2)
					matrix rsq[`i',2] = e(r2_a)
					matrix rownames rsq = `panel_vals'
				}		
		}
		drop pred_`depvar' 
		matrix rownames ec = `panel_vals'
		matrix colnames rsq = r2 r2_a
	}
		
	if ("`fulldisplay'" == "") {
			matrix ec = J(`total_levels',1,.)
			matrix rsq = J(`total_levels',2,.)
			foreach i in `panel_vals'{
				qui xtset `d_idvar' `d_tvar' 
				//disp "Error correction estimation for panel variable group: `i'"
				qui local lagd= `lag' - 1
				if (`lagd' == 0){
					qui reg D.`depvar' L.`errorcorrect' D.(`indepvars') if `idvar'==`i', nohe
					qui predict `residuals'`i', residual
					qui replace `residuals'`i'=. if `idvar'~=`i'
					matrix ec[`i',1] = e(b)[1,1]
					matrix err_`i' = r(table)
					matrix rsq[`i',1] = e(r2)
					matrix rsq[`i',2] = e(r2_a)
					matrix rownames rsq = `panel_vals'
				}
				if (`lagd' ~=0){
					qui reg D.`depvar' L.`errorcorrect' L(1/`lagd').D.`depvar' L(0/`lagd').D.(`indepvars') if `idvar'==`i', nohe
					qui predict `residuals'`i', residual
					qui replace `residuals'`i'=. if `idvar'~=`i'
					matrix ec[`i',1] = e(b)[1,1]
					matrix err_`i' = r(table)
					matrix rsq[`i',1] = e(r2)
					matrix rsq[`i',2] = e(r2_a)
					matrix rownames rsq = `panel_vals'
				}	
			}
			drop pred_`depvar' 
			matrix rownames ec = `panel_vals'
			matrix colnames rsq = r2 r2_a
			matrix rownames rsq = `panel_vals'
	}
	
	qui gen `residuals'=.
	foreach i in `panel_vals'{
		qui replace `residuals'=`residuals'`i' if `idvar'==`i'
		cap drop `residuals'`i'
	}
	
	label variable `errorcorrect' "Error Correction Terms"
	label variable `residuals' "Residuals"
	
	if ("`errorcorrect'" == ""){
		drop `errorcorrect'
	}
	
	if ("`residuals'" == ""){
		drop residuals
	}
	if ("`errorcorrect'" == "" & "`residuals'" == ""){
		restore
	}
	if ("`errorcorrect'" ~= "" | "`residuals'" ~= ""){
		restore, not
	}
	

	// Re-enter all return values 
	ereturn clear
	marksample touse 
	
// 	local beta beta 
// 	local variance variance

	//matrix list beta 

	mata: assignbV("`b'", "`V'", "beta", "variance")
	
	matrix colnames `b' = `cnames'
	matrix colnames `V' = `cnames'
	matrix rownames `V' = `cnames'
	ereturn clear 
	marksample touse 
	ereturn post `b' `V', depname(`depvar')
	
	
	//ereturn post , esample(`touse')
	ereturn local cmd = "xtpb"
	ereturn local panelvar = "`d_idvar'" 
	ereturn local timevar = "`d_tvar'" 
	ereturn local depvar = "`depvar'" 
	ereturn local properties = ""
	
	ereturn scalar N = `N' 
	ereturn scalar Tavg = `avg' 
	ereturn scalar Tmin = `min' 
	ereturn scalar Tmax = `max'
	cap ereturn scalar bcialpha = `bcialpha'
	
	foreach i in `panel_vals'{
		ereturn matrix ec_`i' = err_`i'
	}
	ereturn matrix ec_rsq = rsq
	ereturn matrix ec = ec
	cap ereturn matrix BCI = bci 
	//ereturn matrix beta = beta
	//ereturn matrix Var = variance
	ereturn matrix table = reg_tab
	
	
} // moremata check	

	return clear
	sreturn clear 
end

//Program to define bootstrap suboptions
program simulation, sclass 
	version 15.1
	syntax anything  [, ///
							csrobust ///
							btx(string) ///
							btx_lagorder(integer -1) ///
							bcialpha(real 0.95) ///
							Seed(integer 123456) ///
							]	
							
						sreturn local csrobust `csrobust'
						sreturn local btx `btx'
						sreturn local btx_lagorder `btx_lagorder'
						sreturn local bcialpha `bcialpha'
						sreturn local seed `seed'
						sreturn local simulrep `anything'
end 		


mata: 
void assignbV(string scalar b, // O
			 string scalar V, // O
			 string scalar beta, // I
			 string scalar variance // I 
			 )
	{
		B = st_matrix("beta")
		Vee = st_matrix("variance")
		st_matrix(b, B)
		st_matrix(V, Vee)
	}
end

mata: //PBuncorrected
// function to compute uncorrected PB estimator
void PBuncorrected(string scalar depvar, // I
				string scalar indepvars, // I
				string scalar idvar, // I
				string scalar touse, // I
				string scalar constant, // I   
				real scalar lag_order, // I 
				string scalar new_matrix, // O
				string scalar observations, // O
				string scalar groups, // O
				string scalar minvalue, // O
				string scalar maxvalue, // O
				string scalar avg, // O
				string scalar omg // O
				)
{
		 
		id = st_data(.,idvar,touse) 
		y = st_data(., depvar, touse)
		X = st_data(., indepvars, touse)

		uniqueid = uniqrows(id)
		N = (rows(uniqueid))
		k = cols(X)
		A = J(k, k, 0)
		B = J(k, 1, 0)
		total_obs = J(N, 1, 0)
		p = lag_order
		
		// calling mata function to get bhat
		Bhat = pbestim(X, y, N, k, A, B, id, uniqueid, total_obs, p, half=0)
		
		// calling compute_omega to compute omega
		omega = compute_omega(X, y, N, A, B, id, uniqueid, total_obs, p, Bhat)
		
		// computing min, max and average
		min_value = colmin(total_obs)
		max_value = colmax(total_obs)	
		average = mean(total_obs)
		total_observations = average * N
		
		st_matrix(new_matrix, Bhat)
		st_numscalar(observations, total_observations)
		st_numscalar(groups, N)
		st_numscalar(minvalue, min_value)
		st_numscalar(maxvalue, max_value)
		st_numscalar(avg, average)
		st_matrix(omg, omega)
	
}

end

mata: //compute_omega
// Function to compute variance matrix Omega
function compute_omega(X, // I
					   y, // I
					   N, // I
					   A, // I
					   B, // I
					   id, // I
					   uniqueid, // I
					   total_obs, // I
					   p, // I
					   Bhat // I
					   ){
	// computing omegax
	k = cols(X)
	omegax = A / N
	omegav=J(k, k, 0)
	
	
	// loop for computing variance
	for (i=1; i<=N; i++){
		indic = (id :== uniqueid[i])
		Xi=select(X,indic)
		yi=select(y,indic)
		Ti = rows(Xi)
		total_obs[i, 1] = Ti
		
		// create a vector of 1's 
		obs = Ti - p
		tau = J(obs, 1, 1)
			
		// create an identity matrix
		im = I(obs)
		Mtaui = im - (tau * tau')/(Ti - p)
		
		// computing Xtl and ytl
		Xtli = Mtaui * Xi[p+1..Ti, .]
		ytli = Mtaui * yi[p+1..Ti, .]
		
		// creating zeros matrix for dyi, dxi, yli, and xli
		dyi = J(Ti-p, p, 0)
		dxi = J(Ti-p, p * k, 0)
		xli = J(Ti-p, p * k + k, 0)
		yli = J(Ti-p, p, 0)
		
		// computing dyi and dxi
		for (j=1; j<=p; j++){
			dyi[., j] = yi[p+1+1-j..Ti+1-j, 1] - yi[p+1-j..Ti-j, 1]
			dxi[., (j-1)*k+1..j*k] = Xi[p+1+1-j..Ti+1-j, .] - Xi[p+1-j..Ti-j, .]
			yli[., j] = yi[p+1-j..Ti-j, 1]
		}
		
		// computing xli
		for (j=1; j<=p+1; j++){
			xli[., (j-1)*k+1..j*k] = Xi[p+1+1-j..Ti+1-j, .]
		}
		
		// computing zi and hi
		dzi = dxi, dyi
		hi = xli, yli
			
		// computing htil
		htil = Mtaui * hi
			
		//computing pi and Mi
		htil_inv = invsym(htil' * htil)
		pi = htil * htil_inv * htil'
		
		//computing Mi
		Mi = pi - pi * dzi * invsym(dzi' * pi * dzi) * dzi' * pi
			
		//computing Vi
		Vi = Mi * (ytli - Xtli * Bhat')
			
		//computing w
		Wi = Xtli' * Mi * Vi
			
		//computing omegav
		omegav = omegav + (Wi * Wi') / N
	
	}
	
	// computing inverse of omegax
	invomegax = invsym(omegax)
		
	// computing omega
	omega = invomegax * omegav * invomegax/N
	
	return (omega)
	
}
end

mata: //compute_omegajk
// Function to compute jackknife variance matrix Omegajk
function compute_omegajk(X, // I
					   y, // I
					   N, // I
					   A, // I
					   id, // I
					   uniqueid, // I
					   total_obs, // I
					   p, // I
					   Bhat, // I
					   kappa //I
					   ){
					   
	// computing omegax
	k = cols(X)
	omegax = A / N
	omegav=J(k, k, 0)
	
	
	// loop for computing variance
	for (i=1; i<=N; i++){
		indic = (id :== uniqueid[i])
		Xi=select(X,indic)
		yi=select(y,indic)
		Ti = rows(Xi)
		total_obs[i, 1] = Ti
		
		// create a vector of 1's 
		obs = Ti - p
		tau = J(obs, 1, 1)
	

		// full sample
		Mi=crm(Xi,yi,p)
		
		//sample a
		Xia=selhs(Xi,1)
		yia=selhs(yi,1)
		Mia=crm(Xia,yia,p)
		// sample b
		Xib=selhs(Xi,1)
		yib=selhs(yi,1)
		Mib=crm(Xib,yib,p)
				
		// create an identity matrix
		Tia = rows(Xia)
		Tib = rows(Xib)
		obsa= Tia-p
		obsb= Tib-p
		taua = J(obsa, 1, 1)
		taub = J(obsb, 1, 1)
		
		im =  I(obs)
		ima = I(obsa)
		imb = I(obsb)
		Mtaui = im - (tau * tau')/(Ti - p)
		Mtauia = ima - (taua * taua')/(Tia - p)
		Mtauib = imb - (taub * taub')/(Tib - p)
		
		// computing Xtl and ytl
		Xtli = Mtaui * Xi[p+1..Ti, .]
		ytli = Mtaui * yi[p+1..Ti, .]

		Xtlia = Mtauia * Xia[p+1..Tia, .]
		ytlia = Mtauia * yia[p+1..Tia, .]
		Xtlib = Mtauib * Xib[p+1..Tib, .]
		ytlib = Mtauib * yib[p+1..Tib, .]

		//computing Vi
		Vi = Mi * (ytli - Xtli * Bhat')
		Via = Mia * (ytlia - Xtlia * Bhat')
		Vib = Mib * (ytlib - Xtlib * Bhat')
			
		//computing w

		Wi = (1+kappa)*Xtli' * Mi * Vi
		Wi = Wi-2*kappa*Xtlia' * Mia * Via	
		Wi = Wi-2*kappa*Xtlib' * Mib * Vib	
		
		//computing omegav
		omegav = omegav + (Wi * Wi') / N
	
	}
	
	// computing inverse of omegax
	invomegax = invsym(omegax)
		
	// computing omega
	omega = invomegax * omegav * invomegax/N
	
	return (omega)
	
}
end

mata: //PBsimul
// Function for stochastic simulations and bootstrapping 
void PBsimul(string scalar depvar,  // I
			string scalar indepvars, // I
			string scalar idvar, // I
			string scalar tvar, // I
			string scalar touse, // I
			string scalar constant, // I  
			real scalar lag_order, // I
			real scalar R, // I
			string scalar BTU, // I
			string scalar b_matrix, // O
			string scalar observations, // O
			string scalar groups, // O
			string scalar minvalue, // O
			string scalar maxvalue, // O
			string scalar avg, // O
			string scalar omg, // O
			string scalar CI, // O
			string scalar CIss, // O
			string scalar CIjack, // O			
			real scalar bcialpha, //I
			string scalar bc, // I
			string scalar btx, // I
			real scalar lag_orderx //I
			)
{
	id = st_data(.,idvar,touse) 
	idtime = st_data(.,tvar,touse)
	
	y = st_data(., depvar, touse)
	X = st_data(., indepvars, touse)
	
	uniqueid = uniqrows(id)
	p = lag_order
	px= lag_orderx
	N = (rows(uniqueid))
	k = cols(X)
	A = J(k, k, 0)
	B = J(k, 1, 0)
	total_obs = J(N, 1, 0)
	G = J(N, 2+(p*k)+p-1, 0)
	Gx = J(N*k, 1+(px-1)*k, 0)
	Gxy = J(N*k, 1+(px-1)*(k+1), 0)
	util = J(0, 1, 0)
	utilx = J(0, k, 0)
	utilxy = J(0, k, 0)

	
	
	// get PB estimates
	bhat = pbestim(X, y, N, k, A, B, id, uniqueid, total_obs, p, half=0)
	omegau = compute_omega(X, y, N, A, B, id, uniqueid, total_obs, p, bhat)

	//tsi=diagonal(omega)
	//bhat' :/ tsi :^ J(k,1,0.5)
	
	// calling SR_Param_EST Function
	// https://www.stata.com/statalist/archive/2006-10/msg00633.html
	SR_Param_EST(X, y, N, k, id, uniqueid, total_obs, p, G, util, bhat, btx, Gx, Gxy, utilx, px, utilxy)


	ball =   J(k,R,0)
	tstatsCI = J(k,R,0)




	
	for (r=1; r<=R; r++){
		
		if (r == 1){
			display("Bootstrapping progress:")
		}
		if (round((r/R)*100) == 20 & round(((r+1)/R)*100) ~= 20) {
			display("----------10%")
		}
		if (round((r/R)*100) == 40 & round(((r+1)/R)*100) ~= 40){
			display("----------20%")
		}
		if (round((r/R)*100) == 60 & round(((r+1)/R)*100) ~= 60){
			display("----------30%")
		}
		if (round((r/R)*100) == 80 & round(((r+1)/R)*100) ~= 80){
			display("----------40%")
		}
		if ((r/R) == 1){
			display("----------50%")
		}
		
		Xrs    = J(0, k, 0)
		yrs = GEN_Simulated_Data(X, y, N, k, id,idtime, BTU, uniqueid, total_obs, ///
								p, G, util, bhat, Gx, Gxy, utilx, px, utilxy, btx, Xrs)
		

		
		br = pbestim(Xrs, yrs, N, k, A, B, id, uniqueid, total_obs, p, half=0)
		omegar = compute_omega(Xrs, yrs, N, A, B, id, uniqueid, total_obs, p, br)
		vi=diagonal(omegar)
		ball[.,r] = br'
		tstatsCI[.,r]=abs( (br'- bhat' ):/ (vi :^ J(k,1,0.5)))
	}	
	meanba = ball * J(R,1,1)/R
	biasest = meanba - bhat'
	//biasest
	//bhat'
	//meanba

	bhat_ssimul = bhat'- biasest
	
	// compute bootstrap CI for uncorrected PB estimator
    Q=mm_quantile(tstatsCI', 1, bcialpha)
	vs=diagonal(omegau)
	se= vs :^ J(k,1,0.5)
	CInt=J(k,2,0)
	CInt[.,1]=(bhat- Q :* se')'
	CInt[.,2]=(bhat+ Q :* se')'
	
	// compute asy se for bootstrap bias-corrected estimator
	bhat = pbestim(X, y, N, k, A, B, id, uniqueid, total_obs, p, half=0) // to get A, which serves as input to compute_omega
	omegass = compute_omega(X, y, N, A, B, id, uniqueid, total_obs, p, bhat_ssimul')


	// compute bootstrapped CI for simul-based bias-correction
	if (bc == "bootstrap") {
		tstatsCIss = J(k,R,0)
		for (r=1; r<=R; r++){
		
// 		if (r == 1){
// 			display("Bootstrapping progress:")
// 		}
		if (round((r/R)*100) == 20 & round(((r+1)/R)*100) ~= 20){
			display("----------60%")
		}
		if (round((r/R)*100) == 40 & round(((r+1)/R)*100) ~= 40){
			display("----------70%")
		}
		if (round((r/R)*100) == 60 & round(((r+1)/R)*100) ~= 60){
			display("----------80%")
		}
		if (round((r/R)*100) == 80 & round(((r+1)/R)*100) ~= 80){
			display("----------90%")
		}
		if ((r/R) == 1){
			display("----------Complete!")
		}
		
			Xrs    = J(0, k, 0)
			yrs = GEN_Simulated_Data(X, y, N, k, id,idtime, BTU, uniqueid, total_obs, ///
									 p, G, util, bhat, Gx, Gxy, utilx, px, utilxy, btx, Xrs)
			br = pbestim(Xrs, yrs, N, k, A, B, id, uniqueid, total_obs, p, half=0)
			omegar = compute_omega(Xrs, yrs, N, A, B, id, uniqueid, total_obs, p, br)
			vi=diagonal(omegar)
			
			tstatsCIss[.,r]=abs( (br'-biasest- bhat' ) :/ (vi :^ J(k,1,0.5)))
			}	
	    Qss=mm_quantile(tstatsCIss', 1, bcialpha)
		vs=diagonal(omegass)
		se= vs :^ J(k,1,0.5)
		CIntss=J(k,2,0)
		CIntss[.,1]=(bhat_ssimul'- Qss :* se')'
		CIntss[.,2]=(bhat_ssimul'+ Qss :* se')'
		
	}
	if (bc == "jackknife") {
	
	   // compute jackknife estimator
	   Bhat_tmp = pbestim(X, y, N, k, Afs, B, id, uniqueid, total_obs, lag_order, half=0)
	   Bhatl_tmp = pbestim(X, y, N, k, A, B, id, uniqueid, total_obs, lag_order, half=1)
	   Bhatr_tmp = pbestim(X, y, N, k, A, B, id, uniqueid, total_obs, lag_order, half=2)
	   
	   //computing jackknife estimator
	   K = 1/3
	   Bhatjk = Bhat_tmp - K * ((Bhatl_tmp + Bhatr_tmp) / 2 - Bhat_tmp)
	   omegajk = compute_omegajk(X, y, N, Afs, id, uniqueid, total_obs, lag_order, Bhatjk, K)		


	
		tstatsCIjack = J(k,R,0)
		for (r=1; r<=R; r++){
		
		// 		if (r == 1){
// 			display("Bootstrapping progress:")
// 		}
		if (round((r/R)*100) == 20 & round(((r+1)/R)*100) ~= 20){
			display("----------60%")
		}
		if (round((r/R)*100) == 40 & round(((r+1)/R)*100) ~= 40){
			display("----------70%")
		}
		if (round((r/R)*100) == 60 & round(((r+1)/R)*100) ~= 60){
			display("----------80%")
		}
		if (round((r/R)*100) == 80 & round(((r+1)/R)*100) ~= 80){
			display("----------90%")
		}
		if ((r/R) == 1){
			display("----------Complete!")
		}

			Xrs    = J(0, k, 0)
			yrs = GEN_Simulated_Data(X, y, N, k, id,idtime, BTU, uniqueid, total_obs, ///
								     p, G, util, Bhatjk, Gx, Gxy, utilx, px, utilxy, btx, Xrs)
			Bhat_r = pbestim(Xrs, yrs, N, k, Afs, B, id, uniqueid, total_obs, lag_order, half=0)
	        Bhatl_r = pbestim(Xrs, yrs, N, k, A, B, id, uniqueid, total_obs, lag_order, half=1)
	        Bhatr_r = pbestim(Xrs, yrs, N, k, A, B, id, uniqueid, total_obs, lag_order, half=2)
	   	    Bhatjkr = Bhat_r - K * ((Bhatl_r + Bhatr_r) / 2 - Bhat_r)
 	        omegajkr = compute_omegajk(Xrs, yrs, N, Afs, id, uniqueid, total_obs, lag_order, Bhatjkr, K)
			vi=diagonal(omegajkr)
			
			tstatsCIjack[.,r]=abs( (Bhatjkr'-Bhatjk' ) :/ (vi :^ J(k,1,0.5)))
			}
			
	    Qjack=mm_quantile(tstatsCIjack', 1, bcialpha)
		vs=diagonal(omegajk)
		se= vs :^ J(k,1,0.5)
		
		CIntjack=J(k,2,0)
		CIntjack[.,1]=(Bhatjk- Qjack :* se')'
		CIntjack[.,2]=(Bhatjk+ Qjack :* se')'
		
	}	
	
	
	// computing min, max and average
	min_value = colmin(total_obs)
	max_value = colmax(total_obs)	
	average = mean(total_obs)
	total_observations = average * N
		
	st_matrix(b_matrix, bhat_ssimul')
	st_matrix(omg, omegass)
	st_numscalar(observations, total_observations)
	st_numscalar(groups, N)
	st_numscalar(minvalue, min_value)
	st_numscalar(maxvalue, max_value)
	st_numscalar(avg, average)
	st_matrix(CI, CInt)
	st_matrix(CIss, CIntss)
	st_matrix(CIjack, CIntjack)
}
				
end

mata: //SR_Param_EST
// function for estimation of short-run parameters
void SR_Param_EST(X, // I
				  y, // I
				  N, // I
				  k, // I
				  id, // I
				  uniqueid, // I
				  total_obs, // I
				  p, // I
				  G, // I/O
				  util, // I/O
				  bhat, // I
				  btx, //I
				  Gx, // I/O
				  Gxy, // I/O
				  utilx, // I/O
				  px,  // I
				  utilxy // I/O
				  )
{

	for (i=1; i<=N; i++){
		    indic = (id :== uniqueid[i])
			Xi=select(X,indic)
			yi=select(y,indic)
			Ti = rows(Xi)
			
			total_obs[i, 1] = Ti		
			
			yi1 = yi[p..(Ti - 1), .]
			Xi1 = Xi[p..(Ti - 1), .]
			
			// computing episoloni
			ei = yi1 - Xi1 * bhat'
			
			// create a vector of 1's 
			Tis = Ti - p
			Tisx = Ti - px
			
			tau = J(Tis, 1, 1)
			taux = J(Tisx, 1, 1)
			
			// creating zeros matrix for dyi, dxi
			dyi  = J(Tis, p-1, 0)
			dyi0 = J(Tis, 1, 0)
			dXi  = J(Tis, p * k, 0)
			
			dXi0 = J(Tisx, k, 0)
			dyix  = J(Tisx, px-1, 0)
			dXix  = J(Tisx, (px-1) * k, 0)
			
			// computing dyi and dxi
			dyi0 = yi[p+1..Ti, .] - yi[p..Ti-1, .]
			dXi0 = Xi[px+1..Ti, .] - Xi[px..Ti-1, .]
			dXi[., 1..k] = Xi[p+1..Ti,.] - Xi[p..Ti-1,.]
			
			for (j=1; j<=p-1; j++){
				dyi[., j] = yi[p-j+1..Ti-j, .] - yi[p-j..Ti-j-1, .]
				dXi[., j*k+1..(j+1)*k] = Xi[p-j+1..Ti-j,.] - Xi[p-j..Ti-j-1,.]
			}
			
			for (j=1; j<=px-1; j++){
				dyix[., j] = yi[px-j+1..Ti-j, .] - yi[px-j..Ti-j-1, .]
				dXix[., (j-1)*k+1..j*k] = Xi[px-j+1..Ti-j,.] - Xi[px-j..Ti-j-1,.]
			}
 
			// Constructing R
			R = J(Tis, 2+(p*k)+p-1, 0)
			R[., 1] = tau
			R[., 2] = ei
			R[., 3..p*k+2] = dXi
			if (p > 1){
				R[., p*k+3..p-1+p*k+3-1] = dyi
			}
				
			// Computing ghati
			ghati = invsym(R' * R) * R' * dyi0
			
			// Computing uhati
			uhati = dyi0 - R * ghati
			
			// Computing utili
			utili = J(Ti, 1, 0)
			utili[1..p, 1] = J(p, 1, 0)
			utili[p+1..Ti, 1] = uhati

			
			// Computing utili
			util = util\utili
			
			// Computing G
			G[i, .] = ghati'

			if (btx == "varx") {
				Rx =J(Tisx, 1+(px-1)*k, 0)
				Rx[., 1] = taux
				if (px > 1){
					Rx[., 2..(px-1)*k+1] = dXix
				}
			
				ghatix = invsym(Rx' * Rx) * Rx' * dXi0
				uhatix = dXi0 - Rx * ghatix
				
				utilix = J(Ti, k, 0)
				utilix[1..px, 1..k] = J(px, k, 0)
				utilix[px+1..Ti, 1..k] = uhatix
				utilx = utilx\utilix
				Gx[(i-1)*k+1..i*k, .] = ghatix'
			}
			
			if (btx == "varxy") {

				Rxy =J(Tisx, 1+(px-1)*k+px-1, 0)
				Rxy[., 1] = taux

				if (px > 1){
					Rxy[., 2..(px-1)*k+1] = dXix
					Rxy[., (px-1)*k+2..px*k] = dyix
				}			
				ghatixy = invsym(Rxy' * Rxy) * Rxy' * dXi0
				uhatixy = dXi0 - Rxy * ghatixy
				
				utilixy = J(Ti, k, 0)
				utilixy[1..px, 1..k] = J(px, k, 0)

				utilixy[px+1..Ti, 1..k] = uhatixy
				utilxy = utilxy\utilixy
				Gxy[(i-1)*k+1..i*k, .] = ghatixy'
			}
			
			
	}
	//Gx
}

end

mata: //GEN_Simulated_Data
// function for generating simulated sample (conditional on X)
function GEN_Simulated_Data(X, // I
							y, // I
							N, // I
							k, // I
							id, // I
							idtime, // I
							BTU, // I
							uniqueid, // I
							total_obs, // I
							p, // I
							G, // I
							util, 	// I
							bhat, 	// I
							Gx, 	// I
							Gxy, 	// I
							utilx, 	// I
							px, 	// I
							utilxy, // I
							btx,    // I
							Xrs		// I/O
							){
	yrs = J(0, 1, 0)
	// creating a vector with id's from 1 to max time period
	max = colmax(idtime)

	a_timeid = J(max, 1, 0)
	for (j=1; j<=max; j++){
		a_timeid[j, 1] = j
	}
	
	// create a big vector a with dimensions same as a_timeid with values +1 and -1
	t = rows(a_timeid)
	a_prob = J(2, 1, 0.5)
	a = rdiscrete(t, 1, a_prob)
	v = J(t, 1, 1.5)
	a = a - v
	a = a * 2
	

	for (i=1; i<=N; i++){
		indic = (id :== uniqueid[i])
		Xi=select(X,indic)
		yi=select(y,indic)
		//Xi
		//yi
		Ti = rows(Xi)
		gi=G[i,.]
		gix=Gx[(i-1)*k+1..i*k,.]
		gixy=Gxy[(i-1)*k+1..i*k,.]
		uir = J(Ti, 1, 0)

		
		prob = J(2, 1, 0.5)
		ai = rdiscrete(Ti, 1, prob)
		vi = J(Ti, 1, 1.5)
		ai = ai - vi
		ai = ai * 2
		//ai= J(Ti, 1, 1) //just for checking!
		ui = select(util,indic)
		ui_timeid = select(idtime,indic)
		
		if (BTU == "iid") {
			for (j=1; j<=Ti; j++){
				uir[j, 1] = ui[j, 1] * ai[j, 1] // btu=iid option - default
			}
		}

		if (BTU == "csrobust") {
			for (j=1; j<=Ti; j++){
				val = ui_timeid[j, 1]
				indica = (a_timeid :== val)
				at = select(a,indica)
				uir[j, 1] = at * ui[j, 1]
			}	
		}
		
		if (btx == "varx") {
			uirx = J(Ti, k, 0)
			uix = select(utilx,indic)
			if (BTU == "iid") {
				for (j=1; j<=Ti; j++){
					uirx[j, .] = uix[j, .] * ai[j, 1] 
				}
			}
			if (BTU == "csrobust") {
				for (j=1; j<=Ti; j++){
					val = ui_timeid[j, 1]
					indica = (a_timeid :== val)
					at = select(a,indica)
					uirx[j, .] = uix[j, .] * at
					//at
				}	
			} 
		} //if btx
		if (btx == "varxy") {
			uirxy = J(Ti, k, 0)
			uixy = select(utilxy,indic)
			if (BTU == "iid") {
				for (j=1; j<=Ti; j++){
					uirxy[j, .] = uixy[j, .] * ai[j, 1] 
				}
			}
			if (BTU == "csrobust") {
				for (j=1; j<=Ti; j++){
					val = ui_timeid[j, 1]
					indica = (a_timeid :== val)
					at = select(a,indica)
					uirxy[j, .] = at * uixy[j, .]
				}	
			}
		} //if btx

		
		yir = J(Ti, 1, 0)
		for (t=1; t<=p; t++){
			yir[t,1] = yi[t,1]
		}
		Tis = Ti - p
		Tisx = Ti - px
		
		dyir  = J(Tis, p-1, 0)
		dXir  = J(Tis, p * k, 0)
		Xir   = J(Ti, k, 0)
		
		if (btx== "fixed") {
			dXir[., 1..k] = Xi[p+1..Ti,.] - Xi[p..Ti-1,.]
			Xir=Xi
			for (j=1; j<=p-1; j++){
				dXir[., j*k+1..(j+1)*k] = Xi[p-j+1..Ti-j,.] - Xi[p-j..Ti-j-1,.]
			} //j		
		}
		else {
			Xir[1..px, .] = Xi[1..px,.] // initial values
			
			if (px<p+1) { // generate X up to p+1 
				for (t=px+1; t<=p+1; t++) {
					if (btx== "varx") {
						h=J(1+(px-1)*k,1,0)
						h[1,1]=1
						q=J((px-1)*k,1,0)
						
						if (px > 1){
							for (j=1; j<=px-1; j++){	
								q[(j-1)*k+1..j*k,1]=(Xir[t-j,.]-Xir[t-j-1,.])'
							}
							h[2..1+(px-1)*k,1]=q
						}
						Xir[t,.] = Xir[t-1, .]+ h'*gix' + uirx[t, .]
					}
					if (btx== "varxy") {
						h=J(1+(px-1)*k+px-1,1,0)
						h[1,1]=1
						q=J((px-1)*k,1,0)
						qy=J(px-1,1,0)
						if (px > 1){
							for (j=1; j<=px-1; j++){	
								q[(j-1)*k+1..j*k,1]=(Xir[t-j,.]-Xir[t-j-1,.])'
								qy[(j-1)+1..j,1]=(yir[t-j,.]-yir[t-j-1,.])'
							}
							h[2..1+(px-1)*k,1]=q
							h[1+(px-1)*k+1..1+(px-1)*k+px-1,1]=qy
						}
						Xir[t,.] = Xir[t-1, .]+ h'*gixy' + uirxy[t, .]
		
					
					}					
					
						
				} //t
			} 
		}
		
		// Computing yir

		
		for (t=p+1; t<=Ti; t++){ //
			if (btx== "fixed") {
				z = dXir[t-p, .]'
			}
			if (btx== "varx") {
			
				if (t>px) {
					h=J(1+(px-1)*k,1,0)
					h[1,1]=1
					q=J((px-1)*k,1,0)
					if (px > 1){
						for (j=1; j<=px-1; j++){
					
							q[(j-1)*k+1..j*k,1]=(Xir[t-j,.]-Xir[t-j-1,.])'
						}
						h[2..1+(px-1)*k,1]=q
					}
					Xir[t,.] = Xir[t-1, .]+ h'*gix' + uirx[t, .]
				}
				
				z=J( p * k,1, 0)
				z[1..k,1]=(Xir[t,.]-Xir[t-1,.])'
				for (j=1; j<=p-1; j++){	
						z[j*k+1..(j+1)*k,1]=(Xir[t-j,.]-Xir[t-j-1,.])'
				}
			}
			if (btx== "varxy") {
				//z = dXir[t-p, .]'
				
				if (t>px) {
					h=J(1+(px-1)*k+px-1,1,0)
					h[1,1]=1
					q=J((px-1)*k,1,0)
					qy=J(px-1,1,0)
					if (px > 1){
						for (j=1; j<=px-1; j++){
					
							q[(j-1)*k+1..j*k,1]=(Xir[t-j,.]-Xir[t-j-1,.])'
							qy[(j-1)+1..j,1]=(yir[t-j,.]-yir[t-j-1,.])'
						}
						h[2..1+(px-1)*k,1]=q
						h[1+(px-1)*k+1..1+(px-1)*k+px-1,1]=qy
					}
					Xir[t,.] = Xir[t-1, .]+ h'*gixy' + uirxy[t, .]
				}
				
				z=J( p * k,1, 0)
				z[1..k,1]=(Xir[t,.]-Xir[t-1,.])'
				for (j=1; j<=p-1; j++){	
						z[j*k+1..(j+1)*k,1]=(Xir[t-j,.]-Xir[t-j-1,.])'
				}
				
			}
			
			if (p > 1){
				q = J(p-1,1,0)
				for (j=1; j<=p-1; j++){
					q[j, 1] = yir[t-j, .] - yir[t-j-1, .]
				}
			}
			ecm = yir[t-1,.] - Xir[t-1,.]*bhat'
			h=J(2+p*k+p-1,1,0)
			h[1,1]=1
			h[2,1]=ecm
			h[3..2+p*k,1]=z
			if (p > 1){
				h[2+p*k+1..2+p*k+p-1]=q
			}
			yir[t,1] = yir[t-1, 1]+ gi*h + uir[t, 1]
			//t
			//Xir 
			
		} //t

		yrs = yrs\yir
		Xrs = Xrs\Xir
//Xir-Xi

		
	} //i
//yrs>1
    return (yrs)
	
}

end

mata: //PBjack
// function for computation of jackknife corrected PB estimator
void PBjack(string scalar depvar,    // I
				string scalar indepvars, // I
				string scalar idvar,     // I
				string scalar touse,     // I
				string scalar constant,  // I  
				real scalar lag_order,   // I
				string scalar new_matrix, // O
				string scalar observations, // O
				string scalar groups, // O
				string scalar minvalue, // O
				string scalar maxvalue, // O
				string scalar avg, // O
				string scalar omg // O
				) 
{			     
		 
		id = st_data(.,idvar,touse) 
		y = st_data(., depvar, touse)
		X = st_data(., indepvars, touse)

		uniqueid = uniqrows(id)
		N = (rows(uniqueid))
		k = cols(X)
		A = J(k, k, 0)
		B = J(k, 1, 0)
		total_obs = J(N, 1, 0)
		
		
	   Bhat = pbestim(X, y, N, k, Afs, B, id, uniqueid, total_obs, lag_order, half=0)
	   Bhatl = pbestim(X, y, N, k, A, B, id, uniqueid, total_obs, lag_order, half=1)
	   Bhatr = pbestim(X, y, N, k, A, B, id, uniqueid, total_obs, lag_order, half=2)
	   
	   //computing jackknife estimator
	   K = 1/3
	   Bhatjk = Bhat - K * ((Bhatl + Bhatr) / 2 - Bhat)

	   omegajk = compute_omegajk(X, y, N, Afs, id, uniqueid, total_obs, lag_order, Bhatjk, K)		



		// computing min, max and average
		min_value = colmin(total_obs)
		max_value = colmax(total_obs)	
		average = mean(total_obs)
		total_observations = average * N
		
		//st_matrix(new_matrixB, omega)
		st_matrix(new_matrix, Bhatjk)
		st_numscalar(observations, total_observations)
		st_numscalar(groups, N)
		st_numscalar(minvalue, min_value)
		st_numscalar(maxvalue, max_value)
		st_numscalar(avg, average)
		st_matrix(omg, omegajk)

 }
			 
end

mata: //selhs
// function to select subsample (full, first half, second half)
function selhs(X, // I
				half // I
				  ){
 	T = rows(X)
			if (half == 0) Xs=X; //full sample
			
			if (half == 1){ // first half subsample
				// check if number of rows are even
				if (mod(T, 2) == 0) Xs=X[1..(T / 2), .]
				else Xs=X[1..(T + 1)/ 2, .]
			}

			if (half == 2){ // second half subsample
				// check if number of rows are even
				if (mod(T, 2) == 0)	Xs=X[(T / 2)+1..T, .]
				else Xs=X[(T + 1)/ 2..T, .]
			}			
	return (Xs)
 }
end

mata: //crm
// function to create Mi matrix 
function crm(Xi, // I
			  yi, // I
			  p	  // I
				  ){
	Ti = rows(Xi)
	k  = cols(Xi)
	// create a vector of 1's 
	obs = Ti - p
	tau = J(obs, 1, 1)
			
	// create an identity matrix
	im = I(obs)
	Mtaui = im - (tau * tau')/(Ti - p)
			
	// computing Xtl and ytl
	Xtli = Mtaui * Xi[p+1..Ti, .]
	ytli = Mtaui * yi[p+1..Ti, .]
			
	// creating zeros matrix for dyi, dxi, yli, and xli
	dyi = J(Ti-p, p, 0)
	dxi = J(Ti-p, p * k, 0)
	xli = J(Ti-p, p * k + k, 0)
	yli = J(Ti-p, p, 0)
			
	// computing dyi and dxi
	for (j=1; j<=p; j++){
		dyi[., j] = yi[p+1+1-j..Ti+1-j, 1] - yi[p+1-j..Ti-j, 1]
		dxi[., (j-1)*k+1..j*k] = Xi[p+1+1-j..Ti+1-j, .] - Xi[p+1-j..Ti-j, .]
		yli[., j] = yi[p+1-j..Ti-j, 1]
	}
			
	// computing xli
	for (j=1; j<=p+1; j++){
		xli[., (j-1)*k+1..j*k] = Xi[p+1+1-j..Ti+1-j, .]
	}
			
	// computing zi and hi
	dzi = dxi, dyi
	hi = xli, yli
			
	// computing htil
	htil = Mtaui * hi
			
	//computing pi and Mi
	htil_inv = invsym(htil' * htil)
	pi = htil * htil_inv * htil'
			
	// computing Mi
	Mi = pi - pi * dzi * invsym(dzi' * pi * dzi) * dzi' * pi 
	return (Mi)
 }
end

mata: //pbestim
// function for PB estimator - point estimates only
function pbestim(X, // I
				  y, // I
				  N, // I
				  k, // I
				  A, // O
				  B, // O
				  id, // I
				  uniqueid, // I
				  total_obs, // I
				  lag_order, // I
				  half // I
				  ){
 
 A = J(k, k, 0) // initialize A
 B = J(k, 1, 0) // initialize B
 
	for (i=1; i<=N; i++){
		    indic = (id :== uniqueid[i])
			Xi=select(X,indic)
			yi=select(y,indic)
			
			
			Xi=selhs(Xi,half)
			yi=selhs(yi,half)
			Ti = rows(Xi)
		
			total_obs[i, 1] = Ti
			p = lag_order
			
			// create a vector of 1's 
			obs = Ti - p
			tau = J(obs, 1, 1)
			
			// create an identity matrix
			im = I(obs)
			Mtaui = im - (tau * tau')/(Ti - p)
			
			// computing Xtl and ytl
			Xtli = Mtaui * Xi[p+1..Ti, .]
			ytli = Mtaui * yi[p+1..Ti, .]
			
			Mi= crm(Xi,yi,p)
			
			//computing Ai and Bi
			Ai = Xtli' * Mi * Xtli
			Bi = Xtli' * Mi * ytli
			
			//computing A and B
			A = A + Ai
			B = B + Bi
		}
			
		// computing Bhat
		Bhat = (invsym(A) * B)'
		return (Bhat)
 }

end
