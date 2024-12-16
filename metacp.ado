capture program drop metacp
program define metacp
    version 13.0
    syntax varlist(numeric),method(string) [ corr_matrix(varlist) ] 
    tempvar num_vars t_sum z_mean z_value chi_squared_sum z_sum invchi_squared_sum T_sum t min_val C_sum c h C
    egen `num_vars' = rownonmiss(`varlist')
	

if "`method'" == "logitp" {
	
	* Define the numerical value of pi
	local pi_value = 3.141592653589793
	
	* Calculate the constant C
	local C = sqrt(`num_vars' * `pi_value'^2 * (5 * `num_vars' + 2) / (3 * (5 * `num_vars' + 4)))
	
	* Loop through each variable
	foreach var of varlist `varlist' {
		* Calculate the t-value for each observation
		gen t_value_`var' = ln(`var' / (1 - `var'))
	}
	
	* Sum the t-values for each observation
	egen `t_sum' = rowtotal(t_value_*)
	
	* Calculate the combined p-value
	gen _logitp = 2 * (1 - t(2 * `num_vars', abs(-`t_sum' /`C')))
	
	* Display the combined p-value
	di "Combined p-value: " _logitp
}
else if "`method'" == "meanp" {	
	
	
	* Define the numerical value of pi
	local pi_value = 3.141592653589793
	

	egen `z_mean' = rowmean(`varlist')   
	gen `z_value' = (0.5 - `z_mean')* sqrt(12 * `num_vars')  

	* Calculate combined p-value
	gen _meanp = 1 - normal(`z_value') 
	di "Combined p-value:", _meanp
}
else if "`method'" == "fisherp" {


	foreach var of varlist `varlist' {
		* Calculate the t-value for each observation
		gen log_`var' = ln(`var')
	}

	egen `chi_squared_sum' = rowtotal(log_*)   

	* Calculate combined p-value
	gen _fisher_p = 1 - chi2(2 * `num_vars', -2 * `chi_squared_sum') 
	di "Combined p-value:", _fisher_p
}
else if "`method'" == "stoufferp" {
	 * Define the numerical value of pi
	local pi_value = 3.141592653589793
	foreach var of varlist `varlist' {
		* Calculate the t-value for each observation
		gen z_score_`var' = invnorm(1 - `var')
	}

	egen `z_sum' = rowtotal(z_score_*)    

	* Calculate combined p-value
	gen _stouffer_p = 1 - normal(`z_sum' / sqrt(`num_vars')) 
	di "Combined p-value:", _stouffer_p
}
else if "`method'" == "invchi2p" {	
	* Loop through each variable in the varlist
	foreach var of varlist `varlist' {
		* Calculate the inverse chi-squared statistic for each observation
		gen invchi_squared_`var' = invchi2(1, 1 - `var')
	}

	* Calculate the sum of inverse chi-squared statistics across variables
	egen `invchi_squared_sum' = rowtotal(invchi_squared_*)   

	* Calculate combined p-value
	gen _inv_chi2_p = 1 - chi2(`num_vars', `invchi_squared_sum') 

	* Display combined p-value
	di "Combined p-value:", _inv_chi2_p
}

else if "`method'" == "binomialp" {

	local alpha = 0.05
	local a = 0.5
	tempvar num_vars num_of_successes 
	

	
	local `num_of_successes' = 0
	
	* Loop through each variable in the varlist
	foreach var of varlist `varlist' {
		if `var' < `alpha' {
			* Increment count of successes
			local `num_of_successes' = `num_of_successes' + 1
		}
	}
	
	* Initialize combined p-value
	local pmf = 0
	
	* Calculate combined p-value
	forvalues x = `num_of_successes'/`num_vars' {
		* Calculate probability mass function for each possible success count
		local pmf = `pmf' + binomial(`x', `num_vars', `a')
	}
	
	* Generate combined p-value
	gen _binomial_p = `pmf'
	
	* Display combined p-value
	di "Combined p-value:", _binomial_p
}

else if "`method'" == "cauchyp" {
	
	local pi_val = 3.141592653589793
	local gamma_val = 1
	
	* Loop through each variable in the varlist
	foreach var of varlist `varlist' {
		* Calculate the inverse chi-squared statistic for each observation
		gen T_`var' = tan((0.5 - `var') * `pi_val')
	}

	* Calculate the sum of inverse chi-squared statistics across variables
	egen `T_sum' = rowtotal(T_*)   
	
	gen `t' = `T_sum' / `num_vars'    
	* Calculate combined p-value calculateda 1 - cauchy.cdf(t)
	gen _cauchy_p = 1 - (0.5 + atan(`t' / `gamma_val') / `pi_val')
	
	* Display combined p-value
	di "Combined p-value:", _cauchy_p
	
}	

else if "`method'" == "minp" {

	* Find the minimum value across variables in each observation
	egen `min_val' = rowmin(`varlist')
	
	* Calculate combined p-value
	gen _min_p = 1 - (1 - `min_val') ^ `num_vars'
	
	* Display combined p-value
	di "Combined p-value:", _min_p
}

else if "`method'" == "MCM" {
	
	gen _MCM_p = 2 * min(_cauchy_p, _min_p, 0.05)
	* Display combined p-value
	di "Combined p-value by MCM(MinP-Cauchy-MinP):", _MCM_p
}	
*For the CMC method the user should first calculate the '_cauchy_p' and the '_min_p'.
*The CMC method should be run like this 'metacp __cauchy_p _min_p,method(CMC)'
else if "`method'" == "CMC" {
	local pi_val = 3.141592653589793
	local gamma_val = 1
	
	* Loop through each variable in the varlist
	foreach var of varlist `varlist' {
		gen C_`var' = tan((0.5 - `var') * `pi_val')
	}

	* Calculate the sum of inverse chi-squared statistics across variables
	egen `C_sum' = rowtotal(C_*)   
	
	* Calculate combined p-value
	gen `c' = `C_sum' / `num_vars'
	
	* Calculate combined p-value by CMC method (Cauchy-MinP-Cauchy)
	gen _CMC_p = 1 - (0.5 + atan(`c' / `gamma_val') / `pi_val')
	
	* Display combined p-value
	di "Combined p-value by CMC(Cauchy-MinP-Cauchy):", _CMC_p
}	

else if "`method'" == "bonferronip" {
		egen `min_val' = rowmin(`varlist')
	
	 * Calculate combined p-value
	gen _bonferroni_p = min(1, `min_val' * `num_vars')
	
	* Display combined p-value
	di "Combined p-value:", _bonferroni_p
}

else if "`method'" == "bonferronig_p" {
    * Calculate the correlation matrix for the variables in varlist
    corr `varlist'
    matrix corr_matrix = r(C)
    
    local nrows = rowsof(corr_matrix)
    local ncols = colsof(corr_matrix)

    * Initialize max_corr with a very small value
    scalar max_corr = -1

    * Loop to find the maximum correlation value (excluding the diagonal)
    forval i = 1/`nrows' {
        forval j = `=`i' + 1'/`ncols' {  
            local corr_ij = corr_matrix[`i', `j']
            if `corr_ij' > max_corr {
                scalar max_corr = `corr_ij'
            }
        }
    }

    * Display the maximum correlation value
    di "Maximum correlation value:", max_corr
    
    
    scalar g_adjusted = (`num_vars' + 1) - (1 + (`num_vars' - 1) * max_corr)

    egen `min_val' = rowmin(`varlist')
    
    * Calculate combined p-value using Bonferroni correction
    gen _bonferronig_p = min(1, `min_val' * g_adjusted)
    
    * Display combined p-value
    di "Combined p-value:", _bonferronig_p
}
/*
else if "`method'" == "cheverud_nyholtp" {
    local end_mata end 
	mkmat `corr_matrix', matrix(correlation_matrix)
	matrix eigenvalues lambda V = correlation_matrix
	mat list lambda
	mata: 
		lambda = st_matrix("lambda")

		st_matrix("variance", diagonal(variance(lambda')))
    `end_mata'
	
	mat eig_variance = variance

	* Extract the value of the standard deviation
	scalar sample_variance = eig_variance[1,1]
	di "Sample variance of Eigenvalues: " sample_variance
	gen effective_number_of_tests = 1 + (`num_vars' - 1) * (1 - sample_variance / `num_vars')
	egen `min_val' = rowmin(`varlist')
	gen _cheverud_nyholt_p = min(1, `min_val' * effective_number_of_tests)
	di "CN: " _cheverud_nyholt_p
}

else if "`method'" == "li_jip" {
	local end_mata end 
	mkmat `corr_matrix', matrix(correlation_matrix)
	matrix eigenvalues lambda V = correlation_matrix
	mat list lambda
	mata: 
		lambda = st_matrix("lambda")
		abs_eig = abs(lambda)
		st_matrix("abs_eig", abs_eig)
	`end_mata'
	
	mat list abs_eig
	scalar effective_number_of_tests = 0
	local rows = rowsof(abs_eig)
	local cols = colsof(abs_eig)
	forval i = 1/`rows' {
		local h = 0
		forval j = 1/`cols' {
			if abs_eig[`i', `j'] >= 1 {
				local h = 1 + (abs_eig[`i', `j'] - floor(abs_eig[`i', `j']))
			} 
			else {
				local h = abs_eig[`i', `j'] - floor(abs_eig[`i', `j'])
			}
			scalar effective_number_of_tests = effective_number_of_tests + `h'
		}
	}
	di "Effective number of tests: " effective_number_of_tests
	egen `min_val' = rowmin(`varlist')
	gen _li_ji_p = min(1, `min_val' * effective_number_of_tests)
	di "Li_Ji: " _li_ji_p
}

else if "`method'" == "gaop" {
	local end_mata end 

	mkmat `corr_matrix', matrix(correlation_matrix)
	matrix eigenvalues lambda V = correlation_matrix
	mat list lambda
	mata: 
		lambda = st_matrix("lambda")
		sum_eig = rowsum(lambda)
		st_matrix("sum_eig", sum_eig)
	`end_mata'
	mat list sum_eig
	scalar total_sum = sum_eig[1,1]
	local rows = rowsof(lambda)
	local cols = colsof(lambda)
	local C = 0.995
	scalar effective_number_of_tests = 0
	local x_sum = 0
	forval x = 1/`cols' {
		local x_sum = `x_sum' + lambda[1, `x']
		if (`x_sum' / total_sum) > `C' {
			scalar effective_number_of_tests = `x'
		}
	}
	di "Effective number of tests: " effective_number_of_tests
	egen min_val = rowmin(`varlist')
	gen _gao_p = min(1, min_val * effective_number_of_tests)
	di "GAO: " _gao_p
}
else if "`method'" == "galweyp" {
	local end_mata end 

	mkmat `corr_matrix', matrix(correlation_matrix)
	matrix eigenvalues lambda V = correlation_matrix
	mat list lambda
	mata: 
		lambda = st_matrix("lambda")
		lambda_prime = J(rows(lambda), cols(lambda), 0)
		for (i=1; i<=rows(lambda); i++) {
			for (j=1; j<=cols(lambda); j++) {
				if (lambda[i, j] > 0) {
					lambda_prime[i, j] = lambda[i, j]
				}
			}
		}
		st_matrix("lambda_prime", lambda_prime)
		sum_lambda_prime = rowsum(lambda_prime)
		st_matrix("sum_lambda_prime", sum_lambda_prime)
		sum_eig = rowsum(sqrt(lambda_prime))^2
		st_matrix("sum_eig", sum_eig)
	`end_mata'
	mat list lambda_prime
	mat list sum_lambda_prime
	mat list sum_eig
	scalar sum_squared_lambda_prime = sum_eig[1,1]
	scalar sum_lambda_prime = sum_lambda_prime[1,1]
	scalar effective_number_of_tests = sum_squared_lambda_prime / sum_lambda_prime
	di "Effective number of tests: " effective_number_of_tests
	egen min_val = rowmin(`varlist')
	gen _galwey_p = min(1, min_val * effective_number_of_tests)
	di "Galwey: " _galwey_p
}*/
else if "`method'" == "EBM" {
	corr `varlist', cov
	matrix covar_matrix = r(C)
	local m = rowsof(covar_matrix)
	
	*matrix accum R = var2-var12, nocons dev
	*matrix R = corr(R)
	*mat list R

	scalar df_fisher = 2 * `m' 
	scalar Expected = 2 * `m' 

	local cov_sum = 0

	forval r = 1/`m' {
		forval c = `r'/`m' {
			local cov_sum = `cov_sum' + covar_matrix[`r', `c']
		}
	}
	di `cov_sum'
	scalar Var = 4 * `m' + 2 * `cov_sum'
	local c = Var / (2 * Expected)
	scalar df_brown = 2 * Expected^2 / Var
	di `c'
	 if df_brown > df_fisher {
		scalar df_brown = df_fisher
		local c = 1
	}

	foreach var of varlist `varlist'{
		* Calculate the t-value for each observation
		gen _blog_`var' = ln(`var')
	}

	egen `chi_squared_sum' = rowtotal(_blog_*)   

	

	gen _bfisher_p = 1 - chi2(df_fisher, -2 * `chi_squared_sum') 
	di "Fisher p-value:", _bfisher_p
	
	gen chi_squared_sum_brown = `chi_squared_sum' / `c'
	gen _p_brown = 1 - chi2(df_brown, -2 * chi_squared_sum_brown)
	di "Brown p-value:", _p_brown
	
}	
else if "`method'" == "KostsMethodEBM" {
	scalar a1 = 3.263
	scalar a2 = 0.710
	scalar a3 = 0.027
	corr `varlist'
	matrix corr_matrix = r(C)
	local m = rowsof(corr_matrix)
	
	*matrix accum R = _all, nocons dev
	*matrix R = corr(R)
	*mat list R

	local nrows = rowsof(corr_matrix)
	local ncols = colsof(corr_matrix)

	* Initialize covar_matrix with missing values
	matrix covar_matrix = J(`nrows', `ncols', .)

	* Loop to compute covariance
	forval i = 1/`nrows' {
		forval j = 1/`ncols' {
			local corr_ij = corr_matrix[`i', `j']
			
			* Calculate covariance (assuming Kost's polynomial fit)
			scalar covar = a1 * `corr_ij' + a2 * (`corr_ij'^ 2) + a3 * (`corr_ij'^ 3)
			
			* Assign covariance values to the matrix
			matrix covar_matrix[`i', `j'] = covar
			matrix covar_matrix[`j', `i'] = covar
		}
	}
	
	mat list covar_matrix
	local m = rowsof(covar_matrix)
	scalar df_fisher = 2 * `m' 
	scalar Expected = 2 * `m' 

	local cov_sum = 0

	forval r = 1/`m' {
		forval c = `r'/`m' {
			local cov_sum = `cov_sum' + covar_matrix[`r', `c']
		}
	}
	di `cov_sum'
	scalar Var = 4 * `m' + 2 * `cov_sum'
	local c = Var / (2 * Expected)
	scalar df_brown = 2 * Expected^2 / Var
	di `c'
	 if df_brown > df_fisher {
		scalar df_brown = df_fisher
		local c = 1
	}

	foreach var of varlist `varlist'{
		* Calculate the t-value for each observation
		gen _klog_`var' = ln(`var')
	}

	egen `chi_squared_sum' = rowtotal(_klog_*)   

	* Calculate combined p-value
	gen _kfisher_p = 1 - chi2(df_fisher, -2 * `chi_squared_sum') 
	di "Fisher p-value:", _kfisher_p
	
	gen chi_squared_sum_kost_brown = `chi_squared_sum' / `c'
	gen _p_kost_brown = 1 - chi2(df_brown, -2 * chi_squared_sum_kost_brown)
	di "Brown p-value:", _p_kost_brown
}
else if "`method'" == "BrownbyYang" {
	scalar c1 = 3.9081 
	scalar c2 = 0.0313
	scalar c3 = 0.1022
	scalar c4 = -0.1378
	scalar c5 = 0.0941
	
	corr `varlist'
	matrix corr_matrix = r(C)
	local m = rowsof(corr_matrix)
	
	*matrix accum R = _all, nocons dev
	*matrix R = corr(R)
	*mat list R

	local nrows = rowsof(corr_matrix)
	local ncols = colsof(corr_matrix)

	* Initialize covar_matrix with missing values
	matrix delta_matrix = J(`nrows', `ncols', .)

	* Loop to compute covariance
	forval i = 1/`nrows' {
		forval j = 1/`ncols' {
			local corr_ij = corr_matrix[`i', `j']
			local biased_corrected_cor = `corr_ij' *(1 + (1 - `corr_ij'^2 / 2 * (`ncols' -3)))

			scalar f_r = c1 * (`biased_corrected_cor'  ^ 2) + c2 * (`biased_corrected_cor' ^ 4) + c3 * (`biased_corrected_cor' ^ 6 ) + c4 * (`biased_corrected_cor' ^ 8) + c5 * (`biased_corrected_cor' ^ 10)

			
			local bias = (c1 / `ncols')*(1 - `biased_corrected_cor'^2)^2

			scalar unbiased_delta = f_r - `bias'
			matrix delta_matrix[`i', `j'] = unbiased_delta
			matrix delta_matrix[`j', `i'] = unbiased_delta
		}
	}

	mat list delta_matrix
	local m = rowsof(delta_matrix)
	scalar df_fisher = 2 * `m' 
	scalar Expected = 2 * `m' 

	scalar mean = 2 * `m' 

	local delta_sum = 0

	forval r = 1/`m' {
		forval c = `r'/`m' {
			local delta_sum = `delta_sum' + delta_matrix[`r', `c']
		}
	}
	di `delta_sum'
	scalar Var = 4 * `m' + 2 * `delta_sum'
	local v = 2*(mean^2 / Var)
	scalar gamma_value = Var / mean
	scalar scale = `v' / 2
	di `v'
	 

	foreach var of varlist `varlist'{
		* Calculate the t-value for each observation
	   gen _ylog_`var' = ln(`var')  
	}

	egen `chi_squared_sum' = rowtotal(_ylog_*)   
	
	* Calculate combined p-value
	gen  _yfisher_p = 1 - chi2(df_fisher, -2 * `chi_squared_sum')
	di "Fisher p-value:", _yfisher_p
	
	gen _p_yang = 1 - gammap( 2 * gamma_value, -2 * `chi_squared_sum')
	di "Brown p-value by Yang:", _p_yang
}	

	

else {
        di "Invalid Method"
    }
end
