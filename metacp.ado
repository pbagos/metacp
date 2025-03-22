capture program drop metacp
program define metacp
    version 13.0
    syntax varlist(numeric),method(string) stat(string) [ corr_matrix(varlist) weight_matrix(varlist) ] 
    
	
    tempvar p1 p2 chi1_squared_sum chi2_squared_sum combined_p1 combined_p2 num_vars t_sum z_mean z_value chi_squared_sum z_sum invchi_squared_sum T_sum t min_val C_sum c h C chi_squared_sum_kost_brown chi_squared_sum_brown
    tempvar num_vars1 num_vars2 t1_sum t2_sum chi1_squared_sum chi2_squared_sum z1_mean z2_mean z1_value z2_value invchi1_squared_sum invchi2_squared_sum T1_sum T2_sum t1 t2 C1_sum C2_sum c1 c2
	tempvar sum_weight min_val1 min_val2 effective_number_of_tests1 effective_number_of_tests2 combinedf_p1 combinedf_p2 chi1_squared_sum_brown chi2_squared_sum_brown combinedb_p1 combinedb_p2 chi1_squared_sum_kost_brown chi2_squared_sum_kost_brown combinedk_p1 combinedk_p2 combinedy_p1 combinedy_p2
	tempvar z_sum 
	egen `num_vars' = rownonmiss(`varlist')
	
if "`method'" == "stouffer" {
	 * Define the numerical value of pi
	local pi_value = 3.141592653589793

	foreach var of varlist `varlist' {
		if "`stat'" == "pvalue" {
			gen z_score_`var' = invnorm(1 - `var')
		}
		else{
			gen z_score_`var' = `var'
		}
	}
	
	egen `z_sum' = rowtotal(z_score_*)    

	* Calculate combined p-value
	gen _stouffer_p = 1 - normal(`z_sum' / sqrt(`num_vars')) 
	di "Combined p-value:", _stouffer_p
}
else if "`method'" == "wstouffer" {
	local weight_matrix = "`weight_matrix'"
	egen `sum_weight' = rowtotal(`weight_matrix')

		
	// Count number of weight variables
	local k: word count `varlist'

	// Initialize sum variables
	gen _z_w_sum = 0
	gen _weight2_sum = 0

	// Convert p-values to z-scores if needed
	foreach var of varlist `varlist' {
		if "`stat'" == "pvalue" {
			gen zi_`var' = invnorm(1 - `var')
		}
		else {
			gen zi_`var' = `var'
		}
	}

	// Compute weighted sum and sum of squared weights
	gen _combined_z = .

	forvalues i = 1/`k' {
		local weight: word `i' of `weight_matrix'  // Extract i-th weight
		di `i'
		local var: word `i' of `varlist'  // Corresponding z-score variable
		di `var'
		replace _z_w_sum = _z_w_sum + zi_`var' * `weight'
		replace _weight2_sum = _weight2_sum + (`weight' ^ 2)
	}

	// Compute final combined z-score
	replace _combined_z = _z_w_sum / sqrt(_weight2_sum)

	// Compute the p-value
	gen _wstouffer_p = 1 - normal(_combined_z)   
        
}

	
else if "`method'" == "corstouffer" {

	corr `varlist', cov	
	matrix corr_matrix = r(C)

	local k = colsof(corr_matrix)  // Number of variables
	local sum_r = 0  // Initialize sum accumulator

	forval i = 1/`=`k' - 1' {  
		forval j = 1/`=`i' - 1' {  
			local sum_r = `sum_r' + corr_matrix[`i', `j']
		}
	}

	local result = `k' + 2 * `sum_r'
	display "Result: " `result'

	foreach var of varlist `varlist' {

		if "`stat'" == "pvalue" {
			gen zi_`var' = invnorm(1 - `var')
		}
		else{
			gen zi_`var' = `var'

		}
	}
	egen `z_sum' = rowtotal(zi_*)    

	gen _corstouffer_p = 1 - normal(`z_sum' / sqrt(`result'))
	di "Combined p-value:", _corstouffer_p

}	
if "`method'" == "wcorstouffer" {
	local weight_matrix = "`weight_matrix'"
	egen `sum_weight' = rowtotal(`weight_matrix')

    // Compute the correlation matrix
    corr `varlist'
    matrix corr_matrix = r(C)
	// Count number of weight variables
	local k: word count `varlist'

	// Initialize sum variables
	gen _z_w_sum = 0
	gen _w_cor_sum = 0
  
    foreach var of varlist `varlist' {
        // Convert p-values to z-scores if needed
        if "`stat'" == "pvalue" {
            gen zi_`var' = invnorm(1 - `var')
        }
        else {
            gen zi_`var' = `var'
        }
	}
	
	
	forvalues i = 1/`k' {
		local weight: word `i' of `weight_matrix'  // Extract i-th weight
		di `i'
		local var: word `i' of `varlist'  // Corresponding z-score variable
		di `var'
		replace _z_w_sum = _z_w_sum + zi_`var' * `weight'
	}
	// Compute weighted sum and sum of squared weights
	gen _combined_z = .

	
	forvalues i = 1/`k' {
        local weight_i: word `i' of `weight_matrix'
        forvalues j = 1/`k' {
            local weight_j: word `j' of `weight_matrix'

            // Extract correlation coefficient
            matrix corr_val = corr_matrix[`i', `j']
            scalar corr_ij = corr_val[1,1]

            replace _w_cor_sum = _w_cor_sum + (`weight_i' * `weight_j' * corr_ij)
        }
    }

    // Compute final combined Z-score
    replace _combined_z = _z_w_sum / sqrt(_w_cor_sum)

    // Compute the p-value
    gen _wcorstouffer_p = 1 - normal(_combined_z)
}
	
if "`stat'" == "zscore" {
    foreach var of varlist `varlist' {
       capture gen p1_`var' = 1 - normal(`var')   // Right-tailed p-value
    }
	foreach var of varlist `varlist' {
       capture gen p2_`var' = normal(`var')       // Left-tailed p-value
    }
	
	
	

	if "`method'" == "logit" {
		local pi_value = 3.141592653589793
		local C = sqrt(`num_vars' * `pi_value'^2 * (5 * `num_vars' + 2) / (3 * (5 * `num_vars' + 4)))

		if "`stat'" == "zscore" {

			foreach var of varlist `varlist' {
				* Calculate the t-value for each observation
				gen t1_value_`var' = ln(p1_`var' / (1 - p1_`var'))
			}
			foreach var of varlist `varlist' {
				* Calculate the t-value for each observation
				gen t2_value_`var' = ln(p2_`var' / (1 - p2_`var'))
			}
			* Sum the t-values for each observation
			egen `t1_sum' = rowtotal(t1_value_*)
			gen `combined_p1' = 2 * (1 - t(2 * `num_vars', abs(-`t1_sum' /`C')))

			egen `t2_sum' = rowtotal(t2_value_*)
			gen `combined_p2' = 2 * (1 - t(2 * `num_vars', abs(-`t2_sum' /`C')))

			gen _logitp = 2 * min(`combined_p1', `combined_p2')
			* Display the combined p-value
			di "Combined p-value: " _logitp
		  }

		else {
			* Calculate the constant C
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
	}
	else if "`method'" == "meanp" {	
		local pi_value = 3.141592653589793
		if "`stat'" == "zscore" {

			egen `z1_mean' = rowmean(p1_*)
			egen `z2_mean' = rowmean(p2_*)
			
			gen `z1_value' = (0.5 - `z1_mean')* sqrt(12 * `num_vars')
			gen `z2_value' = (0.5 - `z2_mean')* sqrt(12 * `num_vars')  

			* Calculate combined p-value
			gen `combined_p1' = 1 - normal(`z1_value') 
			gen `combined_p2' = 1 - normal(`z2_value')
			
			gen _meanp = 2 * min(`combined_p1', `combined_p2')

				
		}
		else{
			egen `z_mean' = rowmean(`varlist')   
			gen `z_value' = (0.5 - `z_mean')* sqrt(12 * `num_vars')  

			* Calculate combined p-value
			gen _meanp = 1 - normal(`z_value') 
		}
		
		di "Combined p-value:", _meanp
	}
	else if "`method'" == "fisher" {
		if "`stat'" == "zscore" {
			foreach var of varlist `varlist' {
				gen log1_`var' = ln(p1_`var')
			}
			foreach var of varlist `varlist' {
				gen log2_`var' = ln(p2_`var')
			}
			egen `chi1_squared_sum' = rowtotal(log1_*)   
			egen `chi2_squared_sum' = rowtotal(log2_*)
			
			gen `combined_p1' = 1 - chi2(2 * `num_vars', -2 * `chi1_squared_sum')
			gen `combined_p2' = 1 - chi2(2 * `num_vars', -2 * `chi2_squared_sum')
			
			gen _fisher_p = 2 * min(`combined_p1', `combined_p2')

		}
		else{
			foreach var of varlist `varlist' {
				* Calculate the t-value for each observation
				gen log_`var' = ln(`var')
			}

			egen `chi_squared_sum' = rowtotal(log_*)   

			* Calculate combined p-value
			gen _fisher_p = 1 - chi2(2 * `num_vars', -2 * `chi_squared_sum')
		}
		di "Combined p-value:", _fisher_p
	}

	else if "`method'" == "lancaster" {
		local weight_matrix = "`weight_matrix'"
		egen `sum_weight' = rowtotal(`weight_matrix')
		if "`stat'" == "zscore" {
			local k: word count `varlist'

			gen invchi1_sum = 0
			gen invchi2_sum = 0
			
			forvalues i = 1/`k' {
				local weight: word `i' of `weight_matrix'  // Extract i-th weight
				di `i'
				local var: word `i' of `varlist'  // Corresponding z-score variable
				di `var'
				gen invchi1_`var' = invchi2(`weight', 1 - p1_`var')
				replace invchi1_sum = invchi1_sum + invchi1_`var'
			}
			forvalues i = 1/`k' {
				local weight: word `i' of `weight_matrix'  // Extract i-th weight
				di `i'
				local var: word `i' of `varlist'  // Corresponding z-score variable
				di `var'
				gen invchi2_`var' = invchi2(`weight', 1 - p2_`var')
				replace invchi2_sum = invchi2_sum + invchi2_`var'
		
			
			}
			
			
			// Calculate the combined p-values based on the summed inverse chi-square values
			gen `combined_p1' = 1 - chi2(`sum_weight', invchi1_sum) 
			gen `combined_p2' = 1 - chi2(`sum_weight', invchi2_sum)

			// Lancaster combined p-value: take the minimum of both p-values
			gen _lancaster_p = 2 * min(`combined_p1', `combined_p2')
		}
		else {
			// Loop through each variable in the varlist to calculate the inverse chi-squared statistic
			forvalues i = 1/`k' {
				local weight: word `i' of `weight_matrix'  // Extract i-th weight
				di `i'
				local var: word `i' of `varlist'  // Corresponding z-score variable
				di `var'
				replace invchi_squared_`var' = invchi2(`weight', 1 -`var')			
			}

			// Sum up the inverse chi-square statistics across variables
			egen `invchi_squared_sum' = rowtotal(invchi_squared_*)   

			// Calculate the total weight (sum of all weights in the weight matrix)
			egen `sum_weight' = rowtotal(`weight_matrix')

			// Calculate the combined p-value based on the summed inverse chi-square
			gen _lancaster_p = 1 - chi2(`sum_weight', `invchi_squared_sum')
		}
	}




	else if "`method'" == "invchi" {
		if "`stat'" == "zscore" {
			foreach var of varlist `varlist' {
				gen invchi1_squared_`var' = invchi2(1, 1 - p1_`var')
			}
			foreach var of varlist `varlist' {
				gen invchi2_squared_`var' = invchi2(1, 1 - p2_`var')
			}
			egen `invchi1_squared_sum' = rowtotal(invchi1_squared_*)
			egen `invchi2_squared_sum' = rowtotal(invchi2_squared_*)
			
			gen `combined_p1' = 1 - chi2(`num_vars', `invchi1_squared_sum') 
			gen `combined_p2' = 1 - chi2(`num_vars', `invchi2_squared_sum')
			
			gen _inv_chi2_p = 2 * min(`combined_p1', `combined_p2')


		}
		else{
		* Loop through each variable in the varlist
			foreach var of varlist `varlist' {
				* Calculate the inverse chi-squared statistic for each observation
				gen invchi_squared_`var' = invchi2(1, 1 - `var')
			}

			* Calculate the sum of inverse chi-squared statistics across variables
			egen `invchi_squared_sum' = rowtotal(invchi_squared_*)   

			* Calculate combined p-value
			gen _inv_chi2_p = 1 - chi2(`num_vars', `invchi_squared_sum') 
		}

		* Display combined p-value
		di "Combined p-value:", _inv_chi2_p

	}

	else if "`method'" == "binomial" {
    local alpha 0.05

		if "`stat'" == "zscore" {
			gen combined_p1 = .
			gen combined_p2 = .
			gen _binomial_p = .

			forval i = 1/`=_N' {
				// Reset counts for each observation
				local num_of_successes1 = 0
				local num_of_successes2 = 0
				local num_vars1 = 0
				local num_vars2 = 0

				// Loop through each variable in the varlist
				foreach var of varlist `varlist' {
					if p1_`var'[`i'] < `alpha' {
						local num_of_successes1 = `num_of_successes1' + 1
					}
					local num_vars1 = `num_vars1' + 1

					if p2_`var'[`i'] < `alpha' {
						local num_of_successes2 = `num_of_successes2' + 1
					}
					local num_vars2 = `num_vars2' + 1
				}

				// Compute binomial test p-values
				local combined_p1 = 2 * binomialp(`num_vars1', `num_of_successes1', 0.05)
				local combined_p2 = 2 * binomialp(`num_vars2', `num_of_successes2', 0.05)

				// Store results in dataset
				replace combined_p1 = `combined_p1' in `i'
				replace combined_p2 = `combined_p2' in `i'
				replace _binomial_p = min(`combined_p1', `combined_p2') in `i'

				// Debugging output
				di "Combined p-value for observation `i': " _binomial_p[`i']
			}
		} 
		else {
			gen _binomial_p = .

			forval i = 1/`=_N' {
				local num_of_successes = 0
				local num_vars = 0

				foreach var of varlist `varlist' {
					if `var'[`i'] < `alpha' {
						local num_of_successes = `num_of_successes' + 1
					}
					local num_vars = `num_vars' + 1
				}

				// Compute binomial test p-value
				local combined_p = 2 * binomialp(`num_vars', `num_of_successes', 0.05)

				// Store the result in the dataset
				replace _binomial_p = `combined_p' in `i'

				// Debugging output
				di "Combined p-value for observation `i': " _binomial_p[`i']
			}
		}
	
	}
	else if "`method'" == "cct" {
		
		local pi_val = 3.141592653589793
		local gamma_val = 1
		if "`stat'" == "zscore" {
			foreach var of varlist `varlist'{
				gen T1_`var' = tan((0.5 - p1_`var') * `pi_val')
			}
			foreach var of varlist `varlist'{
				gen T2_`var' = tan((0.5 - p2_`var') * `pi_val')
			}
			egen `T1_sum' = rowtotal(T1_*)
			egen `T2_sum' = rowtotal(T2_*)
			
			gen `t1' = `T1_sum' / `num_vars'
			gen `t2' = `T2_sum' / `num_vars'
			
			gen `combined_p1' = 1 - (0.5 + atan(`t1' / `gamma_val') / `pi_val')
			gen `combined_p2' = 1 - (0.5 + atan(`t2' / `gamma_val') / `pi_val')

			gen _cauchy_p = 2 * min(`combined_p1', `combined_p2')
				
		}
		else{
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
		}
		* Display combined p-value
		di "Combined p-value:", _cauchy_p
		
	}	

	else if "`method'" == "minp" {
		if "`stat'" == "zscore" {
			egen `min_val1' = rowmin(p1_*)  // Minimum of p1 for each observation
			egen `min_val2' = rowmin(p2_*)  // Minimum of p2 for each observation

			* Calculate the combined p-values
			gen `combined_p1' = 1 - ((1 - `min_val1') ^ `num_vars')
			gen `combined_p2' = 1 - ((1 - `min_val2') ^ `num_vars')

			* Calculate the final minimum p-value
			gen _min_p = 2 * min(`combined_p1', `combined_p2')
		}
		else{
		
		* Find the minimum value across variables in each observation
		egen `min_val' = rowmin(`varlist')
		
		* Calculate combined p-value
		gen _min_p = 1 - (1 - `min_val') ^ `num_vars'
		}
		* Display combined p-value
		di "Combined p-value:", _min_p
	}

	else if "`method'" == "mcm" {
		
		gen _MCM_p = 2 * min(_cauchy_p, _min_p, 0.5)
		* Display combined p-value
		di "Combined p-value by MCM(MinP-Cauchy-MinP):", _MCM_p
	}	
	*For the CMC method the user should first calculate the '_cauchy_p' and the '_min_p'.
	*The CMC method should be run like this 'metacp __cauchy_p _min_p,method(CMC)'
	else if "`method'" == "cmc" {
		local pi_val = 3.141592653589793
		local gamma_val = 1

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

	else if "`method'" == "bonferroni" {
		if "`stat'" == "zscore" {
				
			egen `min_val1' = rowmin(p1_*)
			egen `min_val2' = rowmin(p2_*)

			gen `combined_p1' = min(1, `min_val1' * `num_vars')
			gen `combined_p2' = min(1, `min_val2' * `num_vars')

			gen _bonferroni_p = 2 * min(`combined_p1', `combined_p2')
		}
		else{
			egen `min_val' = rowmin(`varlist')
		
			* Calculate combined p-value
			gen _bonferroni_p = min(1, `min_val' * `num_vars')
		}
		* Display combined p-value
		di "Combined p-value:", _bonferroni_p
	}

	else if "`method'" == "adbonferroni" {
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

	   
		if "`stat'" == "zscore" {
			
			egen `min_val1' = rowmin(p1_*)
			egen `min_val2' = rowmin(p2_*)

			gen `combined_p1' = min(1, `min_val1' *  g_adjusted)
			gen `combined_p2' = min(1, `min_val2' *  g_adjusted)

			gen _bonferronig_p = 2 * min(`combined_p1', `combined_p2')
		}
		else{
			egen `min_val' = rowmin(`varlist')
		
			* Calculate combined p-value using Bonferroni correction
			gen _bonferronig_p = min(1, `min_val' * g_adjusted)
		}
		
		* Display combined p-value
		di "Combined p-value:", _bonferronig_p
	}

	else if "`method'" == "cn" {
		local end_mata end 
		mkmat `corr_matrix', matrix(correlation_matrix)
		matrix eigenvalues lambda V = correlation_matrix
		mat list lambda
		mata: lambdaVar(" lambda ")
		mat eig_variance = variance

		* Extract the value of the standard deviation
		scalar sample_variance = eig_variance[1,1]
		di "Sample variance of Eigenvalues: " sample_variance
		
		gen `effective_number_of_tests1' = 1 + (`num_vars' - 1) * (1 - sample_variance / `num_vars')

		
		if "`stat'" == "zscore" {

			egen `min_val1' = rowmin(p1_*)
			egen `min_val2' = rowmin(p2_*)

			gen `combined_p1' = min(1, `min_val1' *  effective_number_of_tests)
			gen `combined_p2' = min(1, `min_val2' *  effective_number_of_tests)

			gen _cheverud_nyholt_p = 2 * min(`combined_p1', `combined_p2')
		}
		else{

			egen `min_val' = rowmin(`varlist')
			gen _cheverud_nyholt_p = min(1, `min_val' * effective_number_of_tests)
		}
		di "CN: " _cheverud_nyholt_p
	}

	else if "`method'" == "lj" {
		local end_mata end 
		mkmat `corr_matrix', matrix(correlation_matrix)
		matrix eigenvalues lambda V = correlation_matrix
		mat list lambda
		mata: absEig(" lambda ")
		
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
		if "`stat'" == "zscore" {
			egen `min_val1' = rowmin(p1_*)
			egen `min_val2' = rowmin(p2_*)

			gen `combined_p1' = min(1, `min_val1' *  effective_number_of_tests)
			gen `combined_p2' = min(1, `min_val2' *  effective_number_of_tests)

			gen _li_ji_p = 2 * min(`combined_p1', `combined_p2')
			}
		else{
			 
			egen `min_val' = rowmin(`varlist')
			gen _li_ji_p = min(1, `min_val' * effective_number_of_tests)
		}
			
		di "Li_Ji: " _li_ji_p
	}

	else if "`method'" == "gao" {

		mkmat `corr_matrix', matrix(correlation_matrix)
		matrix eigenvalues lambda V = correlation_matrix
		mat list lambda
		mata: sumEig(" lambda ")

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
		if "`stat'" == "zscore" {
			egen `min_val1' = rowmin(p1_*)
			egen `min_val2' = rowmin(p2_*)

			gen `combined_p1' = min(1, `min_val1' *  effective_number_of_tests)
			gen `combined_p2' = min(1, `min_val2' *  effective_number_of_tests)

			gen _gao_p = 2 * min(`combined_p1', `combined_p2')
			}
		else{
			 
			egen `min_val' = rowmin(`varlist')
			gen _gao_p = min(1, `min_val' * effective_number_of_tests)
		}
		
		di "GAO: " _gao_p
	}
	else if "`method'" == "galwey" {
		 

		mkmat `corr_matrix', matrix(correlation_matrix)
		matrix eigenvalues lambda V = correlation_matrix
		mat list lambda
		mata: lambdaPrime(" lambda ")
			
		mat list lambda_prime
		mat list sum_lambda_prime
		mat list sum_eig
		scalar sum_squared_lambda_prime = sum_eig[1,1]
		scalar sum_lambda_prime = sum_lambda_prime[1,1]
		scalar effective_number_of_tests = sum_squared_lambda_prime / sum_lambda_prime
		di "Effective number of tests: " effective_number_of_tests
		
		if "`stat'" == "zscore" {
			egen `min_val1' = rowmin(p1_*)
			egen `min_val2' = rowmin(p2_*)
			
			gen `combined_p1' = min(1, `min_val1' *  effective_number_of_tests)
			gen `combined_p2' = min(1, `min_val2' *  effective_number_of_tests)

			gen _galwey_p = 2 * min(`combined_p1', `combined_p2')
		}
		else{
			 
			egen `min_val' = rowmin(`varlist')
			gen _galwey_p = min(1, `min_val' * effective_number_of_tests)
		}
			
		di "Galwey: " _galwey_p
	}
	else if "`method'" == "ebm" {

		if "`stat'" == "zscore" {
		corr p1_*, cov
		matrix covar_matrix1 = r(C)
		corr p2_*, cov
		matrix covar_matrix2 = r(C)

		local m = rowsof(covar_matrix1)
		
		*matrix accum R = var2-var12, nocons dev
		*matrix R = corr(R)
		*mat list R

		scalar df_fisher = 2 * `m' 
		scalar Expected = 2 * `m' 

		local cov_sum1 = 0
		local cov_sum2 = 0


		forval r = 1/`m' {
			forval c = `=`r' + 1'/`m' {
				local cov_sum1 = `cov_sum1' + covar_matrix1[`r', `c']
				local cov_sum2 = `cov_sum2' + covar_matrix2[`r', `c']

			}
		}
		di `cov_sum1'
		di `cov_sum2'

		scalar Var1 = 4 * `m' + 2 * `cov_sum1'
		scalar Var2 = 4 * `m' + 2 * `cov_sum2'

		local c1 = Var1 / (2 * Expected)
		local c2 = Var2 / (2 * Expected)

		scalar df_brown1 = 2 * Expected^2 / Var1
		scalar df_brown2 = 2 * Expected^2 / Var2

		di `c1'
		di `c2'

		 if df_brown1 > df_fisher {
			scalar df_brown1 = df_fisher
			local c1 = 1

		}
		if df_brown2 > df_fisher {
			scalar df_brown2 = df_fisher
			local c2 = 1

		}

			foreach var of varlist `varlist'{	
				gen _blog1_`var' = ln(p1_`var')
			}
			foreach var of varlist `varlist'{	
				gen _blog2_`var' = ln(p2_`var')
			}

			egen `chi1_squared_sum' = rowtotal(_blog1_*) 
			egen `chi2_squared_sum' = rowtotal(_blog2_*) 
			gen `combinedf_p1' = 1 - chi2(df_fisher, -2 * `chi1_squared_sum') 
			gen `combinedf_p2' = 1 - chi2(df_fisher, -2 * `chi2_squared_sum')
			gen _bfisher_p =  2 * min(`combinedf_p1', `combinedf_p2')
			di "Fisher p-value:", _bfisher_p
			
			gen `chi1_squared_sum_brown' = `chi1_squared_sum' / `c1'
			gen `chi2_squared_sum_brown' = `chi2_squared_sum' / `c2'
			
			gen `combinedb_p1' = 1 - chi2(df_brown1, -2 * `chi1_squared_sum_brown') 
			gen `combinedb_p2' = 1 - chi2(df_brown2, -2 * `chi2_squared_sum_brown')
			gen _p_brown = 2 * min(`combinedb_p1', `combinedb_p2')
			di "Brown p-value:", _p_brown

		}
		else{
			corr `varlist', cov
			matrix covar_matrix = r(C)

			local m = rowsof(covar_matrix1)
			
			*matrix accum R = var2-var12, nocons dev
			*matrix R = corr(R)
			*mat list R

			scalar df_fisher = 2 * `m' 
			scalar Expected = 2 * `m' 

			local cov_sum1 = 0
			local cov_sum2 = 0


			forval r = 1/`m' {
				forval c = `=`r' + 1'/`m' {
					local cov_sum = `cov_sum' + covar_matrix1[`r', `c']
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
				
				gen _blog_`var' = ln(`var')
			}

			egen `chi_squared_sum' = rowtotal(_blog_*) 
		

			gen _bfisher_p = 1 - chi2(df_fisher, -2 * `chi_squared_sum') 
			di "Fisher p-value:", _bfisher_p
			
			gen chi_squared_sum_brown = `chi_squared_sum' / `c'
			gen _p_brown = 1 - chi2(df_brown, -2 * chi_squared_sum_brown)
			di "Brown p-value:", _p_brown
		}
		
	}	
	else if "`method'" == "kost" {
		scalar a1 = 3.263
		scalar a2 = 0.710
		scalar a3 = 0.027
		
		if "`stat'" == "zscore" {
			corr p1_*
			matrix corr_matrix1 = r(C)
			corr p2_*
			matrix corr_matrix2 = r(C)

			local nrows = rowsof(corr_matrix1)
			local ncols = colsof(corr_matrix1)

			* Initialize covar_matrix with missing values
			matrix covar_matrix = J(`nrows', `ncols', .)

			* Loop to compute covariance
			forval i = 1/`ncols' {
				forval j = `=`i' + 1'/`ncols' {
					local corr_ij1 = corr_matrix1[`i', `j']
					local corr_ij2 = corr_matrix1[`i', `j']

					
					* Calculate covariance (assuming Kost's polynomial fit)
					scalar covar1 = a1 * `corr_ij1' + a2 * (`corr_ij1'^ 2) + a3 * (`corr_ij1'^ 3)
					scalar covar2 = a1 * `corr_ij2' + a2 * (`corr_ij2'^ 2) + a3 * (`corr_ij2'^ 3)

					* Assign covariance values to the matrix
					matrix covar_matrix1[`i', `j'] = covar1
					matrix covar_matrix1[`j', `i'] = covar1
					
					matrix covar_matrix2[`i', `j'] = covar2
					matrix covar_matrix2[`j', `i'] = covar2
				}
			}
			
			mat list covar_matrix1
			mat list covar_matrix2

			local m = colsof(covar_matrix1)
			scalar df_fisher = 2 * `m' 
			scalar Expected = 2 * `m' 

			local cov_sum1 = 0
			local cov_sum2 = 0


			forval r = 1/`m' {
				forval c = `=`r' + 1'/`m' {
					local cov_sum1 = `cov_sum1' + covar_matrix1[`r', `c']
					local cov_sum2 = `cov_sum2' + covar_matrix2[`r', `c']

				}
			}
			di `cov_sum1'
			di `cov_sum2'

			scalar Var1 = 4 * `m' + 2 * `cov_sum1'
			scalar Var2 = 4 * `m' + 2 * `cov_sum2'

			local c1 = Var1 / (2 * Expected)
			local c2 = Var2 / (2 * Expected)

			scalar df_brown1 = 2 * Expected^2 / Var1
			scalar df_brown2 = 2 * Expected^2 / Var2

			di `c1'
			di `c2'
			 if df_brown1 > df_fisher {
				scalar df_brown1 = df_fisher
				local c1 = 1
			}
			 if df_brown2 > df_fisher {
				scalar df_brown2 = df_fisher
				local c2 = 1
			}
			foreach var of varlist `varlist'{
				gen _klog1_`var' = ln(p1_`var')
			}
			foreach var of varlist `varlist'{
				gen _klog2_`var' = ln(p2_`var')
			}

			egen `chi1_squared_sum' = rowtotal(_klog1_*) 
			egen `chi2_squared_sum' = rowtotal(_klog2_*) 
			gen `combinedf_p1' = 1 - chi2(df_fisher, -2 * `chi1_squared_sum') 
			gen `combinedf_p2' = 1 - chi2(df_fisher, -2 * `chi2_squared_sum')
			gen _kfisher_p =  2 * min(`combinedf_p1', `combinedf_p2')
			di "Fisher p-value:", _bfisher_p
			
			gen `chi1_squared_sum_kost_brown' = `chi1_squared_sum' / `c1'
			gen `chi2_squared_sum_kost_brown' = `chi2_squared_sum' / `c2'
			
			gen `combinedk_p1' = 1 - chi2(df_brown1, -2 * `chi1_squared_sum_kost_brown') 
			gen `combinedk_p2' = 1 - chi2(df_brown2, -2 * `chi2_squared_sum_kost_brown')
			gen _p_kost_brown = 2 * min(`combinedk_p1', `combinedk_p2')
			di "Brown p-value:", _p_brown

		}
		else{
			corr `varlist'
			matrix corr_matrix = r(C)	

			local nrows = rowsof(corr_matrix)
			local ncols = colsof(corr_matrix)

			* Initialize covar_matrix with missing values
			matrix covar_matrix = J(`nrows', `ncols', .)

			* Loop to compute covariance
			forval i = 1/`ncols' {
				forval j = `=`i' + 1'/`ncols' {
					local corr_ij = corr_matrix[`i', `j']
					
					* Calculate covariance (assuming Kost's polynomial fit)
					scalar covar = a1 * `corr_ij' + a2 * (`corr_ij'^ 2) + a3 * (`corr_ij'^ 3)
					
					* Assign covariance values to the matrix
					matrix covar_matrix[`i', `j'] = covar
					matrix covar_matrix[`j', `i'] = covar
				}
			}
			
			mat list covar_matrix
			local m = colsof(covar_matrix)
			scalar df_fisher = 2 * `m' 
			scalar Expected = 2 * `m' 

			local cov_sum = 0

			forval r = 1/`m' {
				forval c = `=`r' + 1'/`m' {
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
			
			gen `chi_squared_sum_kost_brown' = `chi_squared_sum' / `c'
			gen _p_kost_brown = 1 - chi2(df_brown, -2 * chi_squared_sum_kost_brown)
			di "Brown p-value:", _p_kost_brown
		}
	}
	else if "`method'" == "yang" {
		scalar c1 = 3.9081 
		scalar c2 = 0.0313
		scalar c3 = 0.1022
		scalar c4 = -0.1378
		scalar c5 = 0.0941
		
		if "`stat'" == "zscore" {
			foreach var of varlist `varlist' {
				gen p_`var' = 2 * (1 - normal(abs(`var')))
			}
		
			* Ensure p_* variables exist before running correlation
			ds p_*  // Lists all variables starting with "p_"
			if "`r(varlist)'" != "" {
				corr `r(varlist)'  // Runs correlation on all p_* variables
			} 
			else {
				display "No variables matching p_* exist."
			}
		}
		else{
			corr `varlist'
		}
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
		forval i = 1/`ncols' {
			forval j = `=`i' + 1'/`ncols' {
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
		local m = colsof(delta_matrix)
		scalar df_fisher = 2 * `m' 
		scalar Expected = 2 * `m' 

		scalar mean = 2 * `m' 

		local delta_sum = 0

		forval r = 1/`m' {
			forval c = `=`r' + 1'/`m' {
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
			gen _ylog_`var' = ln(p_`var')
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
	}

end

version 13
mata:
void lambdaVar(matrix lambda)
{
	lambda = st_matrix("lambda")

	st_matrix("variance", diagonal(variance(lambda')))
}
end

mata: 
void absEig(matrix lambda)
{
		lambda = st_matrix("lambda")
		abs_eig = abs(lambda)
		st_matrix("abs_eig", abs_eig)
}
end	

mata: 
void sumEig(matrix lambda)
{
		lambda = st_matrix("lambda")
		sum_eig = rowsum(lambda)
		st_matrix("sum_eig", sum_eig)
}
end	

mata: 
void lambdaPrime(matrix lambda)
{
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
}
end	
