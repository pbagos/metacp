A versatile software package for combining dependent or independent p-values

metacp is designed to perform meta-analysis of the p-values of the variables(studies) specified by the user. The program takes as input a text file of n rows and k columns, plus the headers, which contain the p-values or z-scores to be combined. Different rows can correspond to different SNPs or genes from a GWAS, or probes/genes from transcriptomics studies. The different columns can correspond to different statistical tests, different traits, different SNPs within a gene and so on, depending on the setting. In case the same SNP/gene appears in multiple rows, a meta-analysis can be performed prior to combining the different p-values from different columns. Nevertheless, the applicability is not limited in these cases since the program can be used even with different types of data (from social sciences, education, economics and so on) provided that the user specifies correctly the variables to be combined.

The command the user should call to run the program in the terminal/command line should have the form:
python metacp.py <input file> <output file> <p-values/z-scores> <Yes/No> < correlation matrix file> < weight matrix file>  <method1> <method2> <method3> ...

where the 1st argument is the name of the input text file the user should provide, the 2rd argument is the name of output file, the 3rf argument (p-values/z-scores) is the form of the statistic given in the input, that can be p-values or z-scores, in the 4th  argument the user specifies if meta-analysis should be performed (Yes/No) and in the 5th argument the user can provide a text file with the correlation matrix (in case of dependent p-values). If this argument is omitted the program calculates by default a correlation matrix from the observations according to the specifics of the requested method. In the 6th argument the user can provide a text file with the weight matrix (for methods that allow weighting). In the following arguments (7th ...) the user specifies which statistical method(s), he/she wants to use for the combination of the p-values/z-scores in the input file.

The available methods are:

logit: This method computes the combined p-value using a logistic transformation of the individual p-values. It sums the logit-transformed p-values and then applies a t-distribution to calculate the combined p-value.

meanp: This method calculates the combined p-value by taking the mean of the individual p-values and transforming it using a normal distribution.

fisher: This method uses Fisher's method to combine p-values by summing the log-transformed p-values and applying a chi-squared distribution to the result.

stouffer: This method combines p-values using Stouffer's Z-score method. It converts each p-value to a Z-score, sums these Z-scores, and then transforms the sum using a normal distribution.

wstouffer: This method combines p-values using Stouffer's Z-score method with weights. It converts each p-value to a Z-score, sums these Z-scores, and then transforms the sum divided by the square root of the sum of the square of weights using a normal distribution.

invchi: This method combines p-values using the inverse chi-squared method, summing the inverse chi-squared statistics and transforming the sum using a chi-squared distribution.

binomial: This method uses a binomial test to combine p-values, counting the number of p-values below a specified threshold and using the binomial distribution to calculate the combined p-value.

cct: This method combines p-values using a Cauchy distribution. It transforms each p-value using the tangent function, sums these transformed values, and then transforms the sum using the inverse Cauchy cumulative distribution function.

minp: This method calculates the combined p-value by taking the minimum of the individual p-values and adjusting it for the number of tests.

mcm: This method combines p-values using the MinP-Cauchy-MinP (MCM) method. It first calculates the combined p-values using both the Cauchy and MinP methods and then takes the minimum of these values.

cmc: This method combines p-values using the Cauchy-MinP-Cauchy (CMC) method. It transforms the individual p-values using the tangent function, sums these values, and then transforms the sum using the inverse Cauchy cumulative distribution function.

bonferroni: This method applies the Bonferroni correction to combine p-values by adjusting the minimum p-value for the number of tests.

adbonferroni: This method uses a generalized Bonferroni correction, adjusting the minimum p-value for the number of tests and the maximum correlation among the tests.

cn: This method combines p-values using the Cheverud-Nyholt method, which adjusts for the effective number of tests based on eigenvalues of the correlation matrix.

lj: This method uses the Li-Ji method to combine p-values, adjusting for the effective number of tests by considering the eigenvalues of the correlation matrix.

gao: This method combines p-values using Gao's method, adjusting for the effective number of tests by considering a threshold cumulative eigenvalue sum.

galwey: This method combines p-values using Galwey's method, adjusting for the effective number of tests based on the ratio of the sum of squared eigenvalues to the sum of eigenvalues.

ebm: This method combines p-values using the Effective Brownian Motion (EBM) method, adjusting for correlations among tests by considering the covariance matrix.

kost: This method combines p-values using Kost's adjustment to the EBM method, incorporating a polynomial fit for the correlation structure.

yang: This method uses the Brown-Yang adjustment for combining p-values, accounting for correlations among tests with a polynomial fit and bias correction. 

corstouffer:  This method combines p-values using Stouffer's Z-score method for dependent tests. It converts each p-value to a Z-score, sums these Z-scores, and then transforms the sum, divided by the sum of variance of all z-scores using a normal distribution.

wcorstouffer:This method combines p-values using Weighted Stouffer's Z-score method for dependent tests. It converts each p-value to a Z-score, sums these Z-scores, and then transforms the sum divided by the sum of weights multiplied by the correlation coefficient using a normal distribution.


For example, if the command is:
python metacp.py SNP_p_values.txt output.txt p-values Yes correlation_matrix.txt weight_matrix.txt fisher meanp bonferroni wstouffer

the user has provided an input file named SNP_p_values.txt, which contains p-values to be combined, the user wants meta-analysis to be performed and the user provides a text file with the correlation matrix and a text file with the weight matrix for these p-values and wants to use meanp method, Fisher’s method, Bonferroni's method and Weighted Stouffer’s method to combine the p-values in the input file. The combined p-values will be printed in a text file provided by the user called output.txt.

Another example:
python metacp.py SNP_p_values.txt output.txt p-values No fisher meanp bonferroni

the user has provided an input file named SNP_p_values.txt, which contains p-values to be combined, the user does not want metanalysis to be performed and he does not provide a text file with the correlation matrix(the program will calculate by default a correlation matrix). meanp method, Fisher’s method and Bonferroni's to combine the p-values in the input file. The combined p-values will be printed in a text file provided by the user called output.txt.

