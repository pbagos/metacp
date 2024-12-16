# metacp
A versatile software package for combining dependent or independent p-values

metacp is designed to perform meta-analysis of the p-values of the variables(studies) specified by the user. The program takes as input a text file of n rows and k columns, plus the headers, which contain the p-values or z-scores to be combined. Different rows can correspond to different SNPs or genes from a GWAS, or probes/genes from transcriptomics studies. The different columns can correspond to different statistical tests, different traits, different SNPs within a gene and so on, depending on the setting. In case the same SNP/gene appears in multiple rows, a meta-analysis can be performed prior to combining the different p-values from different columns. Nevertheless, the applicability is not limited in these cases since the program can be used even with different 
types of data (from social sciences, education, economics and so on) provided that the user specifies correctly the variables to be combined.

The command the user should call to run the program in the terminal/command line(cmd) should have the form:

python metacp.py  input.txt output.txt p-values/z-scores  Yes/No correlation_matrix.txt/- method1 method2 method3 ...

where the 1st argument is the program in the .py(metacp.py) file, the 2nd argument(input.txt) is the input text file the user should provide, the 3rd argument(output.txt) is the ouput text file the user should provide where the results are printed, the 4th argument(p-values/z-scores) is the form the statistical tests have been measured, that can be p-values or z-scores, in the 5th argument the user specifies if metanalysis should be performed(Yes/No) and in the 6th argument the user can provide a text file with the correlation matrix(in case of dependent p-values). If this argument is ommitted the program calculates by default a correlation matrix. In the following arguments(7th...) the user specifies which method(s), he/she wants to use for the combination of the p-values/z-scores in the input file.


For example, if the command is:

python metacp.py SNP_p_values.txt output.txt p-values Yes correlation_matrix.txt fisher meanp bonferroni

the user has provided an input file named SNP_p_values.txt, which contains p-values to be combined, the user wants metanalysis to be performed and the user provides a text file with the correlation matrix for these p-values and wants to use meanp method, Fisher’s method and Bonferroni's to combine the p-values in the input file.
The combined p-values will be printed in a text file provided by the user called output.txt.

Another example:

python metacp.py SNP_p_values.txt output.txt p-values No fisher meanp bonferroni

the user has provided an input file named SNP_p_values.txt, which contains p-values to be combined, the user does not want metanalysis to be performed and he does not provide a text file with the correlation matrix(the program will calculate by default a correlation matrix). meanp method, Fisher’s method and Bonferroni's to combine the p-values in the input file.
The combined p-values will be printed in a text file provided by the user called output.txt.

