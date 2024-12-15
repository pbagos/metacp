# metacp
Combination of correlated p-values

In order to run the code in the metacp.py file you should have an input .txt file, like this provided as an example (SNP_p_values.txt), which consists of n rows and m columns, plus the headers, which contain the p-values or z-scores to be combined. 
Different rows can correspond to different SNPs or genes from a GWAS, or probes/genes from transcriptomics studies. The different columns can correspond to different statistical tests, different traits, different SNPs within a gene and so on, depending on the setting. 
In case the same SNP/gene appears in multiple rows, a meta-analysis can be performed prior to combining the different p-values from different columns.

The command the user should call to run the program in the terminal/command line(cmd) should have the form:

metacp.py  input.txt   p-values/z-scores  Yes/No  Yes/No

, where the 1st argument is the program in the .py(metacp.py) file, the 2nd argument(input.txt) is the input text file the user should provide, the 3rd argument(p-values/z-scores) is the form the statistical tests have been measured, that can be p-values or z-scores, in the 4th argument the user specifies if metanalysis should be performed(Yes/No) and in the 5th argument the user specifies if a correlation matrix(in case of dependent p-values) is provided or it should be calculated by the program(Yes/No).


For example, if the command is:

metacp.py  SNP_p_values.txt  p-values  Yes  Yes

the user has provided an input file named SNP_p_values.txt, which contains p-values to be combined, the user wants metanalysis to be performed and will provide a text file with the correlation matrix for these p-values.

After the user calls the program with the command above, the program will ask the user to enter the methods he/she wishes for the combination of the p-values/z-scores in the input file, separated by commas.

For example:   
 
logitp,meanp,fisher

The user wants to use logitp method, meanp method and Fisher’s method to combine the p-values/z-scores in the input file.

 The combined p-values will be printed in a text file called output.txt.
