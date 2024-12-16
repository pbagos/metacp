import math
import numpy as np
from scipy.stats import t
from scipy.stats import chi2
from scipy.stats import norm
from scipy.stats import gamma
from scipy.stats import pearsonr
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.special import chdtrc as chi2_cdf
from scipy.stats import moyal
import argparse
from collections import defaultdict

def logit(p_values):
    # Convert input to numpy array for numerical operations
    p_values = np.array(p_values)

    # Ensure that all probabilities are within the valid range (0, 1)
    p_values = np.clip(p_values, 1e-16, 1 - 1e-16)

    k = len(p_values)
    C = np.sqrt(k * np.pi ** 2 * (5 * k + 2) / (3 * (5 * k + 4)))
    # Compute the logit transform
    t_value = -np.sum(np.log((p_values) / (1 - p_values))) / C
    df = 2 * k
    combined_p = 2 * (1 - t.cdf(np.abs(t_value), df))

    return combined_p


def meanp(p_values):
    # Convert input to numpy array for numerical operations
    p_values = np.array(p_values)

    # Ensure that all p-values are within the valid range (0, 1)
    p_values = np.clip(p_values, 1e-09, 1 - 1e-09)
    k = len(p_values)
    # Compute the mean of the p-values
    z_value = (0.5 - np.mean(p_values)) * math.sqrt(12 * k)
    combined_p = 1 - norm.cdf(z_value)  # right-side test
    return combined_p


def fisher_method(p_values):
    # Convert input to numpy array for numerical operations
    p_values = np.array(p_values)

    # Ensure that all p-values are within the valid range (0, 1)
    p_values = np.clip(p_values, 1e-16, 1 - 1e-16)
    # Compute the combined test statistic using Fisher's method
    chi_squared = -2 * np.sum(np.log(p_values))
    # Calculate the degrees of freedom
    df = 2 * len(p_values)

    # Calculate the combined p-value using the chi-squared distribution
    fisher_p = 1 - chi2.cdf(chi_squared, df)

    return fisher_p


def stouffer(p_values, one_tailed=True):
    # Convert input to numpy array for numerical operations
    p_values = np.array(p_values)

    # Ensure that all p-values are within the valid range (0, 1)
    p_values = np.clip(p_values, 1e-16, 1 - 1e-16)
    if one_tailed:
        z_scores = norm.ppf(1 - p_values)  # norm.ppf =inverse cumulative distribution function of normal distribution
    else:
        z_scores = norm.ppf(1 - p_values / 2)

    # z_scores = norm.ppf(1 - p_values) # norm.ppf =inverse cumulative distribution function of normal distribution
    combined_z = np.sum(z_scores / math.sqrt(len(p_values)))
    combined_p = 1 - norm.cdf(combined_z)  # norm.cdf = cumulative distribution function of normal distribution

    return combined_p


def inverse_chi2(p_values):
    # Convert input to numpy array for numerical operations
    p_values = np.array(p_values)

    # Ensure that all p-values are within the valid range (0, 1)
    p_values = np.clip(p_values, 1e-16, 1 - 1e-16)

    df = len(p_values)

    chi_squared = np.sum(chi2.ppf((1 - p_values), 1))  # ppf = inverse of cdf
    inv_chi2_combined_p = 1 - chi2.cdf(chi_squared, df)  # chi2.cdf = cumulative distribution function of chi2

    return inv_chi2_combined_p


def binomial_test(p_values):
    p_values = np.array(p_values)

    k = len(p_values)
    alpha = 0.05
    # Count the number of significant p-values
    r = sum((p < alpha) for p in p_values)
    # Calculate the binomial probability of observing at most num_successes successes
    combined_p = 0
    for x in range(r, k + 1):
        combined_p += math.factorial(k) / (math.factorial(x) * math.factorial(k - x)) * (alpha ** x) * ((1 - alpha) ** (k - x))

    return combined_p


def cauchy_cdf(x, x0=0, gamma=1):
    """Cumulative distribution function (CDF) of the Cauchy distribution."""

    return 0.5 + (np.arctan((x - x0) / gamma) / np.pi)


def cauchy_method(p_values):
    # Convert input to numpy array for numerical operations
    p_values = np.array(p_values)

    # Ensure that all p-values are within the valid range (0, 1)
    p_values = np.clip(p_values, 1e-15, 1 - 1e-15)

    k = len(p_values)

    T = np.tan((0.5 - p_values) * np.pi)

    t = np.sum(T) / k
    # Calculate the combined p-value using the Cauchy distribution
    combined_p = 1 - cauchy_cdf(t)

    return combined_p


def minP(p_values):
    # Convert input to numpy array for numerical operations
    p_values = np.array(p_values)

    # Ensure that all p-values are within the valid range (0, 1)
    p_values = np.clip(p_values, 1e-15, 1 - 1e-15)

    k = len(p_values)

    # Calculate the combined p-value as the minimum p-value
    combined_p = 1 - (1 - np.min(p_values)) ** k

    return combined_p


def CMC(p_values):
    p_value_cauchy = cauchy_method(p_values)
    p_value_minp = minP(p_values)
    combined_p = cauchy_method((p_value_cauchy, p_value_minp))  # pCMC = CCT{pCCT , pMinP}
    return combined_p


def MCM(p_values):
    p_value_cauchy = cauchy_method(p_values)
    p_value_minp = minP(p_values)
    combined_p = 2 * min(p_value_cauchy, p_value_minp, 0.5)  # pMCM = 2 min{pCCT , pMinP, 0.5}

    return combined_p


def HMP(p_values):
    p_values = np.array(p_values)
    L = len(p_values)
    harmonic_mean = L / np.sum(1 / p_values)
    combined_p = moyal.cdf(-2 * np.log(harmonic_mean), 1, L)
    return combined_p


def EmpiricalBrownsMethod(data_matrix, extra_info=False):
    covar_matrix = CalculateCovariances(data_matrix)
    return CombinePValuesEBM(covar_matrix, data_matrix, extra_info)


def KostsMethod(data_matrix, extra_info=False):
    covar_matrix, _ = CalculateKostCovarianceEBM(data_matrix)
    return CombinePValuesEBM(covar_matrix, data_matrix, extra_info)


def BrownsMethodbyYang(data_matrix, extra_info=False):
    delta_matrix, _ = CalculateCovarianceYang(data_matrix)
    return CombinePValuesYang(delta_matrix, data_matrix, extra_info)


# Input: raw data vector (of one variable) with no missing samples. May be a list or an array.
# Output Transformed data vector w.
def TransformData(data_vector):
    m = np.mean(data_vector)
    sd = np.std(data_vector)
    s = [(d - m) / sd for d in data_vector]
    W = lambda x: -2 * np.log(ECDF(s)(x))
    return np.array([W(x) for x in s])


# Input: An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample.
# Should be of type numpy.array.
# Note: Method does not deal with missing values within the data.
# Output: An m x m matrix of pairwise covariances between transformed raw data vectors
def CalculateCovariances(data_matrix):

    covar_matrix = np.cov(data_matrix)

    return covar_matrix


# Input: A m x m numpy array of covariances between transformed data vectors and a vector of m p-values to combine.
# Output: A combined P-value.
# If extra_info == True: also returns the p-value from Fisher's method, the scale factor c, and the new degrees of freedom from Brown's Method
def CombinePValuesEBM(covar_matrix, p_values, extra_info=False):
    m = covar_matrix.shape[0]
    print ("m", m)
    df_fisher = 2.0 * m
    Expected = 2.0 * m
    cov_sum = 0
    for i in range(m):
        for j in range(i + 1, m):
            cov_sum += covar_matrix[i, j]

    print("cov sum", cov_sum)
    Var = 4.0 * m + 2 * cov_sum
    c = Var / (2.0 * Expected)
    print(c)
    df_brown = 2.0 * Expected ** 2 / Var
    if df_brown > df_fisher:
        df_brown = df_fisher
        c = 1.0

    x = 2.0 * sum([-np.log(np.clip(p, 1e-15, 1 - 1e-15)) for p in p_values])
    # print "x", x
    p_brown = chi2_cdf(df_brown, 1.0 * x / c)
    p_fisher = chi2_cdf(df_fisher, 1.0 * x)

    if extra_info:
        return p_brown, p_fisher, c, df_brown
    else:
        return p_brown


def CombinePValuesYang(delta_matrix, p_values, extra_info=False):
    m = int(delta_matrix.shape[0])
    mean = 2.0 * m
    df_fisher = 2.0 * m
    delta_sum = 0
    for i in range(m):
        for j in range(i + 1, m):
            delta_sum += delta_matrix[i, j]

    Var = 4.0 * m + delta_sum
    v = 2 * (mean ** 2 / Var)
    gamma_value = Var / mean
    T = 2.0 * sum([-np.log(np.clip(p, 1e-15, 1 - 1e-15)) for p in p_values])
    p_fisher = chi2_cdf(df_fisher, 1.0 * T)
    p_yang = 1 - gamma.cdf(T, v / 2, scale=2 * gamma_value)

    if extra_info:
        return p_yang, p_fisher
    else:
        return p_yang


# Input: An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample. Should be of type numpy.array
#       A vector of m P-values to combine. May be a list or of type numpy.array.
# Output: A combined P-value using Kost's Method.
#        If extra_info == True: also returns the p-value from Fisher's method, the scale factor c, and the new degrees of freedom from Brown's Method

def KostsMethoDYang(data_matrix):
    delta_matrix = CalculateCovarianceYang(data_matrix)
    return CombinePValuesYang(delta_matrix, data_matrix)


# Input correlation between two n x n data vectors.
# Output: Kost's approximation of the covariance between the -log cumulative
# distributions. This is calculated with a cubic polynomial fit.
def KostPolyFitEBM(cor):
    a1, a2, a3 = 3.263, 0.710, .027  # Kost cubic coeficients
    return a1 * cor + a2 * cor ** 2 + a3 * cor ** 3


def PolyFitYang(cor):
    c1, c2, c3, c4, c5 = 3.9081, 0.0313, 0.1022, -0.1378, 0.0941
    return c1 * cor ** 2 + c2 * cor ** 4 + c3 * cor ** 6 + c4 * cor ** 8 + c5 * cor ** 10


# Input: An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample. Should be of type numpy.array.
# Note: Method does not deal with missing values within the data.
# Output: An m x m matrix of pairwise covariances between the data vectors calculated using Kost's polynomial fit and numpy's pearson correlation function.
def CalculateKostCovarianceEBM(data_matrix):
    n = data_matrix.shape[1]  # Number of columns
    covar_matrix = np.zeros((n, n))
    max_cor = -np.inf
    for i in range(n):
        for j in range(i + 1, n):
            cor, p_val = pearsonr(data_matrix[:, i], data_matrix[:, j])  # Correlate columns
            covar = KostPolyFitEBM(cor)
            covar_matrix[i, j] = covar
            covar_matrix[j, i] = covar

            max_cor = max(max_cor, cor)

    return covar_matrix, max_cor


def CalculateCovarianceYang(data_matrix):
    c1 = 3.9081
    n = data_matrix.shape[1]  # Number of columns
    delta_matrix = np.zeros((n, n))
    corr_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            cor, _ = pearsonr(data_matrix[:, i], data_matrix[:, j])  # Correlate columns
            corr_matrix[i, j] = cor
            biased_corrected_cor = cor * (1 + (1 - cor ** 2 / 2 * (n - 3)))

            f_r = PolyFitYang(biased_corrected_cor)
            bias = (c1 / n) * (1 - biased_corrected_cor ** 2) ** 2

            unbiased_delta = f_r - bias
            delta_matrix[i, j] = unbiased_delta
            delta_matrix[j, i] = unbiased_delta

    return delta_matrix, corr_matrix


def Bonferronis_correction_g(p_values):

    g = len(p_values)
    _, ICC = CalculateKostCovarianceEBM(p_values)
    print(ICC)
    g_adjusted = (g + 1) - (1 + (g - 1) * ICC)

    return g_adjusted

def effective_number_of_tests_cheverud_nyholt(eigenvalues):
    k = len(eigenvalues)
    sample_variance = np.var(eigenvalues)
    print("var", sample_variance)
    effective_number_of_tests = 1 + (k - 1) * (1 - sample_variance / k)
    return effective_number_of_tests


def h_function(x):
    if x >= 1:
        h = 1 + (x - math.floor(x))
    else:
        h = x - math.floor(x)
    return h


def effective_number_of_tests_li_ji(eigenvalues):
    abs_eigenvalues = np.abs(eigenvalues)
    effective_number_of_tests = np.sum([h_function(abs_lambda) for abs_lambda in abs_eigenvalues])

    return effective_number_of_tests


def effective_number_of_tests_gao(eigenvalues, C=0.995):
    k = len(eigenvalues)
    total_sum = np.sum(eigenvalues)

    for x in range(k):
        x_sum = np.sum(eigenvalues[:x + 1])  # Sum of eigenvalues up to index x
        if x_sum / total_sum > C:
            effective_number_of_tests = x + 1  # Add 1 to convert from index to count
            break

    return effective_number_of_tests


def effective_number_of_tests_galwey(eigenvalues):
    k = len(eigenvalues)
    lambda_prime = np.maximum(0, eigenvalues)

    # Step 3: Compute sum of squares of Î»'_i
    sum_squared_lambda_prime = (np.sum(np.sqrt(lambda_prime))) ** 2

    # Step 4: Divide by the sum of all Î»'_i
    effective_number_of_tests = sum_squared_lambda_prime / np.sum(lambda_prime)

    return effective_number_of_tests


def bonferroni_method_with_effective_tests(p_values, effective_number_of_tests):
    p_values = np.array(p_values)

    # Ensure that all p-values are within the valid range (0, 1)
    min_p_value = min(p_values)
    combined_p = min(1, min_p_value * effective_number_of_tests)

    return combined_p


def default():
    print("Invalid choice")

def transform_to_pvalues(z_scores):
    # Transform z-scores to two-sided p-values
    p_values = [2 * (1 - norm.cdf(abs(z))) for z in z_scores]
    return p_values


def read_from_file_transform(file_path, data_type):
    SNP_values = []

    with open(file_path, 'r') as file:
        content = file.readlines()

    for line in content:
        elements = line.strip().split()
        SNP, values = elements[0], list(map(float, elements[1:]))
        if data_type == 'z-scores':
            values = transform_to_pvalues(values)
        SNP_values.append((SNP, values))


    return SNP_values


def read_matrix_from_file(input_file_path):
    matrix = []
    with open(input_file_path, 'r') as file:
        for line in file:
            values = line.strip().split()
            row = [float(value) for value in values]
            matrix.append(row)
    return matrix



def combine_pvalues(combine_function, SNP_values, metanalysis_needed, file):
    p_values = []
    if metanalysis_needed:
        SNP_dict = {}
        for SNP, values in SNP_values:
            SNP_dict.setdefault(SNP, []).append(values)

        for SNP, values in SNP_dict.items():
            if len(values) > 1:
                grouped_values = list(zip(*values))
                combined_results = [combine_function(np.array(grouped_value)) for grouped_value in grouped_values]
                SNP_dict[SNP] = np.float_(combined_results)
            else:
                # Directly use the single set of values
                SNP_dict[SNP] = np.float_(values[0])

            p_values.append({'SNP': SNP, 'values': SNP_dict[SNP]})

        for p_value in p_values:
            print(f"Combining p-values for SNP {p_value['SNP']}: {p_value['values']}")  # Debug print
            comb_p = combine_function(p_value['values'])
            file.write(f"Combined p-values for SNP {p_value['SNP']}: {comb_p}\n")
    else:
        for SNP, values in SNP_values:
            values = np.float_(values)  # Ensure values are converted to float
            p_values.append({'SNP': SNP, 'values': values})

        for p_value in p_values:
            print(f"Combining p-values for SNP {p_value['SNP']}: {p_value['values']}")  # Debug print
            comb_p = combine_function(np.array(p_value['values']))
            file.write(f"Combined p-values for SNP {p_value['SNP']}: {comb_p}\n")

    return p_values


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process text files and data type.')
    parser.add_argument('input_file_path', type=str, help='Path to the input file')
    parser.add_argument('output_file_path', type=str, help='Path to the output file')
    parser.add_argument('data_type', type=str, choices=['p-values', 'z-scores'], help="Type of data in the file ('p-values' or 'z-scores')")
    parser.add_argument('meta_choice', type=str, choices=['Yes', 'No'], help="Do you want the program to perform meta-analysis? ('Yes' or 'No')")
    parser.add_argument('correlation_matrix_path', type=str, nargs='?', default=None,
                        help="Path to the correlation matrix file (optional)")
    parser.add_argument('methods', type=str, nargs='+', choices=['logit', 'meanp', 'fisher','stouffer','invchi','binomial','binomial', 'cct', 'minp', 'mcm', 'cmc', 'bonferroni', 'ebm', 'kost', 'yang']
                        , help="Which method(s) would you like to use to combine your p-values(or z-scores)?")
    args = parser.parse_args()


    SNP_values = read_from_file_transform(args.input_file_path, args.data_type.lower())
    metanalysis_needed = args.meta_choice.lower() == 'yes'

    correlation_matrix = None
    if args.correlation_matrix_path:
        print(f"Using correlation matrix file: {args.correlation_matrix_path}")

    else:
        print("No correlation matrix file provided.")

    print(f"File Path: {args.input_file_path}")
    print(f"Output File Path: {args.output_file_path}")
    print(f"Data Type: {args.data_type}")
    print(f"Meta-analysis: {args.meta_choice}")
    print(f"Correlation Matrix: {args.correlation_matrix_path}")
    print(f"Methods: {args.methods}")


    method_functions = {
        'logit': logit,
        'meanp': meanp,
        'fisher': fisher_method,
        'stouffer': stouffer,
        'invchi': inverse_chi2,
        'binomial': binomial_test,
        'cct': cauchy_method,
        'minp': minP,
        'cmc': CMC,
        'mcm': MCM,
        'hmp': HMP,
        'ebm': EmpiricalBrownsMethod,
        'kost': KostsMethod,
        'yang': BrownsMethodbyYang,
        'bonferroni': bonferroni_method_with_effective_tests
    }
    selected_methods = args.methods

    method_names = {
        'logit': "Logitp Method",
        'meanp': "Meanp Method",
        'fisher': "Fisher's Method",
        'stouffer': "Stouffer's Method",
        'invchi': "Inverse Chi2 Method",
        'binomial': "Binomial Test",
        'cct': "Cauchy Method (CCT)",
        'minp': "MinP Method (Tippett's Method)",
        'cmc': "CMC (CCT-MinP-CCT)",
        'mcm': "MCM (MinP-CCT-MinP)",
        'hmp': "Harmonic Mean P-value (HMP)",
        'ebm': "Empirical Brown's Method (EBM)",
        'kost': "Empirical Brown's Method (EBM) by Kost",
        'yang': "Brown's Method by Yang",
        'bonferroni': "Bonferroni with effective number of tests"
    }

    print("\n****** Combining methods for dependent and independent P-values ******\n")
    for key, value in method_names.items():
        print(f"{key}: *** {value} ***")


    with open(args.output_file_path, "w") as file:
        for method in selected_methods:
            if method in method_functions:
                file.write(f"\nSelected Method: *** {method_names[method]} ***\n")
                combine_function = method_functions[method]

                if method in ['ebm', 'kost', 'yang', 'bonferroni']:
                    values = np.array([values for _, values in SNP_values])
                    print(values)
                    SNP = np.array([values for SNP, _ in SNP_values])
                    g_adjusted = Bonferronis_correction_g(values)

                    if method in ['ebm', 'kost', 'yang']:

                        # Get combined values for each SNP
                        combined_values = file.write(str(combine_function(values, extra_info=False)))


                    elif method in ['bonferroni']:

                        if correlation_matrix is None:
                            corr_matrix = np.corrcoef(values, rowvar=False)
                        else:
                            corr_matrix = read_matrix_from_file(args.correlation_matrix_path)


                        eigenvalues, _ = np.linalg.eig(corr_matrix)
                        print(eigenvalues)

                        sorted_eigenvalues = np.sort(eigenvalues)[::-1]

                        cn = effective_number_of_tests_cheverud_nyholt(sorted_eigenvalues)
                        gao = effective_number_of_tests_gao(sorted_eigenvalues)
                        gal = effective_number_of_tests_galwey(sorted_eigenvalues)
                        li_ji = effective_number_of_tests_li_ji(sorted_eigenvalues)


                        file.write(f"Effective Number of Tests (Cheverud-Nyholt): {cn}\n")


                        for SNP, values in SNP_values:
                            print(values)

                            combined_p_value_cn = combine_function(values, cn)
                            file.write(f"Combined p-value for {SNP}: {combined_p_value_cn}\n")

                        file.write(f"Effective Number of Tests (Gao): {gao}\n")
                        for SNP, values in SNP_values:
                            combined_p_value_gao = combine_function(values, gao)
                            file.write(f"Combined p-value for {SNP}: {combined_p_value_gao}\n")

                        file.write(f"Effective Number of Tests (Galwey): {gal}\n")
                        for SNP, values in SNP_values:
                            combined_p_value_gal = combine_function(values, gal)
                            file.write(f"Combined p-value for {SNP}: {combined_p_value_gal}\n")

                        file.write(f"Effective Number of Tests (Li-Ji): {li_ji}\n")
                        for SNP, values in SNP_values:
                            combined_p_value_li_ji = combine_function(values, li_ji)
                            file.write(f"Combined p-value for {SNP}: {combined_p_value_li_ji}\n")


                        file.write(f"Bonferroni's correction with g* =  {g_adjusted}:\n")
                        for SNP, values in SNP_values:
                            combined_p_bc = combine_function(values, g_adjusted)
                            file.write(f"Combined p-value for {SNP}: {combined_p_bc}\n")



                else:
                    p_values = combine_pvalues(combine_function, SNP_values, metanalysis_needed, file)
            else:
                file.write(f"Invalid method choice: {method}\n")
