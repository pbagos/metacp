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

def logit(values, data_type = "z-scores"):
    values = np.array(values)
    k = len(values)
    C = np.sqrt(k * np.pi ** 2 * (5 * k + 2) / (3 * (5 * k + 4)))
    df = 2 * k
    if data_type == "z-scores":
        p1, p2 = transform_to_pvalues(values)

        t1_value = -np.sum(np.log((p1) / (p2)))
        combined_p1 = 2 * (1 - t.cdf(np.abs(t1_value / C), df))
        t2_value = -np.sum(np.log((p2) / (p1)))
        combined_p2 = 2 * (1 - t.cdf(np.abs(t2_value / C), df))

        p_final = 2 * min(combined_p1, combined_p2)
    else:
        # Ensure that all probabilities are within the valid range (0, 1)
        p_values = np.clip(values, 1e-16, 1 - 1e-16)

        # Compute the logit transform
        t_value = -np.sum(np.log((p_values) / (1 - p_values)))

        p_final = 2 * (1 - t.cdf(np.abs(t_value / C), df))

    return p_final


def meanp(values, data_type="z-scores"):
    k = len(values)
    if data_type == "z-scores":
        p1, p2 = transform_to_pvalues(values)

        # Compute the mean of the p-values
        z1_value = (0.5 - np.mean(p1)) * math.sqrt(12 * k)
        combined_p1 = 1 - norm.cdf(z1_value)
        z2_value = (0.5 - np.mean(p2)) * math.sqrt(12 * k)
        combined_p2 = 1 - norm.cdf(z2_value)
        p_final = 2 * min(combined_p1, combined_p2)

    else:

        p_values = np.array(values)

        # Ensure that all p-values are within the valid range (0, 1)
        p_values = np.clip(p_values, 1e-09, 1 - 1e-09)
        # Compute the mean of the p-values
        z_value = (0.5 - np.mean(p_values)) * math.sqrt(12 * k)
        p_final = 1 - norm.cdf(z_value)

    return p_final


def fisher_method(values, data_type = "z-scores"):
    df = 2 * len(values)
    if data_type == "z-scores":
        p1, p2 = transform_to_pvalues(values)

        chi_squared1 = -2 * np.sum(np.log(p1))
        combined_p1 = 1 - chi2.cdf(chi_squared1, df)
        chi_squared2 = -2 * np.sum(np.log(p2))
        combined_p2 = 1 - chi2.cdf(chi_squared2, df)

        p_final = 2 * min(combined_p1, combined_p2)
    else:
        # Convert input to numpy array for numerical operations
        p_values = np.array(values)

        # Ensure that all p-values are within the valid range (0, 1)
        p_values = np.clip(p_values, 1e-16, 1 - 1e-16)
        # Compute the combined test statistic using Fisher's method
        chi_squared = -2 * np.sum(np.log(p_values))

        # Calculate the combined p-value using the chi-squared distribution
        p_final = 1 - chi2.cdf(chi_squared, df)

    return p_final


def lancaster_method(values, data_matrix, data_type="z-scores", weight_matrix=None):
    weight_matrix = np.array(weight_matrix, dtype=float)
    # Compute the degrees of freedom (sum of the weights or sample sizes)
    df = np.sum(weight_matrix)

    if data_type == "z-scores":
        values = np.array(values)
        p1, p2 = transform_to_pvalues(values)



        # Compute chi-square statistics using the inverse chi-square (quantile) of each p-value
        chi_statistic1 = np.sum(chi2.ppf(1 - p1, weight_matrix))  # Applying weight for each study
        chi_statistic2 = np.sum(chi2.ppf(1 - p2, weight_matrix))  # Applying weight for each study

        # Combined p-values using chi-square distribution
        combined_p1 = 1 - chi2.cdf(chi_statistic1, df)
        combined_p2 = 1 - chi2.cdf(chi_statistic2, df)
        p_final = (2 * min(combined_p1, combined_p2))


    else:
        # If the data type is not z-scores, assume we are directly working with p-values
        p_values = np.array(values)

        # Ensure p-values are within a valid range to avoid numerical issues
        p_values = np.clip(p_values, 1e-16, 1 - 1e-16)

        # Compute the chi-square statistic by applying the inverse chi-square to each p-value
        chi_statistic = np.sum(chi2.ppf(1 - p_values, weight_matrix))  # Applying weight for each study

        # Combined p-value using the chi-square distribution
        p_final = 1 - chi2.cdf(chi_statistic, df)

    return p_final


def stouffer(values, data_type="z-scores"):
    # Convert input to numpy array for numerical operations
    values = np.array(values)

    if data_type == "p-values":
        z_scores = norm.ppf(1 - values)  # norm.ppf =inverse cumulative distribution function of normal distribution
    else:
        z_scores = values

    # z_scores = norm.ppf(1 - p_values) # norm.ppf =inverse cumulative distribution function of normal distribution
    combined_z = np.sum(z_scores / math.sqrt(len(values)))
    combined_p = 1 - norm.cdf(combined_z)  # norm.cdf = cumulative distribution function of normal distribution

    return combined_p

def weighted_stouffer(values, data_matrix, data_type="z-scores", weight_matrix = None):
    values = np.array(values)

    weight_matrix = np.array(weight_matrix, dtype=float)

    if data_type == "p-values":
        z_scores = norm.ppf(1 - values)  # norm.ppf =inverse cumulative distribution function of normal distribution
    else:
        z_scores = values

    # z_scores = norm.ppf(1 - p_values) # norm.ppf =inverse cumulative distribution function of normal distribution
    combined_z = np.sum(weight_matrix * z_scores / math.sqrt(np.sum(weight_matrix**2)))
    combined_p = 1 - norm.cdf(combined_z)  # norm.cdf = cumulative distribution function of normal distribution

    return combined_p

def inverse_chi2(values, data_type="z-scores"):
    values = np.array(values)
    df = len(values)
    if data_type == "z-scores":
        p1, p2 = transform_to_pvalues(values)

        chi_squared1 = np.sum(chi2.ppf((1 - p1), 1))  # ppf = inverse of cdf
        combined_p1 = 1 - chi2.cdf(chi_squared1, df)  # chi2.cdf = cumulative distribution function of chi2
        chi_squared2 = np.sum(chi2.ppf((1 - p2), 1))  # ppf = inverse of cdf
        combined_p2 = 1 - chi2.cdf(chi_squared2, df)
        p_final = 2 * min(combined_p1, combined_p2)
    else:
        # Convert input to numpy array for numerical operations
        # Ensure that all p-values are within the valid range (0, 1)
        p_values = np.clip(values, 1e-16, 1 - 1e-16)

        chi_squared = np.sum(chi2.ppf((1 - p_values), 1))  # ppf = inverse of cdf
        p_final = 1 - chi2.cdf(chi_squared, df)  # chi2.cdf = cumulative distribution function of chi2

    return p_final


def binomial_test(values, data_type="z-scores"):
    values = np.array(values)
    k = len(values)
    alpha = 0.05
    if data_type == "z-scores":
        p1, p2 = transform_to_pvalues(values)
        # Count the number of significant p-values
        r1 = sum((p1 < alpha) for p in p1)
        # Calculate the binomial probability of observing at most num_successes successes
        combined_p1 = 0
        for x1 in range(r1, k + 1):
            combined_p1 += math.factorial(k) / (math.factorial(x1) * math.factorial(k - x1)) * (alpha ** x1) * (
                        (1 - alpha) ** (k - x1))
        r2 = sum((p2 < alpha) for p in p2)
        # Calculate the binomial probability of observing at most num_successes successes
        combined_p2 = 0
        for x2 in range(r2, k + 1):
            combined_p2 += math.factorial(k) / (math.factorial(x2) * math.factorial(k - x2)) * (alpha ** x2) * (
                (1 - alpha) ** (k - x2))

        p_final = 2 * min(combined_p1, combined_p2)

    else:
        p_values = np.array(values)

        # Count the number of significant p-values
        r = sum((p < alpha) for p in p_values)
        # Calculate the binomial probability of observing at most num_successes successes
        p_final = 0
        for x in range(r, k + 1):
            p_final += math.factorial(k) / (math.factorial(x) * math.factorial(k - x)) * (alpha ** x) * ((1 - alpha) ** (k - x))

    return p_final


def cauchy_cdf(x, x0=0, gamma=1):
    """Cumulative distribution function (CDF) of the Cauchy distribution."""

    return 0.5 + (np.arctan((x - x0) / gamma) / np.pi)


def cauchy_method(values, data_type="z-scores"):
    k = len(values)
    values = np.array(values)

    if data_type == "z-scores":
        # Transform to p-values
        p1, p2 = transform_to_pvalues(values)

        T1 = np.tan((0.5 - p1) * np.pi)
        t1 = np.sum(T1) / k
        # Calculate the combined p-value using the Cauchy distribution
        combined_p1 = 1 - cauchy_cdf(t1)

        T2 = np.tan((0.5 - p2) * np.pi)
        t2 = np.sum(T2) / k
        # Calculate the combined p-value using the Cauchy distribution
        combined_p2 = 1 - cauchy_cdf(t2)

        p_final = 2 * min(combined_p1, combined_p2)
    else:

        # Ensure that all p-values are within the valid range (0, 1)
        p_values = np.clip(values, 1e-15, 1 - 1e-15)

        T = np.tan((0.5 - p_values) * np.pi)
        t = np.sum(T) / k
        # Calculate the combined p-value using the Cauchy distribution
        p_final = 1 - cauchy_cdf(t)

    return p_final


def minP(values, data_type="z-scores"):
    values = np.array(values)
    k = len(values)

    if data_type == "z-scores":

        p1, p2 = transform_to_pvalues(values)

        combined_p1 = 1 - (1 - np.min(p1)) ** k
        combined_p2 = 1 - (1 - np.min(p2)) ** k

        p_final = 2 * min(combined_p1, combined_p2)
    else:
        # Ensure that all p-values are within the valid range (0, 1)
        p_values = np.clip(values, 1e-15, 1 - 1e-15)

        # Calculate the combined p-value as the minimum p-value
        p_final = 1 - (1 - np.min(p_values)) ** k

    return p_final



def CMC(values, data_type="z-scores"):
    # Compute individual p-values
    p_value_cauchy = cauchy_method(values, data_type="z-scores")
    p_value_minp = minP(values, data_type="z-scores")

    # Combine p_value_cauchy and p_value_minp into one array
    combined_values = np.array([p_value_cauchy, p_value_minp])

    # Apply the Cauchy method logic directly to the combined array
    k = len(combined_values)
    combined_values = np.clip(combined_values, 1e-15, 1 - 1e-15)  # Ensure p-values are in range [0, 1]

    # Apply the Cauchy transformation
    T = np.tan((0.5 - combined_values) * np.pi)
    t = np.sum(T) / k

    # Calculate the combined p-value using the Cauchy distribution
    combined_p = 1 - cauchy_cdf(t)

    return combined_p


def MCM(values, data_type="z-scores"):
    p_value_cauchy = cauchy_method(values)
    p_value_minp = minP(values)
    combined_p = 2 * min(p_value_cauchy, p_value_minp, 0.5)  # pMCM = 2 min{pCCT , pMinP, 0.5}

    return combined_p


def HMP(values, data_type="z-scores"):
    values = np.array(values)
    k = len(values)

    if data_type == "z-scores":
        p1, p2 = transform_to_pvalues(values)

        combined_p1 = k / np.sum(1 / p1)
        combined_p2 = k / np.sum(1 / p2)

        p_final = 2 * min(combined_p1, combined_p2)
    else:
        p_final = k / np.sum(1 / p_values)

    return p_final


def EmpiricalBrownsMethod(values, data_matrix, data_type="z-scores"):
    values = np.array(values)
    n = data_matrix.shape[1]
    df_fisher = 2.0 * n
    Expected = 2.0 * n

    if data_type == "z-scores":
        p1, p2 = transform_to_pvalues(values)
        print("p1",p1)
        print("p2",p2)

        covar_matrix1 = np.cov(p1)
        covar_matrix2 = np.cov(p2)

        cov_sum1 = np.sum(covar_matrix1)
        cov_sum2 = np.sum(covar_matrix2)

        Var1 = 4.0 * n + 2 * cov_sum1
        Var2 = 4.0 * n + 2 * cov_sum2

        c1 = Var1 / (2.0 * Expected)
        df_brown1 = (2.0 * Expected ** 2) / Var1
        df_brown1 = min(df_brown1, df_fisher)
        c1 = 1.0 if df_brown1 == df_fisher else c1

        c2 = Var2 / (2.0 * Expected)
        df_brown2 = (2.0 * Expected ** 2) / Var2
        df_brown2 = min(df_brown2, df_fisher)
        c2 = 1.0 if df_brown2 == df_fisher else c2

        chi_squared1 = 2.0 * np.sum(-np.log(np.clip(p1, 1e-15, 1 - 1e-15)))
        combined_p1 = 1 - chi2.cdf(chi_squared1 / c1, df_brown1)
        print("combined_p1", combined_p1)

        chi_squared2 = 2.0 * np.sum(-np.log(np.clip(p2, 1e-15, 1 - 1e-15)))
        combined_p2 = 1 - chi2.cdf(chi_squared2 / c2, df_brown2)
        print("combined_p2", combined_p2)
        p_brown_final = 2 * min(combined_p1, combined_p2)

    else:
        covar_matrix = np.cov(data_matrix)

        cov_sum = np.sum(covar_matrix)
        Var = 4.0 * n + 2 * cov_sum
        c = Var / (2.0 * Expected)
        df_brown = (2.0 * Expected ** 2) / Var
        df_brown = min(df_brown, df_fisher)
        c = 1.0 if df_brown == df_fisher else c

        x = 2.0 * np.sum(-np.log(np.clip(values, 1e-15, 1 - 1e-15)))
        p_brown_final = chi2.cdf(x / c, df_brown)

    return p_brown_final


def KostsMethod(values, data_matrix, data_type="z-scores"):
    n = len(values)
    covar_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            cor, _ = pearsonr(data_matrix[:, i], data_matrix[:, j])
            print(f"Correlation between {i} and {j}: {cor}")
            covar = 3.263 * cor + 0.710 * cor ** 2 + 0.027 * cor ** 3
            covar_matrix[i, j] = covar_matrix[j, i] = covar
    return EmpiricalBrownsMethod(covar_matrix,data_matrix, data_type)


def BrownsMethodbyYang(values, data_matrix, data_type="z-scores"):
    c1 = 3.9081
    values = np.array(values)
    n = data_matrix.shape[1]  # Number of columns (SNPs or variables)


    if data_type == "z-scores":
        p_values = 2 * (1 - norm.cdf(abs(values)))
    else:
        p_values = values

    delta_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            cor, _ = pearsonr(data_matrix[:,i], data_matrix[:,j])
            print(f"Correlation between {i} and {j}: {cor}")

            biased_corrected_cor = cor * (1 + (1 - cor ** 2 / 2 * (n - 3)))
            print(f"Biased corrected correlation: {biased_corrected_cor}")

            f_r = (3.9081 * biased_corrected_cor ** 2 + 0.0313 * biased_corrected_cor ** 4 +
                   0.1022 * biased_corrected_cor ** 6 - 0.1378 * biased_corrected_cor ** 8 +
                   0.0941 * biased_corrected_cor ** 10)
            bias = (c1 / n) * (1 - biased_corrected_cor ** 2) ** 2
            print(f"f_r: {f_r}, bias: {bias}")

            delta_matrix[i, j] = delta_matrix[j, i] = f_r - bias

    mean, df_fisher = 2.0 * n, 2.0 * n
    delta_sum = np.sum(delta_matrix)
    print("delta sum =", delta_sum)
    Var = 4.0 * n + delta_sum
    v_gamma = 2 * (mean ** 2 / Var)
    gamma_value = Var / mean

    T = 2.0 * np.sum(-np.log(np.clip(p_values, 1e-15, 1 - 1e-15)))

    p_yang_final = 1 - gamma.cdf(T, v_gamma / 2, scale=2 * gamma_value)

    return p_yang_final

def weighted_correlated_Stouffer(values, data_matrix, data_type="z-scores", weight_matrix=None):
    data_matrix = np.array(data, dtype=float)
    print(f"Input type of values: {type(data_matrix)}")  # Debugging line
    print(f"Input values: {data}")  # Debugging line
    k = data_matrix.shape[1]
    if weight_matrix is None:
        print("DEBUG: No weight matrix provided, using default weights.")
        weight_matrix = np.ones(len(values))  # Default to equal weights

    weight_matrix = np.array(weight_matrix, dtype=float)  # Convert to NumPy array
    correlation_matrix = np.eye(k)  # Initialize with identity matrix (diagonal = 1)
    for i in range(k):
        for j in range(i + 1, k):
            cor, _ = pearsonr(data_matrix[:, i], data_matrix[:, j])
            correlation_matrix[i, j] = cor
            correlation_matrix[j, i] = cor  # Symmetric matrix
            print(f"Correlation between {i} and {j}: {cor}")

    # Convert p-values to z-scores if needed
    if data_type == "p-values":
        z_scores = norm.ppf(1 - values)
    else:
        z_scores = values

    # Compute weighted sum of z-scores
    weighted_z = np.sum(weight_matrix * z_scores)

    # Compute weighted variance
    weighted_variance = np.sum(weight_matrix[:, None] * weight_matrix[None, :] * correlation_matrix)
    weighted_std_dev = np.sqrt(weighted_variance)

    # Compute combined z-score and p-value
    combined_z = weighted_z / weighted_std_dev
    combined_p_value = 1 - norm.cdf(combined_z)

    return combined_p_value

def correlated_Stouffer(values, data_matrix, data_type="zscores"):
    data_matrix = np.array(data, dtype=float)
    values = np.array(values)
    k = data_matrix.shape[1]  # Number of tests (columns)

    if data_type == "p-values":
        z_scores = norm.ppf(1 - values)
    else:
        z_scores = values

    cor_sum = 0
    for i in range(k):
        for j in range(i + 1, k):
            cor, _ = pearsonr(data_matrix[:, i], data_matrix[:, j])  # Correlation between columns (tests)
            print(f"Correlation between {i} and {j}: {cor}")
            cor_sum += cor
    total_variance = k + 2 * cor_sum

    # Compute combined test statistic
    combined_z = np.sum(z_scores) / np.sqrt(total_variance)

    # Compute combined p-value
    combined_p_value = 1 - norm.cdf(combined_z)

    return combined_p_value






def Bonferronis_correction_g(values):
    g = len(values)
    n = values.shape[1]  # Number of columns
    max_cor = -np.inf
    for i in range(n):
        for j in range(i + 1, n):
            cor, _ = pearsonr(values[:, i], values[:, j])  # Correlate columns
            max_cor = max(max_cor, cor)
    ICC = max_cor
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


def bonferroni_method_with_effective_tests(values, effective_number_of_tests, data_type="z-scores"):
    values = np.array(values)

    if data_type == "z-scores":
        p1, p2 = transform_to_pvalues(values)
        min_p_value1 = np.min(p1)
        min_p_value2 = np.min(p2)

        combined_p1 = min(1, min_p_value1 * effective_number_of_tests)
        combined_p2 = min(1, min_p_value2 * effective_number_of_tests)

        p_final = 2 * min(combined_p1, combined_p2)
    else:
        p_values = np.array(values)
        min_p_value = np.min(p_values)
        p_final = min(1, min_p_value * effective_number_of_tests)

    return p_final


def default():
    print("Invalid choice")

def transform_to_pvalues(z_scores):
    # Convert each z-score to p-values for both positive and negative directions
    p1 = np.where(z_scores > 0, 1 - norm.cdf(z_scores), 1)  # p1 for positive direction
    p2 = np.where(z_scores < 0, norm.cdf(z_scores), 1)  # p2 for negative direction
    return p1, p2

def read_from_file(file_path):
    SNP_values = []

    with open(file_path, 'r') as file:
        for line in file:
            elements = line.strip().split()
            if not elements:  # Skip empty lines
                continue
            SNP, values = elements[0], list(map(float, elements[1:]))
            SNP_values.append((SNP, values))

    return SNP_values

def read_matrix_from_file(input_file_path):
    with open(input_file_path, 'r') as file:
        return np.array([[float(value) for value in line.strip().split()] for line in file])


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
            comb_p = combine_function(p_value['values'], args.data_type)
            file.write(f"Combined p-values for SNP {p_value['SNP']}: {comb_p}\n")
    else:
        for SNP, values in SNP_values:
            values = np.float_(values)  # Ensure values are converted to float
            p_values.append({'SNP': SNP, 'values': values})

        for p_value in p_values:
            print(f"Combining p-values for SNP {p_value['SNP']}: {p_value['values']}")
            comb_p = combine_function(np.array(p_value['values']), args.data_type)
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
    parser.add_argument('weights_matrix_path', type=str, nargs='?', default=None,
                        help="Path to the correlation matrix file (optional)")
    parser.add_argument('methods', type=str, nargs='+', choices=['logit', 'meanp', 'fisher','lancaster','stouffer','wstouffer','invchi','binomial', 'cct', 'minp', 'mcm', 'cmc', 'bonferroni', 'ebm', 'kost', 'yang','corstouffer','wcorstouffer']
                        , help="Which method(s) would you like to use to combine your p-values(or z-scores)?")
    args = parser.parse_args()


    SNP_values = read_from_file(args.input_file_path)

    metanalysis_needed = args.meta_choice.lower() == 'yes'

    correlation_matrix = None
    if args.correlation_matrix_path:
        print(f"Using correlation matrix file: {args.correlation_matrix_path}")

    else:
        print("No correlation matrix file provided.")

    weights_matrix = None
    if 'wcorstouffer' in args.methods or 'wstouffer' in args.methods:
        if args.weights_matrix_path:
            print(f"Using weights matrix from: {args.weights_matrix_path}")
            weights_matrix = read_matrix_from_file(args.weights_matrix_path)  # Load user-provided matrix
        else:
            print("ERROR: Weights matrix is required for 'wcorstouffer' but not provided.")


    print(f"File Path: {args.input_file_path}")
    print(f"Output File Path: {args.output_file_path}")
    print(f"Data Type: {args.data_type}")
    print(f"Meta-analysis: {args.meta_choice}")
    print(f"Correlation Matrix: {args.correlation_matrix_path}")
    print(f"Weight Matrix: {args.weights_matrix_path}")
    print(f"Methods: {args.methods}")


    method_functions = {
        'logit': logit,
        'meanp': meanp,
        'fisher': fisher_method,
        'lancaster': lancaster_method,
        'stouffer': stouffer,
        'wstouffer': weighted_stouffer,
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
        'bonferroni': bonferroni_method_with_effective_tests,
        'corstouffer': correlated_Stouffer,
        'wcorstouffer': weighted_correlated_Stouffer
    }
    selected_methods = args.methods

    method_names = {
        'logit': "Logitp Method",
        'meanp': "Meanp Method",
        'fisher': "Fisher's Method",
        'lancaster': "Lancaster's Method",
         'stouffer': "Stouffer's Method",
        'wstouffer': "Weighted Stouffer's Method",
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
        'bonferroni': "Bonferroni with effective number of tests",
        'corstouffer': "Stouffer's Method for dependent tests",
        'wcorstouffer': "Weighted Stouffer's Method for dependent tests"

    }

    print("\n****** Combining methods for dependent and independent P-values ******\n")
    for key, value in method_names.items():
        print(f"{key}: *** {value} ***")


    with open(args.output_file_path, "w") as file:
        for method in selected_methods:
            if method in method_functions:
                file.write(f"\nSelected Method: *** {method_names[method]} ***\n")
                combine_function = method_functions[method]


                if method in ['ebm','kost', 'yang','bonferroni', 'corstoufer', 'wcortouffer', 'lancaster', 'wstouffer']:
                    data = np.array([values for _, values in SNP_values])
                    SNP = np.array([SNP for SNP, _ in SNP_values])
                    g_adjusted = Bonferronis_correction_g(data)

                    if method in ['ebm', 'kost', 'yang', 'corstouffer']:
                        # Get combined values for each SNP
                        for SNP, values in SNP_values:
                            combined_p = combine_function(values, data, args.data_type)
                            file.write(f"Combined p-value for {SNP}: {combined_p}\n")

                    elif method in ['lancaster', 'wstouffer', 'wcorstouffer']:

                        print(f"DEBUG: Loaded weight matrix from file: {weights_matrix}")
                        for SNP, values in SNP_values:
                            combined_p = combine_function(values, data, args.data_type, weights_matrix)
                            file.write(f"Combined p-value for {SNP}: {combined_p}\n")


                    elif method in ['bonferroni']:

                        if correlation_matrix is None:
                            corr_matrix = np.corrcoef(data, rowvar=False)
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
                        for SNP, data in SNP_values:
                            combined_p_value_cn = combine_function(data, cn, args.data_type)
                            file.write(f"Combined p-value for {SNP}: {combined_p_value_cn}\n")

                        file.write(f"Effective Number of Tests (Gao): {gao}\n")
                        for SNP, data in SNP_values:
                            combined_p_value_gao = combine_function(data, gao, args.data_type)
                            file.write(f"Combined p-value for {SNP}: {combined_p_value_gao}\n")

                        file.write(f"Effective Number of Tests (Galwey): {gal}\n")
                        for SNP, data in SNP_values:
                            combined_p_value_gal = combine_function(data, gal, args.data_type)
                            file.write(f"Combined p-value for {SNP}: {combined_p_value_gal}\n")

                        file.write(f"Effective Number of Tests (Li-Ji): {li_ji}\n")
                        for SNP, data in SNP_values:
                            combined_p_value_li_ji = combine_function(data, li_ji,args.data_type)
                            file.write(f"Combined p-value for {SNP}: {combined_p_value_li_ji}\n")


                        file.write(f"Bonferroni's correction with g* =  {g_adjusted}:\n")
                        for SNP, data in SNP_values:
                            combined_p_bc = combine_function(data, g_adjusted)
                            file.write(f"Combined p-value for {SNP}: {combined_p_bc}\n")


                else:

                    p_values = combine_pvalues(combine_function, SNP_values, metanalysis_needed, file)
            else:
                file.write(f"Invalid method choice: {method}\n")
