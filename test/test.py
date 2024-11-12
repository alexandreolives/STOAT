import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact
import time

def chi2_test(df) -> float:
    """Calculate p_value from list of dataframe using chi-2 test"""
    # Check if dataframe has at least 2 columns and more than 0 counts in every cell
    if df.shape[1] >= 2 and np.all(df.sum(axis=0)) and np.all(df.sum(axis=1)):
        try:
            # Perform Chi-Square test
            p_value = chi2_contingency(df)[1]  # chi2_contingency returns a tuple (chi2, p, dof, expected)
        except ValueError as e:
            p_value = "Error"
    else:
        p_value = "N/A"
    return p_value

def fisher_test(df) -> float:
    """Calculate p_value using Fisher exact test"""
    try:
        p_value = fisher_exact(df)[1]  # fisher_exact returns a tuple (oddsratio, p_value)
    except ValueError as e:
        p_value = 'N/A'
    return p_value

# 2x2 Table for testing
data = [[1982, 3018], [2056, 2944]]
df = pd.DataFrame(data, columns=['Group1', 'Group2'])

# Timing and testing Chi2
start_time = time.time()
chi2_p_value = chi2_test(df)
chi2_time = time.time() - start_time

# Timing and testing Fisher's Exact Test
start_time = time.time()
fisher_p_value = fisher_test(df)
fisher_time = time.time() - start_time

print(f"Chi-square p-value: {chi2_p_value}, Time taken: {chi2_time:.6f} seconds")
print(f"Fisher's exact p-value: {fisher_p_value}, Time taken: {fisher_time:.6f} seconds")


















from functools import reduce
import numpy as np
from ._stats_py import power_divergence
from scipy._lib._bunch import _make_tuple_bunch

def margins(a):

    margsums = []
    ranged = list(range(a.ndim))
    for k in ranged:
        marg = np.apply_over_axes(np.sum, a, [j for j in ranged if j != k])
        margsums.append(marg)
    return margsums


def expected_freq(observed):

    # Typically `observed` is an integer array. If `observed` has a large
    # number of dimensions or holds large values, some of the following
    # computations may overflow, so we first switch to floating point.
    observed = np.asarray(observed, dtype=np.float64)

    # Create a list of the marginal sums.
    margsums = margins(observed)

    # Create the array of expected frequencies.  The shapes of the
    # marginal sums returned by apply_over_axes() are just what we
    # need for broadcasting in the following product.
    d = observed.ndim
    expected = reduce(np.multiply, margsums) / observed.sum() ** (d - 1)
    return expected


Chi2ContingencyResult = _make_tuple_bunch(
    'Chi2ContingencyResult',
    ['statistic', 'pvalue', 'dof', 'expected_freq'], []
)


def chi2_contingency(observed, correction=True, lambda_=None):

    observed = np.asarray(observed)
    if np.any(observed < 0):
        raise ValueError("All values in `observed` must be nonnegative.")
    if observed.size == 0:
        raise ValueError("No data; `observed` has size 0.")

    expected = expected_freq(observed)
    if np.any(expected == 0):
        # Include one of the positions where expected is zero in
        # the exception message.
        zeropos = list(zip(*np.nonzero(expected == 0)))[0]
        raise ValueError("The internally computed table of expected "
                         f"frequencies has a zero element at {zeropos}.")

    # The degrees of freedom
    dof = expected.size - sum(expected.shape) + expected.ndim - 1

    if dof == 0:
        # Degenerate case; this occurs when `observed` is 1D (or, more
        # generally, when it has only one nontrivial dimension).  In this
        # case, we also have observed == expected, so chi2 is 0.
        chi2 = 0.0
        p = 1.0
    else:
        if dof == 1 and correction:
            # Adjust `observed` according to Yates' correction for continuity.
            # Magnitude of correction no bigger than difference; see gh-13875
            diff = expected - observed
            direction = np.sign(diff)
            magnitude = np.minimum(0.5, np.abs(diff))
            observed = observed + magnitude * direction

        chi2, p = power_divergence(observed, expected,
                                   ddof=observed.size - 1 - dof, axis=None,
                                   lambda_=lambda_)

    return Chi2ContingencyResult(chi2, p, dof, expected)