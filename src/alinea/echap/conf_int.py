import numpy
from scipy.stats import t as t_distribution


def conf_int(lst, perc_conf=95):
    """
    Confidence interval - given a list of values compute the square root of
    the variance of the list (v) divided by the number of entries (n)
    multiplied by a constant factor of (c). This means that I can
    be confident of a result +/- this amount from the mean..
    """
    n, v = len(lst), numpy.var(lst, ddof=1)
    c = t_distribution.interval(perc_conf * 1.0 / 100, n - 1)[1]

    return numpy.sqrt(v / n) * c