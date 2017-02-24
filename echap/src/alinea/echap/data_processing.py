""" Utilities to process / aggregate data
"""

import numpy
import scipy.stats


def conf_int(lst, perc_conf=95):
    """
    Confidence interval - given a list of values compute the square root of
    the variance of the list (v) divided by the number of entries (n)
    multiplied by a constant factor of (c). This means that I can
    be confident of a result +/- this amount from the mean.
    The constant factor is looked up from a t-statistics table.
    """

    n, v = len(lst), numpy.var(lst, ddof=1)
    c = scipy.stats.t.interval(perc_conf * 1.0 / 100, n - 1)[1]

    return numpy.sqrt(v / n) * c


def aggregate(data, group, fun=(numpy.mean, numpy.nanstd, conf_int)):
    """

    :param data: a pandas dataframe with datat to aggregate
    :param group: a list of grouping columns
    :param fun: a list of aggregating functions
    :return:
    """
    df = data.groupby(group).agg(fun)
    df.columns = ['_'.join(c) for c in df.columns]
    return df.reset_index()


