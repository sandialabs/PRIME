import numpy as np
from   prime_utils import lognorm_pdf, lognorm_cdf

def _incubation(time,incubation_median,incubation_sigma):
    """
    Compute the probability density function of the incubation rate,
    currently modeled as log-normal distribution
    """
    time = np.atleast_1d(time)
    vals = np.zeros_like(time)
    I = np.where(time>=0)
    vals[I] = lognorm_pdf(time[I],incubation_sigma,scale=incubation_median)
    return vals

def _incubation_cdf(time,incubation_median,incubation_sigma):
    """
    Compute the cumulative density function of the incubation rate,
    currently modeled as log-normal distribution
    """
    time = np.atleast_1d(time)
    vals = np.zeros_like(time)
    I = np.where(time>=0)
    vals[I] = lognorm_cdf(time[I],incubation_sigma,scale=incubation_median)
    return vals

def incubation_fcn(time,incubation_median,incubation_sigma,is_cdf=False):
    """
    Computes the incubation rate

    Parameters
    ----------
    time: float, list, or numpy array
        instances in time for the evaluation of the incubation rate model
    incubation_median: float
        median of the incubation rate model
    incubation_sigma: float
        standard deviation of the incubation rate model
    is_cdf: boolean (optional, default False)
        select either the CDF of the incubation rate model (True) or 
        its PDF (False)

    Returns
    -------
    vals: numpy array
        incubation rates corresponding to the time values provided as input parameters
    """
    if is_cdf:
        vals = _incubation_cdf(time,incubation_median,incubation_sigma)
    else:
        vals = _incubation(time,incubation_median,incubation_sigma)
    return vals
