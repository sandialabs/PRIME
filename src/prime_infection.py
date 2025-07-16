import numpy    as np
from   dateutil import parser
import datetime

from prime_utils  import lognorm_pdf, gamma_pdf

def infection_rate(time,qshape,qscale,inftype):
    """
    Infection rate (gamma or log-normal distribution)

    Parameters
    ----------
    time: float, list, or numpy array
        instances in time for the evaluation of the infection_rate model
    qshape: float
        shape parameter
    qscale: float
        scale parameter
    inftype: string
        infection rate type (\"gamma\" for Gamma distribution, otherwise
        the Log-normal distribution)

    Returns
    -------
    vals: numpy array
        infection rates corresponding to the time values provided as input parameters
    """
    time = np.atleast_1d(time)
    vals = np.zeros_like(time)
    I = np.where(time>0)
    if inftype=="gamma":
        vals[I] = gamma_pdf(time[I],qshape,scale=qscale)
    else:
        vals[I] = lognorm_pdf(time[I],qshape,scale=qscale)
    return vals
        
def _infection_multiWave(state,params,n_waves=1,regionID=0):
    """
    Compute infection curve for multi-wave epidemics
      - this function is currently used by the post-processing script to push-forward
        the posterior into a set of infection curves that are consistent with the 
        observed cases
    """
    assert n_waves>0

    day0  = params['day0']
    ndays = params['ndays']

    # compute cases
    t01, N1, qshape1, qscale1 = state[regionID*4:regionID*4+4]
    dates = np.array([parser.parse(day0)+datetime.timedelta(days=int(t01)+i) for i in range(ndays)])
    infections = N1 * 1e6 * infection_rate(np.arange(1.0*ndays),qshape1,qscale1,params['inftype'])

    # infections due to multiple waves
    for i in range(1,n_waves):
        dt1_i, N_i, qshape_i, qscale_i = state[4*i:4*i+4]
        ndays_wave = ndays-int(dt1_i)
        infections_wave = N_i * 1e6 * infection_rate(np.arange(1.0*ndays_wave),qshape_i,qscale_i,params['inftype'])
        # print(len(infections_wave),int(dt1_i),len(infections[int(dt1_i):]),len(infections))
        if int(dt1_i)>=0:
            infections[int(dt1_i):] = infections[int(dt1_i):]+infections_wave
        else:
            infections = infections+infections_wave[-int(dt1_i):]            

    # convert dates to timestamp 
    dates = np.array([dates[i].timestamp() for i in range(ndays)])

    return [dates,infections]

def infection(state,params,regionID=0):
    """
    Compute infection curve for multi-wave epidemics
      - this function is currently used by the post-processing script to push-forward
        the posterior into a set of infection curves that are consistent with the 
        observed cases
        
    Parameters
    ----------
    state: python list or numpy array
         model parameters
    params: dictionary
         detailed settings for the epidemiological model

    Returns
    -------
    dates: numpy array
        list of dates for which the infection rates were computed
    infectons: numpy array
        infection rate values corresponding to the list of dates
    """
    # model_type  = params['model_type']
    # wavecounts={"oneWave":1,"twoWave":2,"threeWave":3}
    # return _infection_multiWave(state,params,n_waves=wavecounts[model_type])
    num_waves  = params['num_waves']
    return _infection_multiWave(state,params,n_waves=num_waves,regionID=regionID)
