import numpy as np
from prime_infection  import infection_rate
from prime_incubation import incubation_fcn

inc_median_lo, inc_median_hi = 1.504, 1.755
inc_sigma_lo,  inc_sigma_hi  = 0.271, 0.542

# Gauss-Legendre quadrature
nquad = 30
gauss_x,gauss_w = np.polynomial.legendre.leggauss(nquad) # x in [-1,1]

# Convolution option
nconv = 250

def _modelPred_oneWave(state, params, is_cdf=False):
    """
    model predictions
        - it computes predictions for a list of days given (presumably when we have data)
          and adds a number of extra days for forward predictions
        - when used to compute log-likelihood the number of extra days should be zero
          most of the time; this is left unchecked here... just FYI
        
    Convolution
        - User can now choose to use a convolution intead of the
        default quadrature. 
        - Just add "useconv" = 1 to the params dictionary
    """
    # parameters
    days_since_day0 = params['days_since_day0'][0] # need to fix this at some point
    inftype         = params['inftype']
    days_extra      = params['days_extra']

    # fixed or uncertain incubation
    if "incubation_type" in params:
        if params["incubation_type"] == "stochastic":
            incubation_median = np.exp(inc_median_lo + (inc_median_hi-inc_median_lo)* np.random.uniform())
            incubation_sigma  = inc_sigma_lo +( inc_sigma_hi - inc_sigma_lo ) * np.random.uniform()
        elif params["incubation_type"] == "fixed" :
            incubation_median = params['incubation_median']
            incubation_sigma  = params['incubation_sigma']
        else:
            print("_modelPred_oneWave(): Unknown option for incubation_type!")
            quit()
    else:
        incubation_median = params['incubation_median']
        incubation_sigma  = params['incubation_sigma']
    
    # add convolution option
    if "useconv" in params:
        useconv = params["useconv"]
    else:
        useconv = 0
    #----------------------------------------------------
    t0, N, qshape, qscale = state[:4]

    #----------------------------------------------------
    ndays = days_since_day0.shape[0]
    days_since_day0 = days_since_day0 - t0
    people_with_symptoms = np.zeros(ndays+days_extra)
  
    # infection/incubation function
    f = lambda tau: infection_rate(tau, qshape, qscale, inftype) 
    g = lambda tau: incubation_fcn(tau, incubation_median, incubation_sigma, is_cdf=is_cdf)

    if useconv == 0:

        for i in range(days_since_day0.shape[0]):
            dayNo   = days_since_day0[i]
            xi = (gauss_x+1)/2 * dayNo
            wi = gauss_w/2     * dayNo
            inf_inc = f(xi) * g(xi)
            people_with_symptoms[i] = inf_inc.dot(wi)

        for i in range(days_extra):
            dayNo   = days_since_day0[-1]+i+1
            xi=(gauss_x+1)/2*dayNo
            wi=gauss_w/2*dayNo
            inf_inc = f(xi)*g(xi)
            people_with_symptoms[ndays+i] = inf_inc.dot(wi)

    elif useconv == 1:

        iday = days_since_day0[0]  # first day
        fday = days_since_day0[-1]+days_extra # final day

        # alternative using interpolation
        ti    = np.linspace(0, fday, nconv)
        dti   = ti[1] - ti[0]
        fi    = f(ti)
        fi[0] = 0.0   # set to zero since f(0) is undefined
        gi    = g(ti)
        pconv    = dti * np.convolve(fi,gi)[:len(ti)]
        mod_days = days_since_day0.copy()
        mod_days = np.append(mod_days, [days_since_day0[-1]+i+1 for i in range(days_extra)])
        mod_days[mod_days <= 0] = 0 
        people_with_symptoms = np.interp(mod_days,ti,pconv)

    else:

        print("Unknown flag value for params[\"useconv\"]:",useconv)
        quit()

    return people_with_symptoms * N # deleted 1.e6 since now we work with normalized counts

def _modelPred_multiWave(state, params, n_waves=1, n_regions=1, is_cdf=False):
    '''
    Multi-wave model predictions
        - it computes predictions for a list of days given (presumably when we have data)
          and adds a number of extra days for forward predictions
        - this version assumes two infection waves
        - can output pdf or cdf data by specifying the optional "is_cdf" input
        - when used to compute log-likelihood the number of extra days should be zero
          most of the time; this is left unchecked here... just FYI
    '''

    # compute cases for each wave
    # assert n_waves == 1 # for now limit multi-region at 1 wave
 
    Ncases = [_modelPred_oneWave(state[4*j:4*j+4],params,is_cdf=is_cdf) for j in range(n_regions)]

    if n_waves > 1:
        for i in range(1,n_waves):
            for j in range(n_regions):
                ist = i * 4 * n_regions
                state_i = state[ist+j*4:ist+j*4+4]
                state_i[0] = state_i[0]+state[4*j] # shift by the first t0
                Ncases[j] = Ncases[j] +_modelPred_oneWave(state_i,params,is_cdf=is_cdf)

    return Ncases

def modelPred(state, params, n_regions=1, is_cdf=False):
    r"""
    Evaluates the PRIME model for a set of model parameters; specific model settings 
    (e.g. date range, other control knobs, etc) are specified via the \"params\" dictionary

    Parameters
    ----------
    state: python list or numpy array
        model parameters
    params: dictionary
        detailed settings for the epidemiological model
    is_cdf: boolean (optional, default False)
        estimate the epidemiological curve based on the CDF of the incubation
        model (True) or via the formulation that employs the PDF of the
        icubation model (False)

    Returns
    -------
    Ncases: numpy array
        daily counts for people turning symptomatic
    """
    assert "num_waves" in params
    return _modelPred_multiWave(state, params, n_regions=n_regions, n_waves=params["num_waves"], is_cdf=is_cdf)
