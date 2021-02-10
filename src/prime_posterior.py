import numpy as np
from scipy.optimize  import brentq
from scipy.integrate import quad
from scipy.stats     import norm,nbinom
from scipy.stats     import multivariate_normal as mvn
from prime_utils     import normal_logpdf
from prime_model     import modelPred

import datetime
from dateutil import parser

def logpost(state,params):
    """
    Compute log-posterior density values; this function assumes 
    the likelihood is a product of independent Gaussian distributions

    Parameters
    ----------
    state: python list or numpy array
         model parameters
    params: dictionary
         detailed settings for the epidemiological model

    Returns
    -------
    llik: float
        natural logarithm of the likelihood density
    lpri: float
        natural logarithm of the prior density
    """
    # parameters
    new_cases         = params['new_cases']
    prior_types       = params['prior_types']
    prior_info        = params['prior_info']
    model_type        = params['model_type']
    error_model_type  = params['error_model_type']
    error_weight      = params['error_weight']

    assert(params['days_extra']==0)
    # evaluate the model
    people_with_symptoms_cdf = modelPred(state,params,is_cdf=True)
    people_with_symptoms = np.zeros(people_with_symptoms_cdf.shape[0])
    for i in range(1,people_with_symptoms.shape[0]):
        people_with_symptoms[i] = people_with_symptoms_cdf[i]-people_with_symptoms_cdf[i-1]
    # print(people_with_symptoms_cdf,people_with_symptoms)
    # quit()
    # log-likelihood
    ndays = params['days_since_day0'].shape[0]
    llik  = 0.0
    if error_model_type == "add":
        # additive error
        err = np.exp(state[-1])
    elif error_model_type == "addMult":
        # additive & multiplicative error
        err = np.exp(state[-2])+np.exp(state[-1])*people_with_symptoms

    # apply weighting to error terms if specified
    if error_weight is not None:
        err *= error_weight

    # kchowdh: vectorize log norm pdf
    npws = (people_with_symptoms - new_cases)/err

    # apply weighting to error terms if specified
    llik = np.sum(norm._logpdf(npws) - np.log(err))

    # log-prior
    lpri = 0.0
    for i in range(state.shape[0]):
        if prior_types[i]=='g':
            log_pdf_vals = normal_logpdf(state[i],loc=prior_info[i][0],scale=prior_info[i][1])
            lpri = lpri+log_pdf_vals
    return [llik,lpri]

def logpost_negb(state,params):
    """
    Compute log-posterior density values; this function assumes 
    the likelihood is a product of negative-binomial distributions

    Parameters
    ----------
    state: python list or numpy array
         model parameters
    params: dictionary
         detailed settings for the epidemiological model

    Returns
    -------
    llik: float
        natural logarithm of the likelihood density
    lpri: float
        natural logarithm of the prior density
    """
    # parameters
    new_cases    = params['new_cases']
    prior_types  = params['prior_types']
    prior_info   = params['prior_info']
    model_type   = params['model_type']
    error_weight = params['error_weight']
    assert(params['days_extra']==0)
    # compute cases
    # people_with_symptoms = modelPred(state,params)
    people_with_symptoms_cdf = modelPred(state,params,is_cdf=True)
    people_with_symptoms = np.zeros(people_with_symptoms_cdf.shape[0])
    for i in range(1,people_with_symptoms.shape[0]):
        people_with_symptoms[i] = people_with_symptoms_cdf[i]-people_with_symptoms_cdf[i-1]
    # log-likelihood
    alpha_ind = 4
    if model_type == "twoWave":
        alpha_ind = 8
    elif model_type == "threeWave":
        alpha_ind = 12
    alpha = np.exp(state[alpha_ind]) 
    prob = alpha/(alpha+people_with_symptoms)
    llkarray=np.array([np.log(1e-10+nbinom._pmf(obs, n=alpha, p=p)) for obs,p in zip(new_cases,prob)])
     # apply weighting to error terms if specified
    if error_weight is not None:
        llkarray += np.log(error_weight)
    
    llik = np.sum(llkarray[1:])
    # log-prior
    lpri = 0.0
    for i in range(state.shape[0]):
        if prior_types[i]=='g':
            log_pdf_vals = normal_logpdf(state[i],loc=prior_info[i][0],scale=prior_info[i][1])
            lpri = lpri+log_pdf_vals
    return [llik,lpri]

def logpost_poisson(state,params):
    """
    Compute log-posterior density values; this function assumes 
    the likelihood is a product of poisson distributions

    Parameters
    ----------
    state: python list or numpy array
         model parameters
    params: dictionary
         detailed settings for the epidemiological model

    Returns
    -------
    llik: float
        natural logarithm of the likelihood density
    lpri: float
        natural logarithm of the prior density
    """
    # parameters
    new_cases    = params['new_cases']
    prior_types  = params['prior_types']
    prior_info   = params['prior_info']
    error_weight = params['error_weight']
    assert(params['days_extra']==0)
    # compute cases
    # people_with_symptoms = modelPred(state,params)
    people_with_symptoms_cdf = modelPred(state,params,is_cdf=True)
    people_with_symptoms = np.zeros(people_with_symptoms_cdf.shape[0])
    for i in range(1,people_with_symptoms.shape[0]):
        people_with_symptoms[i] = people_with_symptoms_cdf[i]-people_with_symptoms_cdf[i-1]
    # log-likelihood
    # alpha = np.exp(state[4])
    llkarray=np.array([-lbd+k*np.log(lbd+1.e-4) for k,lbd in zip(new_cases,people_with_symptoms)])
    # apply weighting to error terms if specified
    if error_weight is not None:
        llkarray += np.log(error_weight)
    
    llik = np.sum(llkarray[1:])-params["sumLogK"]
    # log-prior
    lpri = 0.0
    for i in range(state.shape[0]):
        if prior_types[i]=='g':
            log_pdf_vals = normal_logpdf(state[i],loc=prior_info[i][0],scale=prior_info[i][1])
            lpri = lpri+log_pdf_vals
    return [llik,lpri]


