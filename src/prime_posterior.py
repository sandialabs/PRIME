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
    daily_counts      = params['daily_counts']
    prior_types       = params['prior_types']
    prior_info        = params['prior_info']
    error_model_type  = params['error_model_type']
    # error_weight      = params['error_weight']

    assert(params['days_extra']==0)

    # evaluate the model
    num_regions = len(daily_counts)
    people_with_symptoms_cdf = modelPred(state, params, n_regions=num_regions, is_cdf=True)
    reg_people_with_symptoms = []
    for i in range(len(daily_counts)):
        people_with_symptoms = np.ediff1d(people_with_symptoms_cdf[i], to_begin=0)
        reg_people_with_symptoms.append(people_with_symptoms)

    # debug text
    # print(people_with_symptoms_cdf,people_with_symptoms)
    # quit()

    # log-likelihood for several regions
    ndays = params['days_since_day0'][0].shape[0]
    llik = 0.0

    cvm0 = 0.0
    if num_regions > 1:
        if num_regions==4:
            tauphi2 = np.exp(state[-6])
            lbdphi = state[-5]
        else:
            tauphi2 = np.exp(state[-4])
            # tauphi2 = tauphi2/(6.7e5)**2
            lbdphi = state[-3]
        lbdphi2 = lbdphi * lbdphi
        if num_regions == 2:
            cvm0 = tauphi2/(1.-lbdphi**2) * np.array([[1,lbdphi],[lbdphi,1]])
        elif num_regions == 3:
            row0 = [1,            0.5*lbdphi,    0.5*lbdphi]
            row1 = [0.5*lbdphi, 1-0.5*lbdphi2,   0.5*lbdphi2]
            row2 = [0.5*lbdphi,   0.5*lbdphi2, 1-0.5*lbdphi2]
            cvm0  = tauphi2/(1.-lbdphi2) * np.array([row0,row1,row2])
        elif num_regions == 4:
            row0 = [1,            0.5*lbdphi,    0.5*lbdphi,  0.0]
            row1 = [0.5*lbdphi, 1-0.5*lbdphi2,   0.5*lbdphi2, 0.0]
            row2 = [0.5*lbdphi,   0.5*lbdphi2, 1-0.5*lbdphi2, 0.0]
            row3 = [0.0, 0.0, 0.0, 0.0]
            cvm0  = tauphi2/(1.-lbdphi2) * np.array([row0,row1,row2,row3])
        else:
            print("Number of regions {} not modeled yet -> exit".format(num_regions))
            sys.exit

    # print(num_regions)
    if num_regions==4:
        siga,   sigm   = np.exp(state[-4:-2])    
        siga_1, sigm_1 = np.exp(state[-2:])    
    else:
        siga, sigm = np.exp(state[-2:])   
                 
    for i in range(ndays):
        pws = [reg_people_with_symptoms[k][i] for k in range(num_regions)]
        if num_regions==4:
            errM = [(siga+sigm*p)**2 for p in pws[:-1]]+[(siga_1+sigm_1*pws[-1])**2] # error model
        else:
            errM = [(siga+sigm*p)**2 for p in pws] # error model
        cvm = cvm0 + np.diag(errM)
        # print(i,cvm,np.linalg.inv(cvm))
        llik += mvn.logpdf(pws, mean=[daily_counts[k][i] for k in range(num_regions)], cov=cvm)
 
    # log-prior
    lpri = 0.0
    # if num_regions > 1 : 
    #     tauphi2 = tauphi2*(6.7e5)**2

    for i in range(state.shape[0]):

        if prior_types[i]=='g':
            log_pdf_vals = normal_logpdf(state[i],loc=prior_info[i][0],scale=prior_info[i][1])
            lpri = lpri+log_pdf_vals


        # gamma pdf prior for tau_phi^2 with shape:prior_info[i][0] and scale:prior_info[i][1]
        if num_regions > 1 and prior_types[i]=='lgm':
            log_pdf_vals = -tauphi2/prior_info[i][1]+prior_info[i][0]*np.log(tauphi2)
            #log_pdf_vals = - prior_info[i][0]*np.log(tauphi2) - prior_info[i][1]/tauphi2
            lpri = lpri + log_pdf_vals

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
    daily_counts = params['daily_counts']
    prior_types  = params['prior_types']
    prior_info   = params['prior_info']
    num_waves    = params['num_waves']
    #assert(params['days_extra']==0)
    # compute cases
    # people_with_symptoms = modelPred(state,params)
    people_with_symptoms_cdf = modelPred(state,params,is_cdf=True)
    people_with_symptoms = np.zeros(people_with_symptoms_cdf.shape[0])
    for i in range(1,people_with_symptoms.shape[0]):
        people_with_symptoms[i] = people_with_symptoms_cdf[i]-people_with_symptoms_cdf[i-1]

    # log-likelihood - negative binomial formulation
    alpha_ind = 4 * num_waves
    alpha = np.exp(state[alpha_ind]) 
    prob = alpha/(alpha+people_with_symptoms)
    llkarray=np.array([np.log(1e-10 + nbinom._pmf(obs, n=alpha, p=p)) for obs,p in zip(daily_counts,prob)])
    
    llik = np.sum(llkarray[1:])
    # log-prior
    lpri = 0.0
    for i in range(state.shape[0]):
        if prior_types[i]=='g':
            log_pdf_vals = normal_logpdf(state[i],loc=prior_info[i][0],scale=prior_info[i][1])
            lpri = lpri + log_pdf_vals
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
    daily_counts = params['daily_counts']
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
    llkarray=np.array([-lbd+k*np.log(lbd+1.e-4) for k,lbd in zip(daily_counts,people_with_symptoms)])

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


