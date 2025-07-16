import sys
import numpy as np
import csv

from   scipy.stats import gamma,lognorm,norm
import matplotlib.pyplot as plt

def normal_logpdf(x,loc,scale):
    return norm._logpdf((x-loc)/scale)-np.log(scale)

def lognorm_pdf(x,s,loc=0,scale=1):
    return lognorm._pdf((x - loc)/scale,s)/scale

def lognorm_cdf(x,s,loc=0,scale=1):
    return lognorm.cdf(x,s,loc=loc,scale=scale)

def gamma_pdf(x,s,loc=0,scale=1):
    return gamma._pdf((x - loc)/scale,s)/scale

def _runningAvgWgts(nDays):
    '''
    Compute average weights
    '''
    disi = np.ones(nDays) / nDays
   
    ka = nDays // 2 + 1

    disb = []
    for i in range(ka,nDays):
        disb.append( np.ones(i) / i )

    return disi,disb

def runningAvg(f,nDays):
    r"""
    Apply nDays running average to the input f

    Parameters
    ----------
    f: numpy array
        array (with daily data for this project) to by filtered 
    nDays: int
        window width for the running average

    Returns
    -------
    favg: numpy array
        filtered data
    """
    disi,disb = _runningAvgWgts(nDays)
    ka = nDays // 2
    npts = f.shape[0]
    favg = np.empty_like(f)
    # interior
    for i in range(ka,favg.shape[0]-ka):
        favg[i] = np.sum(disi*f[i-ka:i+ka+1])
    # boundaries
    for i in range(ka):
        fwidth  = len(disb[i])
        favg[i] = np.sum(disb[i]*f[0:fwidth])
        favg[npts-1-i] = np.sum(disb[i]*f[npts-1:npts-fwidth-1:-1])
    return favg

def test_runningAvg():
    np.random.seed(2020)
    n = 100
    sig = np.random.randn(n)**3 + 3*np.random.randn(n).cumsum()
    plt.plot(sig,label='data')
    sigf = runningAvg(sig,3)
    plt.plot(sigf,label='3')
    sigf = runningAvg(sig,5)
    plt.plot(sigf,label='5')
    sigf = runningAvg(sig,7)
    plt.plot(sigf,label='7')
    sigf = runningAvg(sig,9)
    plt.plot(sigf,label='9')

    plt.legend()
    plt.show()


def prediction_filename(run_setup):
    r"""
    Generate informative name for hdf5 file with prediction data

    Parameters
    ----------
    run_setup: dictionary
        detailed settings for the epidemiological model

    Returns
    -------
    filename: string
        file name ending with a .h5 extension
    """
    fh5 = run_setup["ppopts"]["fpredout"]
    return fh5+".h5"

def output_epicurves(pred,daysPred,newcases,nskip,quantile_list,fileout):
    with open(fileout, mode='w') as output_file:
        csv_writer = csv.writer(output_file,delimiter=',',quotechar='"',quoting=csv.QUOTE_MINIMAL)
        outputData=["#Date"]
        for qk in quantile_list:
            outputData = outputData+["quantile"+'%.3f'%(qk)]
        for j in range(1,pred.shape[0]+1,nskip):
            outputData = outputData+["sample"+str((j-1)//nskip+1)]
        outputData = outputData+["ConfirmedCases"]
        cso = csv_writer.writerow(outputData)
        ndaysData = len(newcases)
        ndaysPred = pred.shape[1]
        for i in range(ndaysPred):
            outputData = [daysPred[i].date()]
            for qk in quantile_list:
                outputData = outputData+["%d"%(np.quantile(pred,qk,axis=0)[i])]
            outputData = outputData+["%d"%(pred[j,i]) for j in range(0,pred.shape[0],nskip)]
            if i < ndaysData:
                outputData = outputData+["%d"%(newcases[i])]
            else:
                outputData = outputData+[-999]
            cso = csv_writer.writerow(outputData)
        output_file.close()

def output_infcurves(infc,datesmean,nskip,quantile_list,fileout):
    with open(fileout, mode='w') as output_file:
        csv_writer = csv.writer(output_file,delimiter=',',quotechar='"',quoting=csv.QUOTE_MINIMAL)
        outputData=["#Date"]
        for qk in quantile_list:
            outputData = outputData+["quantile"+'%.3f'%(qk)]
        for j in range(1,infc.shape[0]+1,nskip):
            outputData = outputData+["sample"+str((j-1)//nskip+1)]
        cso = csv_writer.writerow(outputData)
        for i in range(len(datesmean)):
            outputData = [datesmean[i].date()]
            for qk in quantile_list:
                outputData = outputData+["%d"%(np.quantile(infc,qk,axis=0)[i])]
            outputData = outputData+["%d"%(infc[j,i]) for j in range(0,infc.shape[0],nskip)]
            cso = csv_writer.writerow(outputData)
        output_file.close()

def _linear_error_weight(min_wgt,days):
    '''
    compute linearly increasing weighting
    from  at first data point to 1.0 for most recent data point 
    '''
    ndays = len(days)
    return min_wgt + (1.0 - min_wgt)*np.arange(1,int(ndays)+1) / ndays

def _gaussian_error_weight(min_wgt,tau,days):
    '''
    compute semi-gaussian increasing weighting
    "mean" is at most recent data point. 
    Weight increases from min_wgt to 1
    '''
    day_max = np.max(days)
    return min_wgt + (1.0-min_wgt)*np.exp(-0.5 * ((days-day_max)/tau)**2)

def compute_error_weight(error_info,days):
    r"""
    Compute array with specified weighting for the daily cases data. 
    The weights follow either linear of Gaussian expressions with higher
    weights for recent data and lower weights for older data

    Parameters
    ----------
    error_info: list
        (error_type,min_wgt,[tau]), error type is either 'linear' or 'gaussian',
        min_wgt is the minimum weight and tau is the standard deviation 
        of the exponential term if a Gaussian formulation is chosen. 
    days: int
        lenght of the weights array
    Returns
    -------
    error_weight: numpy array
        array of weights
    """
    error_type = error_info[0]
    assert error_info[1] > 0.0, "error_weight second parameter needs to be positive"
    assert error_info[1] < 1.0, "error_weight second parameter needs to be less than 1.0"

    inv_error_weight = None
    if error_type=="linear":
        inv_error_weight = _linear_error_weight(error_info[1],days)
    elif error_type=="gaussian":
        if len(error_info) < 3:
            sys.exit("Need to specify minimum weight and width for 'gaussian'")
        inv_error_weight = _gaussian_error_weight(error_info[1],error_info[2],days)
    else:
        sys.exit("Only current options for error weighting are 'linear' or 'gaussian'")

    # compute error_weight from reciprocal
    error_weight = 1./inv_error_weight
    return error_weight


