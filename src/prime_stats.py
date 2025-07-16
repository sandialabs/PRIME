import numpy as np
import h5py
from scipy import stats, mgrid, c_, reshape
from scipy.spatial.distance import cdist
from dateutil    import parser

from prime_utils import runningAvg, prediction_filename


def getKDE(spl,nskip=0,nthin=1,npts=100,bwfac=1.0):
    r"""
    Compute 1D and 2D marginal PDFs via Kernel Density Estimate

    Parameters
    ----------
    spl: numpy array
         MCMC chain [number of samples x number of parameters]
    nskip: int
           number of initial samples to skip when sampling the MCMC chain
    nthin: int
           use every 'nthin' samples
    npts: int
          number of grid points
    bwfac: double
          bandwidth factor

    Returns
    -------
    dict: dictionary with results
         'x1D': list of numpy arrays with grids for the 1D PDFs;
         'p1D': list of numpy arrays with 1D PDFs;
         'x2D': list of numpy arrays of x-axis grids for the 2D PDFs;
         'y2D': list of numpy arrays of y-axis grids for the 2D PDFs;
         'p2D': list of numpy arrays containing 2D PDFs
    """
    spl=spl[nskip::nthin,:]
    # Compute 1D pdf's
    x1D=[];
    p1D=[];
    for i in range(spl.shape[1]):
        print(" - Evaluate 1D pdf for variable",i+1)
        spls = spl[:,i]
        x1D.append(np.linspace(spls.min(),spls.max(),npts));
        kern=stats.kde.gaussian_kde(spls);
        p1D.append(kern(x1D[i]));
    # Compute 2D joint pdf's
    npts = 1j*npts
    x2D = []
    y2D = []
    p2D = []
    for j in range(1,spl.shape[1]):
        for i in range(j):
            print(" - Evaluate 2D pdf for variables",i+1,j+1)
            kern=stats.kde.gaussian_kde(c_[spl[:,i],spl[:,j]].T)
            kern.set_bandwidth(bw_method=kern.factor*bwfac)
            xmin = spl[:,i].min(); xmax = spl[:,i].max()
            ymin = spl[:,j].min(); ymax = spl[:,j].max()
            x,y = mgrid[xmin:xmax:npts,ymin:ymax:npts]
            z   = reshape(kern(c_[x.ravel(), y.ravel()].T).T, x.T.shape)
            x2D.append(x);
            y2D.append(y);
            p2D.append(z);
    return {'x1D':x1D,'p1D':p1D,'x2D':x2D,'y2D':y2D,'p2D':p2D}

NWARN=10000 # if using more than this (approximately), things might become very costly

def _compute_centered_matrix(spl):
    """
    centers a matrix
    """
    num_samples = spl.shape[0]
    if len(spl.shape) == 1:
        Amat = np.abs(np.subtract.outer(spl,spl))
    elif len(spl.shape) == 2:
        Amat = cdist(spl, spl, metric='euclidean')
    else:
        print('_compute_centered_matrix works with order 1 and 2 tensors')
        quit()

    # compute means
    Alin = np.mean(Amat, axis=1)
    Acol = np.mean(Amat, axis=0)
    Amn = np.mean(Amat)
    # subtract/add means (linewise, columnwise, overall)
    Amat = Amat - Alin.reshape(num_samples,1)
    Amat = (Amat.T - Acol.reshape(num_samples,1)).T
    Amat = Amat+Amn
    return Amat

def distcorr(spl,rv_sizes=None):
    """
    Compute distance correlation between random vectors,
    possibly of different dimensionalities

    Args:
        spl - 2D array of samples, each row is a sample
        rv_sizes - array/list of rv dimensions (assumed all equal to 1 if not provided)

    Output:
       Returns a 2D array of distance correlations between pairs of random vectors;
       only entries 0<=j<i<no. of random vectors are populated

    References:
       * http://en.wikipedia.org/wiki/Distance_correlation
       * Szekely & Rizzo, "Brownian Distance Covariance", The Annals of Applied Statistics
         2009, Vol. 3, No. 4, 1236–1265, DOI: 10.1214/09-AOAS312
    """
    num_samples = spl.shape[0]
    if (num_samples > NWARN):
        print('Warning ! This might be a lengthy calculation: num_samples='+str(num_samples))

    if rv_sizes is not None:
        assert sum(rv_sizes)==spl.shape[1]
        num_rvs = len(rv_sizes)
        rv_idx = np.insert(np.cumsum(rv_sizes),0,0)
    else:
        num_rvs = spl.shape[1]
        rv_idx = np.arange(num_rvs+1)

    dCor = np.zeros((num_rvs, num_rvs))
    As = [_compute_centered_matrix(spl[:,rv_idx[0]:rv_idx[1]])]
    dSigX = [np.sqrt(np.sum(As[0]**2))/num_samples]
    for i in range(1,num_rvs):

        if rv_idx[i+1]-rv_idx[i]>1:
            As.append(_compute_centered_matrix(spl[:,rv_idx[i]:rv_idx[i+1]]))
        else:
            As.append(_compute_centered_matrix(spl[:,rv_idx[i]]))
        dSigX.append(np.sqrt(np.sum(As[-1]**2))/num_samples)

        if i%10 == 0:
            print('{},'.format(i),end='',flush=True)

        for j in range(i):
            dCov = np.sum(np.multiply(As[i], As[j]))/(num_samples * num_samples)
            dCor[i,j] = dCov/(dSigX[i] * dSigX[j])

    if num_rvs>9:
        print('\n')

    As.clear()
    return np.sqrt(dCor)

def computeAICandBIC(run_setup,verbose=0):
    """
    Compute Akaike Information Criterion (AIC) and Bayesian Information Criterion (BIC)

    Parameters
    ----------
    run_setup: dictionary with run settings; see the Examples section in the manual

    Returns
    -------
    AIC: float
    BIC: float
    """
    #-------------------------------------------------------
    # definitions
    fdata = run_setup["regioninfo"]["regionname"]+".dat"
    fchno = run_setup["regioninfo"]["fchain"]
    
    #-------------------------------------------------------
    # retrieve log likelihoods
    file = h5py.File(fchno, 'r')
    loglik = np.array(file["minfo"][:,1])
    file.close()
    
    if verbose>0:
        print(loglik.shape)
    
    #-------------------------------------------------------
    # compute AIC and BIC
    
    # Read in initial parameter guess to obtain number of parameters
    Nparam = len(run_setup["mcmcopts"]["cini"])
    
    # compute number of data points used (number of days)
    rawdata = np.loadtxt(fdata,dtype=str)
    ndays = rawdata.shape[0]

    AIC = 2 * Nparam - 2 * np.max(loglik)
    BIC = 2 * Nparam * np.log(ndays) - 2 * np.max(loglik)

    return AIC, BIC

def computeCRPS(run_setup):
    """
    Compute Continuous Rank Predictive Score (CRPS)

    Parameters
    ----------
    run_setup: dictionary with run settings; see the Examples section in the manual

    Returns
    -------
    CRPS: float
    """
    #-------------------------------------------------------
    # definitions
    fdata = run_setup["regioninfo"]["regionname"]+".dat"
    fchno = run_setup["regioninfo"]["fchain"]
    day0  = run_setup["regioninfo"]["day0"]
 
    #-------------------------------------------------------
    # extract data from raw data
    rawdata = np.loadtxt(fdata,dtype=str)
    ndays_data = rawdata.shape[0]
    days_since_day0 = np.array([(parser.parse(rawdata[i,0])-parser.parse(day0)).days for i in range(ndays_data)])
    new_cases       = runningAvg(np.array([float(rawdata[i,1]) for i in range(ndays_data)]),7)
    

    #------------------------------------------------------
    # read model predictions
    filename = prediction_filename(run_setup)
    file   = h5py.File(filename, 'r')
    pred   = np.array(file["predictions"])
    file.close()
    
    # only need predictions for days on which data is available
    pred = pred[:,:ndays_data]

    #------------------------------------------------------
    # compute integrals on each day
    CRPS_daily = np.zeros(ndays_data)

    for i in range(ndays_data):
        pred_i = pred[:,i]

        # compute empirical pdf and cdf
        y = np.sort(pred_i)
        ecdf = np.array(range(1,len(pred_i)+1)) / float(len(pred_i))

        # compute heaviside function
        h_i = np.ones(ecdf.shape) * (y > new_cases[i])
     
        # compute integrand 
        integrand = (ecdf - h_i) ** 2

        # average for trapezoid rule integration
        integrand = (integrand[1:] + integrand[:-1]) / 2

        # compute inegral for day i
        CRPS_daily[i] = np.sum( (y[1:] - y[:-1]) * integrand) 

    
    CRPS = np.sum(CRPS_daily) / ndays_data
    return CRPS



