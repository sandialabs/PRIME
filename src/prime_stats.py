import numpy as np
import h5py
from scipy import stats, mgrid, c_, reshape
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

def distcorr(spl):
    r"""
    Compute distance correlation between random vectors

    Parameters
    ----------

    spl: numpy array [number of samples x number of variables]
         first dimension is the number of samples,
         second dimension is the number of random vectors

    Returns
    -------
       Returns a 2D array of distance correlations between pairs of random vectors;
       only entries 0<=j<i<no. of random vectors are populated
       
    References:
       http://en.wikipedia.org/wiki/Distance_correlation
     """
    nspl,nvars = spl.shape
    if (nspl>5000):
        print('Warning ! This might be a lengthy calculation: nspl=',nspl)
    As=[]
    for i in range(nvars):
        Amat = np.matrix(np.zeros((nspl,nspl)))
        for i1 in range(nspl):
            for j1 in range(nspl):
                Amat[i1,j1] = abs(spl[i1,i]-spl[j1,i])
        # compute means
        Alin = np.array([np.mean(Amat[i1,:]) for i1 in range(nspl)])
        Acol = np.array([np.mean(Amat[:,j1]) for j1 in range(nspl)])
        Amn  = np.mean(Amat)
        #print Amat
        #for i in range(nspl):
        #    print i,Alin[i],Acol[i]
        # subtract/add means (linewise, columnwise, overall)
        Amat = Amat - Alin.reshape(nspl,1)
        Amat = (Amat.T - Acol.reshape(nspl,1)).T
        Amat = Amat+Amn
        #print Amat
        As.append(Amat.copy())
    dCor = np.zeros((nvars,nvars))
    dVarX = [np.sqrt(np.sum(np.multiply(As[i],As[i]))/(nspl*nspl))
       for i in range(nvars)]
    print("Variances:")
    print(dVarX)
    for i in range(1,nvars):
        for j in range(i):
            dCov   = np.sqrt(np.sum(np.multiply(As[i],As[j]))/(nspl*nspl))
            dCor[i,j] = dCov/np.sqrt(dVarX[i]*dVarX[j])
    As=[]
    return dCor

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



