import matplotlib.pyplot as plt
import h5py
import sys
import json
import numpy as np
from prime_stats import distcorr

def main(setupfile):
    r"""
    Computes and saves distance correlations based on samples.
    The distance correlation matrix is saved in "distanceCorr.txt"

    Parameters
    ----------
    setupfile: string
        json file (.json) including run setup information and 
        postprocessing information for an MCMC run. It should 
        specify the name of the file containing the MCMC chain
    """
    run_setup=json.load(open(setupfile))
    print(run_setup)
    
    #-------------------------------------------------------
    # definitions
    fchno = run_setup["regioninfo"]["fchain"]
    
    print("Computing distance correlation data from {}".format(fchno))
    
    # retrieve MCMC chain
    file = h5py.File(fchno, 'r')
    chn  = np.array(file["chain"])
    file.close()
    
    # extract samples from the MCMC chain
    nstart   = run_setup["ppopts"]["nstart"]
    nsamples = run_setup["ppopts"]["nsamples"]
    nskip    = (chn.shape[0]-nstart)//nsamples
    chnPred=chn[nstart::nskip,:]
    
    # compute and save distance correlation matrix
    dcorr = distcorr(chnPred)
    np.savetxt("distanceCorr.txt",dcorr)

if __name__ == '__main__':
    main(sys.argv[1])
