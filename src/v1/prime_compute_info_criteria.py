import sys
import numpy as np
import json

from prime_utils import runningAvg
from prime_stats import computeAICandBIC,computeCRPS

def main(setupfile):
    r"""
    This script postprocesses data from PRIME to compute statistical information including:
    - AIC: Akaike Information Criterion
    - BIC: Bayesian Information Criterion
    - CPRS: Continuous Rank Probability Score
    Results are saved in "info_criteria.txt"

    Parameters
    ----------
    setupfile: string
        json file (.json) including run setup information and 
        postprocessing information for an MCMC run. It should 
        specify the name of the file containing the MCMC chain
    """
    #-------------------------------------------------------
    run_setup=json.load(open(setupfile))
    print(run_setup)
    
    #-------------------------------------------------------
    # compute AIC and BIC
    
    AIC,BIC = computeAICandBIC(run_setup)
    print("AIC={}".format(AIC))
    print("BIC={}".format(BIC))
    
    #-------------------------------------------------------
    # compute CRPS
    
    CRPS = computeCRPS(run_setup)
    print("CRPS={}".format(CRPS))
    
    # save to file
    np.savetxt("info_criteria.txt",np.array([AIC,BIC,CRPS]))

if __name__ == '__main__':
    main(sys.argv[1])
