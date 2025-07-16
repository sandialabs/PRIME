#!/usr/bin/env /Users/csafta/miniconda3/bin/python
import sys
import numpy as np
from   dateutil import parser
import json

from prime_mcmc       import ammcmc, save_mcmc_chain
from prime_posterior  import logpost, logpost_negb, logpost_poisson
from prime_utils      import runningAvg, compute_error_weight

def get_opts(setupfile,verb=False,return_run_setup=False):

    run_setup=json.load(open(setupfile))
    if verb:
        print("=====================================================")
        print(run_setup)
        print("=====================================================")

    run_opts = dict()

    #-daily counts
    run_opts["count_data"] = run_setup["regioninfo"]["count_data"]
    run_opts["population_data"] = run_setup["regioninfo"]["population_data"]
    if "running_avg_obs" in run_setup["regioninfo"]:
        run_opts["running_avg_obs"] = run_setup["regioninfo"]["running_avg_obs"]

    run_opts["region_tag"] = run_setup["regioninfo"]["region_tag"]
    run_opts["day0"] = run_setup["regioninfo"]["day0"]

    #------------------------------------------------------------------
    #-incubation model
    assert "num_waves" in run_setup["modelopts"]
    run_opts["num_waves"] = run_setup["modelopts"]["num_waves"]
    run_opts["useconv"] = run_setup["modelopts"]["useconv"]
    run_opts["inc_median"] = run_setup["modelopts"]["incubation_median"]
    run_opts["inc_sigma"] = run_setup["modelopts"]["incubation_sigma"]

    if "incubation_model" in run_setup["modelopts"]:
        run_opts["inc_model"] = run_setup["modelopts"]["incubation_model"]
    else:
        run_opts["inc_model"] = "lognormal"

    if "incubation_type" in run_setup["modelopts"]:
        run_opts["inc_type"] = run_setup["modelopts"]["incubation_type"]
    else:
        run_opts["inc_type"] = "deterministic"
    
    #------------------------------------------------------------------
    #-mcmc model parameters
    run_opts["mcmc_log"] = run_setup["mcmcopts"]["logfile"]
    run_opts["mcmc_nsteps"] = run_setup["mcmcopts"]["nsteps"]
    run_opts["mcmc_nfinal"] = run_setup["mcmcopts"]["nfinal"]
    run_opts["mcmc_gamma"] = run_setup["mcmcopts"]["gamma"]

    run_opts["inicov"] = np.array(run_setup["mcmcopts"]["cvini"])
    run_opts["inistate"] = run_setup["mcmcopts"]["cini"]
    if len(run_opts["inicov"].shape) == 1:
        run_opts['inicov'] = np.diag(run_opts["inicov"])

    run_opts["spllo"] = np.array(run_setup["mcmcopts"]["spllo"])
    run_opts["splhi"] = np.array(run_setup["mcmcopts"]["splhi"])

    #------------------------------------------------------------------
    #-bayes framework
    run_opts["lpf_type"] = run_setup["bayesmod"]["lpf_type"]
    run_opts["error_model_type"] = run_setup["bayesmod"]["error_model_type"]
    run_opts["prior_types"] = run_setup["bayesmod"]["prior_types"]
    run_opts["prior_info"] = run_setup["bayesmod"]["prior_info"]
          
    #------------------------------------------------------------------
    run_opts["days_extra"] = run_setup["ppopts"]["days_extra"]

    if return_run_setup:
        return run_opts, run_setup
    else:
        return run_opts

def get_counts(run_opts,return_raw_data=False):
    """
    Get counts from raw files
    """
    # extract data from raw data
    days_since_day0 = []
    daily_counts = []
    rawdata_all = []
    for ireg, region in enumerate(run_opts["count_data"]):
        rawdata = np.loadtxt(region,delimiter=",",dtype=str)
        rawdata_all.append(rawdata)
        ndays = rawdata.shape[0]
        days_since_day0.append(np.array([(parser.parse(rawdata[i,0])-parser.parse(run_opts["day0"])).days \
                                    for i in range(ndays)]))

        daily_counts.append(np.array([float(rawdata[i,1]) for i in range(rawdata.shape[0])]))
        # scale daily counts
        daily_counts[-1] = daily_counts[-1]/(run_opts["population_data"][ireg] * 1.e6)
        # run averages
        if "running_avg_obs" in run_opts:
            daily_counts[-1] = runningAvg(daily_counts[-1], run_opts["running_avg_obs"])
            print("Taking {}-day running average of observations for {}".format(run_opts["running_avg_obs"],run_opts["region_tag"][ireg]))

    if return_raw_data:
        return days_since_day0, daily_counts, rawdata_all
    else:
        return days_since_day0, daily_counts


def main(setupfile):
    r"""
    Driver script to run MCMC for parameter inference for a multi-wave 
    epidemic model. Currently limited to up to three infection curves.

    To run this script:

    python <path-to-this-directory>/prime_run.py <name-of-json-input-file>

    Parameters
    ----------
    setupfile: string
        json format input file with information on observations data, filtering options,
        MCMC options, and postprocessing options. See "setup_template.json" for a detailed
        example
    """
    #----------- --------------------------------------------
    setupfile=sys.argv[1]
    run_opts = get_opts(setupfile)

    print(run_opts)
    print("=====================================================")
    
    # #-------------------------------------------------------
    # # definitions
    # fdata = run_setup["regioninfo"]["regionname"]+".dat"
    # day0  = run_setup["regioninfo"]["day0"]
    
    #-------------------------------------------------------
    # echo some settings
    print("Running inference with %d waves"%(run_opts["num_waves"]))   
    print("Error model %s"%(run_opts["error_model_type"]))
    assert run_opts["error_model_type"] in ["add","addMult"]

    #-------------------------------------------------------
    # get counts
    days_since_day0, daily_counts = get_counts(run_opts)
    
    #-------------------------------------------------------
    # mcmc
    opts = {"nsteps": run_opts["mcmc_nsteps"], "nfinal": run_opts["mcmc_nfinal"],"gamma": run_opts["mcmc_gamma"],
            "inicov": np.array(run_opts["inicov"]),"inistate": np.array(run_opts["inistate"]),
            "spllo": np.array(run_opts["spllo"]),"splhi": np.array(run_opts["splhi"]),
            "logfile": run_opts["mcmc_log"],"burnsc":5,
            "nburn":1000,"nadapt":100,"coveps":1.e-10,"ofreq":5000,"tmpchn":"tmpchn"
            }

    modelinfo={"num_waves":        run_opts["num_waves"],
               "error_model_type": run_opts["error_model_type"],
               "days_since_day0":  days_since_day0,
               "daily_counts":     daily_counts,
               "incubation_model": run_opts["inc_model"],
               "incubation_median":run_opts["inc_median"],
               "incubation_sigma": run_opts["inc_sigma"], 
               "incubation_type":  run_opts["inc_type"], 
               "inftype":          "gamma",
               "useconv":          run_opts["useconv"],
               "days_extra":       0,
               "prior_types":run_opts["prior_types"],"prior_info": run_opts["prior_info"]}
    
    # Convolution vs Quadrature:
    #   -The user can choose to use a fft convolution instead of 
    #    quadrature to perform the integration of Y(t)
    #   -default is set to zero if the user defines nothing
    #   -To set, add "useconv":1 to the mcmcopts in the *json file 
    if modelinfo["useconv"] == 1:
        print("Using FFT convolution instead of quadrature")
    
    # choose log-posterior function
    logpost_types={"gaussian":logpost,"negative_binomial":logpost_negb,"poisson":logpost_poisson}
    lpf = run_opts["lpf_type"]
    if lpf == "poisson":
        modelinfo["sumLogK"] = sum([sum([np.log(i) for i in range(1,int(k)+1)]) for k in daily_counts if k>0])
    
    # run MCMC
    sol=ammcmc(opts,logpost_types[lpf],modelinfo)
    save_mcmc_chain("".join(run_opts["region_tag"])+"_mcmc.h5",sol)

if __name__ == '__main__':
    main(sys.argv[1])

