import sys
import numpy as np
from   dateutil import parser
import json

from prime_mcmc       import ammcmc
from prime_posterior  import logpost, logpost_negb, logpost_poisson
from prime_utils      import runningAvg, compute_error_weight


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
    #-------------------------------------------------------
    run_setup=json.load(open(sys.argv[1]))
    print("=====================================================")
    print(run_setup)
    print("=====================================================")
    
    #-------------------------------------------------------
    # definitions
    fdata = run_setup["regioninfo"]["regionname"]+".dat"
    fchno = run_setup["regioninfo"]["fchain"]
    day0  = run_setup["regioninfo"]["day0"]
    
    #-------------------------------------------------------
    # determine model type
    model_type = run_setup["mcmcopts"]["model_type"]
    
    if model_type == "oneWave":
        print("Running MCMC with one-wave model")
    elif model_type == "twoWave":
        print("Running MCMC with two-wave model")
    elif model_type == "threeWave":
        print("Running MCMC with three-wave model")
    else:
        sys.exit("Did not recognize 'model_type' specified in {}\n Options are 'oneWave', 'twoWave', or 'threeWave'".format(sys.argv[1]))
    
    #-------------------------------------------------------
    # determine error model type
    error_model_type = run_setup["mcmcopts"]["error_model_type"]
    
    if error_model_type == "add":
        print("Additive error model selected")
    elif error_model_type == "addMult":
        print("Additive & Multiplicative error model selected")
    else:
        sys.exit("Did not recognize 'error_model_type' specified in {}\n Options are 'add' or 'addMult'".format(sys.argv[1]))

    #-------------------------------------------------------
    # ensure that inputs are consistent with specified model
    # and error model type 
    
    def check_param_lengths(n_param,model_type,error_model_type):
            e_message=" needs to be length {} for a '{}' model with '{}' error model".format(n_param,model_type,error_model_type)
            assert len(run_setup["mcmcopts"]["cini"])==n_param,"cini"+e_message    
            assert len(run_setup["mcmcopts"]["cvini"])==n_param,"cvini"+e_message    
            assert len(run_setup["mcmcopts"]["spllo"])==n_param,"spllo"+e_message    
            assert len(run_setup["mcmcopts"]["splhi"])==n_param,"splhi"+e_message    
            assert len(run_setup["bayesmod"]["prior_types"])==n_param,"prior_types"+e_message    
            assert len(run_setup["bayesmod"]["prior_info"])==n_param,"prior_info"+e_message    
    
    
    if model_type=="oneWave":
        if error_model_type=="add":
            check_param_lengths(5,model_type,error_model_type)
        else:
            check_param_lengths(6,model_type,error_model_type)
    elif model_type=="twoWave":
        if error_model_type=="add":
            check_param_lengths(9,model_type,error_model_type)
        else:
            check_param_lengths(10,model_type,error_model_type)
    else:
        if error_model_type=="add":
            check_param_lengths(13,model_type,error_model_type)
        else:
            check_param_lengths(14,model_type,error_model_type)

    #-------------------------------------------------------
    # extract data from raw data
    rawdata = np.loadtxt(fdata,dtype=str)
    ndays = rawdata.shape[0]
    days_since_day0 = np.array([(parser.parse(rawdata[i,0])-parser.parse(day0)).days for i in range(ndays)])
    
    if "running_avg_obs" in run_setup["regioninfo"]:
        new_cases = runningAvg(np.array([float(rawdata[i,1]) for i in range(rawdata.shape[0])]),
                           run_setup["regioninfo"]["running_avg_obs"])
        print("Taking {}-day running average of observations".format(run_setup["regioninfo"]["running_avg_obs"]))
    else:
        new_cases = np.array([float(rawdata[i,1]) for i in range(rawdata.shape[0])])
    
    # get sigma for the incubation model
    incubation_median = run_setup["incopts"]["incubation_median"]
    incubation_sigma  = run_setup["incopts"]["incubation_sigma"]
    
    #-------------------------------------------------------
    # mcmc
    opts = {"nsteps": run_setup["mcmcopts"]["nsteps"],
            "nfinal": run_setup["mcmcopts"]["nfinal"],
            "gamma":  run_setup["mcmcopts"]["gamma"],
            "inicov": np.array(run_setup["mcmcopts"]["cvini"]),
            "spllo":  np.array(run_setup["mcmcopts"]["spllo"]),
            "splhi":  np.array(run_setup["mcmcopts"]["splhi"]),
            "logfile":run_setup["mcmcopts"]["logfile"],
            "nburn":1000,"nadapt":100,"coveps":1.e-10,"burnsc":5,
            "ndr":2,"drscale":[5,4,3],"ofreq":5000,"tmpchn":"tmpchn"
            }
    # This will allow full covariance matrices
    # to be passed for inicov, or just the diagonal [tporton]
    if len(opts['inicov'].shape)==1:
        opts['inicov'] = np.diag(opts['inicov'])
    
    if "incubation_model" in run_setup["incopts"]:
        inc_model = run_setup["incopts"]["incubation_model"]
    else:
        inc_model = "lognormal"
    
    error_weight = None
    if "error_weight" in run_setup["mcmcopts"]:
        print("Applying weighting to error term")
        error_weight = compute_error_weight(run_setup["mcmcopts"]["error_weight"],days_since_day0) 

    modelinfo={"model_type":       model_type,
               "error_model_type": error_model_type,
               "error_weight":     error_weight,
               "days_since_day0":  days_since_day0,
               "new_cases":        new_cases,
               "incubation_model": inc_model,
               "incubation_median":incubation_median,
               "incubation_sigma": incubation_sigma, 
               "inftype":          "gamma",
               "days_extra":       0,
               "prior_types":      run_setup["bayesmod"]["prior_types"],
               "prior_info":       run_setup["bayesmod"]["prior_info"]}
    
    # Convolution vs Quadrature:
    #   -The user can choose to use a fft convolution instead of 
    #    quadrature to perform the integration of Y(t)
    #   -default is set to zero if the user defines nothing
    #   -To set, add "useconv":1 to the mcmcopts in the *json file 
    if "useconv" in run_setup["mcmcopts"]:
        modelinfo["useconv"] = run_setup["mcmcopts"]["useconv"]
        if modelinfo["useconv"] == 1:
            print("Using fft convolution instead of quadrature")
    
    if "incubation_type" in run_setup["mcmcopts"]:
        modelinfo["incubation_type"] = run_setup["mcmcopts"]["incubation_type"]
        print("Using incubation type:",modelinfo["incubation_type"])
    else:
        print("Using fixed incubation type")
    
    # choose log-posterior function
    if "likl_type" in run_setup["mcmcopts"]:
        lliktype = run_setup["mcmcopts"]["likl_type"]
        print("Using %s likelihood"%(lliktype))
        if lliktype=="gaussian":
            lpf = logpost
        elif lliktype=="negative_binomial":
            lpf = logpost_negb
        elif lliktype=="poisson":
            lpf = logpost_poisson
        else:
            print("Unknown likelihood construction")
            quit()
    else:
        lpf = logpost
    
    if "likl_type" in run_setup["mcmcopts"]:
        lliktype = run_setup["mcmcopts"]["likl_type"]
        if lliktype=="poisson":
            modelinfo["sumLogK"] = sum([sum([np.log(i) for i in range(1,int(k)+1)]) for k in new_cases if k>0])
        
    sol=ammcmc(opts,np.array(run_setup["mcmcopts"]["cini"]),lpf,modelinfo)
    
    #-------------------------------------------------------------------------------------------
    # save mcmc output
    import h5py
    f = h5py.File(fchno, 'w')
    dset = f.create_dataset("chain",     data=sol['chain'],     compression="gzip")
    mdIn = f.create_dataset("mode",      data=sol['cmap'],      compression="gzip")
    mdPs = f.create_dataset("modepost",  data=[sol['pmap']],    compression="gzip")
    minf = f.create_dataset("minfo",     data=sol['minfo'],     compression="gzip")
    accr = f.create_dataset("accr",      data=[sol['accr']],    compression="gzip")
    fcov = f.create_dataset("final_cov", data=sol['final_cov'], compression="gzip")
    f.close()

if __name__ == '__main__':
    main(sys.argv[1])

