import sys
import json
import h5py
import numpy as np

def json_parameter_lists_from_chain(json_dict,n_sigma=3):
    """
    Read in chain data and compute parameter stats and bounds 
    """

    # obtain number of parameters to estimate stats and limits for
    model_type = json_dict["mcmcopts"]["model_type"]

    nparams = 0
    if model_type == "oneWave":
        nparams = 4
    elif model_type == "twoWave":
        nparams = 8
    elif model_type == "threeWave":
        nparams = 12
    else:
        sys.exit("ERROR: did not recognize model_type from json file. Should be oneWave, twoWave, or threeWave")

    # get chain file name
    fchno = json_dict["regioninfo"]["fchain"]

    # retrieve MCMC chain
    file = h5py.File(fchno, 'r')
    chn  = np.array(file["chain"])
    file.close()
    
    # sample MCMC chain
    nstart   = json_dict["ppopts"]["nstart"]
    nsamples = json_dict["ppopts"]["nsamples"]
    nskip    = (chn.shape[0]-nstart)//nsamples
    chnSamps=chn[nstart::nskip,:nparams]

    # set parameter statistics initial guesses
    cini = chnSamps.mean(axis=0)
    cvini = chnSamps.var(axis=0)

    # set parameter bounds
    spllo = cini - n_sigma * np.sqrt(cvini)
    splhi = cini + n_sigma * np.sqrt(cvini)
   
    # convert to lists
    cini = list(cini)
    cvini = list(cvini)
    spllo = list(spllo)
    splhi = list(splhi)

    return cini, cvini, spllo, splhi

def clip_json_parameter_lists(spllo,splhi,cini,nWaves=2):
    """
    Clip negative valus for non-negative parameters

    """
    non_negative_state_ind = []
    for i in range(nWaves):
        non_negative_state_ind += [1+4*i,2+4*i,3+4*i]
    
    for ind in non_negative_state_ind:
        spllo[ind] = max(spllo[ind],0.0) 
        splhi[ind] = max(splhi[ind],0.0) 
        cini[ind] = max(cini[ind],0.0) 

    return spllo,splhi,cini 

class PrimeJsonCreator:
    """
    Handles the creation and modifications of a new json dict based on an existing json dict
    """
    def __init__(self, json_in):
        """
        Initialize primeJson object by copying existing json dict.

        :param json_in: dictionary read from json input file
        :type json_in: dict
        """
        self.json_dict = json_in.copy() 

    def set_region_name(self, regionname, filesuffix=""):
        """
        :param regionname: region name to be used in new json dict
        :type regionname: string
        """
        self.json_dict["regioninfo"]["regionname"] = regionname
        
        # set log filenames to include updated regionname
        self.json_dict["regioninfo"]["fchain"] = regionname+"_mcmc"+filesuffix+".h5"
        self.json_dict["mcmcopts"]["logfile"] = "logmcmc"+regionname+filesuffix+".txt" 
        self.json_dict["ppopts"]["fpredout"] = regionname+"_epidemic_curve"+filesuffix
        self.json_dict["ppopts"]["fout_newcases"] = regionname+"_epidemic_curve"+filesuffix
        self.json_dict["infopts"]["finfout"] = regionname+"_infection_curve"+filesuffix
        self.json_dict["infopts"]["fout_inf"] = regionname+"_infection_curve"+filesuffix
        self.json_dict["csvout"]["finfcurve"] = regionname+"_infection_curve"+filesuffix
        self.json_dict["csvout"]["fnewcases"] = regionname+"_epidemic_curve"+filesuffix

    def set_model_type(self, model_type):
        self.json_dict["mcmcopts"]["model_type"] = model_type
    
    def set_gamma(self, gamma):
        self.json_dict["mcmcopts"]["gamma"] = gamma

    def set_parameter_limits(self, spllo, splhi):
        self.json_dict["mcmcopts"]["spllo"] = spllo 
        self.json_dict["mcmcopts"]["splhi"] = splhi 

    def set_parameter_stats(self, cini, cvini):
        self.json_dict["mcmcopts"]["cini"] = cini 
        self.json_dict["mcmcopts"]["cvini"] = cvini 

    def set_prior_distributions(self, prior_types, prior_info):
        self.json_dict["bayesmod"]["prior_types"] = prior_types 
        self.json_dict["bayesmod"]["prior_info"] = prior_info

    def get_json_dict(self,filename):

        """
        Get modified json dict

        :return: json_out, a new dictionary to be written to a new json file
        :rtype: dict
    
        """
        return self.json_dict

    def write_json_dict(self, filename):
        """
        Write json dict to a new json file with a specified filename

        :param filename: file name to be used for new json file
        :type filename: string
        """
        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(self.json_dict, f, ensure_ascii=False, indent=4)

