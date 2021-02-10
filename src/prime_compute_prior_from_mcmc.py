import sys
import os
import numpy as np
import h5py
import json
from scipy.stats     import gamma,lognorm,norm
import matplotlib.pyplot as plt
from dateutil import parser
import datetime
from prime_json_utils import *

sys.path.append('../../scripts')

'''
This script postprocesses data from a PRIME model
output to determine the priors for another prime model
with the same or more waves. It then writes a
new json file with the prior data entered.

Inputs:
-json run setup script
-model_type for desired prior/json script (oneWave, twoWave, or threeWave)
-hdf5 file with MCMC chain data (name from json file)


Outputs:
-New json file with the same name as the input with
the suffix "_<model_type_out> appended to the filename. 

Notes: 
-Needs to be run in directory where a model has been run. 
-Priors set using mean and variance computed from chain
-Upper and lower parameter values set to 3 times the 
standard deviation computed from the chain data
-Error model parameters are not modified by this script
-Default values are in place for the relative 
start day of the 2nd (or 3rd) wave. It is recommended to modify 
these on a case-by-case basis

'''

#-------------------------------------------------------
run_setup_in=json.load(open(sys.argv[1]))
print(run_setup_in)

model_type_out=sys.argv[2]

if model_type_out == "oneWave":
    print("Generating json file with MCMC-informed prior for one-wave model")
    nWaves = 1
elif model_type_out == "twoWave":
    print("Generating json file with MCMC-informed prior for two-wave model")
    nWaves = 2
elif model_type_out == "threeWave":
    print("Generating json file with MCMC-informed prior for three-wave model")
    nWaves = 3
else:
    sys.exit("Did not recognize 'model_type' in input arguments, was: {}\n Options are 'oneWave', 'twoWave', or 'threeWave'".format(sys.argv[2]))


#-------------------------------------------------------
# definitions
region_name = run_setup_in["regioninfo"]["regionname"]
model_type_in = run_setup_in["mcmcopts"]["model_type"]
error_model_type = run_setup_in["mcmcopts"]["error_model_type"]

# estimate first wave parameter statistics and limits from MCMC chain
cini, cvini, spllo, splhi = json_parameter_lists_from_chain(run_setup_in)

# prior distributions for each model type
if model_type_out == "oneWave":
    if model_type_in == "oneWave":
        # Generate oneWave json file from oneWave MCMC chain

        # create prior distribution lists
        prior_types = ["g","u","u","u"]
        prior_info =  [[cini[0],np.sqrt(cvini[0])],[0,1],[0,1],[0,1]]
    else: 
        sys.exit("ERROR: can only generate MCMC-informed prior for a oneWave model using data from a oneWave model")
elif model_type_out == "twoWave":
    if model_type_in == "oneWave":
        # Generate twoWave json file from oneWave MCMC chain

        # append wave and error model prior from second wave
        cini += run_setup_in["mcmcopts"]["cini"]
        cvini += run_setup_in["mcmcopts"]["cvini"]
        spllo += run_setup_in["mcmcopts"]["spllo"]
        splhi += run_setup_in["mcmcopts"]["splhi"]

        # hardcode dt2 stats and limits
        cini[4] = 60
        cvini[4] = 225
        spllo[4] = 20
        splhi[4] = 100
        print("Warning: dt2 stats and limits set automatically. Recommend tuning these parameters (5th entry of cini, cvini, spllo, splhi, prior_types, prior_info) as needed. ")
        
        # create prior distribution lists
        prior_types = ["g","u","u","u","g","u","u","u"]
        prior_info =  [[cini[0],np.sqrt(cvini[0])],[0,1],[0,1],[0,1],[cini[4],np.sqrt(cvini[4])],[0,1],[0,1],[0,1]]
    elif model_type_in == "twoWave":
        # Generate twoWave json file from twoWave MCMC chain

        # append error model prior from second wave
        cini += run_setup_in["mcmcopts"]["cini"][8:]
        cvini += run_setup_in["mcmcopts"]["cvini"][8:]
        spllo += run_setup_in["mcmcopts"]["spllo"][8:]
        splhi += run_setup_in["mcmcopts"]["splhi"][8:]

        # create prior distribution lists
        prior_types = ["g","u","u","u","g","u","u","u"]
        prior_info =  [[cini[0],np.sqrt(cvini[0])],[0,1],[0,1],[0,1],[cini[4],np.sqrt(cvini[4])],[0,1],[0,1],[0,1]]
    else: 
        sys.exit("ERROR: can only generate MCMC-informed prior for a twoWave model using data from a oneWave or twoWave model")
elif model_type_out == "threeWave":
    if model_type_in == "oneWave":
        sys.exit("ERROR: can only generate MCMC-informed prior for a threeWave model using data from a twoWave or threeWave model")
    elif model_type_in == "twoWave":
        # Generate threeWave json filr from twoWave MCMC chain

        # add priors for third wave
        cini += run_setup_in["mcmcopts"]["cini"][4:8]
        cvini += run_setup_in["mcmcopts"]["cvini"][4:8]
        spllo += run_setup_in["mcmcopts"]["spllo"][4:8]
        splhi += run_setup_in["mcmcopts"]["splhi"][4:8]
        
        # append error model prior from second wave
        cini += run_setup_in["mcmcopts"]["cini"][8:]
        cvini += run_setup_in["mcmcopts"]["cvini"][8:]
        spllo += run_setup_in["mcmcopts"]["spllo"][8:]
        splhi += run_setup_in["mcmcopts"]["splhi"][8:]
        
        # hardcode dt3 stats and limits
        cini[8] = 180
        cvini[8] = 225
        spllo[8] = 140
        splhi[8] = 220
        print("Warning: dt3 stats and limits set automatically. Recommend tuning these parameters (9th entry of cini, cvini, spllo, splhi, prior_types, prior_info) as needed. ")

        # create prior distribution lists
        prior_types = ["g","u","u","u","g","u","u","u","g","u","u","u"]
        prior_info =  [[cini[0],np.sqrt(cvini[0])],[0,1],[0,1],[0,1],\
                       [cini[4],np.sqrt(cvini[4])],[0,1],[0,1],[0,1],\
                       [cini[8],np.sqrt(cvini[8])],[0,1],[0,1],[0,1]]
    else:
        # Generate threeWave json file from three wave MCMC chain

        # append error model prior
        cini += run_setup_in["mcmcopts"]["cini"][12:]
        cvini += run_setup_in["mcmcopts"]["cvini"][12:]
        spllo += run_setup_in["mcmcopts"]["spllo"][12:]
        splhi += run_setup_in["mcmcopts"]["splhi"][12:]

        # create prior distribution lists
        prior_types = ["g","u","u","u","g","u","u","u","g","u","u","u"]
        prior_info =  [[cini[0],np.sqrt(cvini[0])],[0,1],[0,1],[0,1],\
                       [cini[4],np.sqrt(cvini[4])],[0,1],[0,1],[0,1],\
                       [cini[8],np.sqrt(cvini[8])],[0,1],[0,1],[0,1]]


# clip negative values for non-negative parameters
spllo,splhi,cini = clip_json_parameter_lists(spllo,splhi,cini,nWaves=nWaves)

# append appropriate error model priors
if error_model_type=="add":
    prior_types += ["u"]
    prior_info += [[0,1]]
elif error_model_type=="addMult":
    prior_types += ["u","u"]
    prior_info += [[0,1],[0,1]]
else:
    sys.exit("error_model_type should be 'add' or 'addMult'")

# create json file for second wave, leaving entry for start date empty
json_out = PrimeJsonCreator(run_setup_in)

# modify entries for 2-wave model json file
json_out.set_region_name(region_name)
json_out.set_model_type(model_type_out)
json_out.set_gamma(0.7)
json_out.set_parameter_limits(spllo,splhi)
json_out.set_parameter_stats(cini,cvini)
json_out.set_prior_distributions(prior_types,prior_info)

# write new json file for second wave model
json_filename = "setup_" + region_name + "_" + model_type_out  + ".json"
if os.path.exists(json_filename):
    json_filename += ".new"
json_out.write_json_dict(json_filename)
