#!/usr/bin/env /Users/csafta/miniconda3/bin/python

import sys
import numpy as np
import h5py
import json
from scipy.stats import gamma,lognorm,norm
from numpy.random import negative_binomial as nb
import matplotlib.pyplot as plt
from   dateutil import parser
import datetime
from prime_utils import prediction_filename

# sys.path.append('../../scripts')
from prime_model      import modelPred
from prime_infection  import infection
from prime_plot       import plot_post_pred, plot_infection_curve

cLen = 10

def addModelDiscr(pred,chn, setup, ireg="corr"):
    num_waves = setup["modelopts"]['num_waves']
    error_model_type  = setup["bayesmod"]['error_model_type']

    # if not specified return
    if "postpred" not in setup["ppopts"]:
        return pred

    # if pp not 1 return
    if setup["ppopts"]["postpred"] != 1:
        return pred

    # posterior-predictive
    if "lpf_type" in setup["bayesmod"] and setup["bayesmod"]["lpf_type"]=="negative_binomial":
        alpha_ind = 4 * num_waves

        for k in range(chn.shape[0]):
            alpha = np.exp(chn[k,alpha_ind])
            p     = alpha/(alpha+pred[k,1:]+1.e-4)
            if np.any(p<0):
                print(k,alpha,p,pred[k])
                quit()
            pred[k,1:] = nb(alpha,p)
    else:
        if error_model_type=="add":
            pred = np.array([pred[k]+np.exp(chn[k,-1])*np.random.normal() for k in range(chn.shape[0])])
        elif error_model_type=="addMult":
            pred = np.array([pred[k]+(np.exp(chn[k,-2])+np.exp(chn[k,-1])*pred[k])*np.random.normal() for k in range(chn.shape[0])])
            # if ireg == "corr":
            #     pred = np.array([pred[k]+(np.exp(chn[k,-4])+np.exp(chn[k,-3])*pred[k])*np.random.normal() for k in range(chn.shape[0])])
            # else:
            #     pred = np.array([pred[k]+(np.exp(chn[k,-2])+np.exp(chn[k,-1])*pred[k])*np.random.normal() for k in range(chn.shape[0])])
        else:
            print("WARNING: Error model not recognized, not adding error to prediction")

    return pred

#-------------------------------------------------------
#run_setup=json.load(open(sys.argv[1]))
#print(run_setup)
from prime_run import get_opts

assert len(sys.argv)>2

setupfile=sys.argv[1]
ireg = int(sys.argv[2])

run_opts, run_setup = get_opts(setupfile, return_run_setup=True)

#-------------------------------------------------------
# definitions
#fdata = run_setup["regioninfo"]["count_data"]
fchno = "".join(run_opts["region_tag"])+"_mcmc.h5"
day0  = run_opts["day0"]

error_model_type  = run_opts["error_model_type"]


#-------------------------------------------------------
# For emcee testing. uncomment if using run_emcee
# fhno = 'emcee_' + run_setup["regioninfo"]["fchain"]

#-------------------------------------------------------
# retrieve MCMC chain
file = h5py.File(fchno, 'r')
chn  = np.array(file["chain"])
file.close()

#-------------------------------------------------------
nc_new=None
if "newdata" in run_setup["ppopts"]:
    rawdata = np.loadtxt(run_setup["ppopts"]["newdata"],delimiter=",",dtype=str)
    nc_new = [np.array([parser.parse(rawdata[i,0]) for i in range(rawdata.shape[0])])]
    nc_new.append(np.array([float(rawdata[i,1]) for i in range(rawdata.shape[0])]))
    day0_new = np.array([(parser.parse(rawdata[i,0])-parser.parse(day0)).days for i in range(rawdata.shape[0])])
    counts_new = np.array([float(rawdata[i,1]) for i in range(rawdata.shape[0])])
    if "running_avg_obs" in run_opts:
        from prime_utils import runningAvg
        nc_new[-1] = runningAvg(nc_new[-1], run_opts["running_avg_obs"])


#-------------------------------------------------------
# extract data from raw data
from prime_run import get_counts
days_since_day0, daily_counts, rawdata = get_counts(run_opts,return_raw_data=True)
# rawdata = np.loadtxt(fdata,delimiter=",",dtype=str)
# days_since_day0 = np.array([(parser.parse(rawdata[i,0])-parser.parse(day0)).days for i in range(rawdata.shape[0])])
# daily_counts = np.array([float(rawdata[i,1]) for i in range(rawdata.shape[0])])
new_cases = np.array([float(rawdata[ireg][i,1]) for i in range(rawdata[ireg].shape[0])])
num_regions = len(daily_counts)

#-------------------------------------------------------
# compute predictions -> push forward or posterior predictive pdf
fh5 = prediction_filename(run_setup)

if run_setup["ppopts"]["runmodel"]:
    modelinfo={"num_waves":        run_setup["modelopts"]["num_waves"],
               "error_model_type": error_model_type,
               "days_since_day0":  days_since_day0,
               "daily_counts":     daily_counts,
               # "days_since_day0":  day0_new,
               # "daily_counts":     counts_new,
               "incubation_model": run_opts["inc_model"],
               "incubation_median":run_opts["inc_median"],
               "incubation_sigma": run_opts["inc_sigma"], 
               "incubation_type":  run_opts["inc_type"], 
               "inftype":          "gamma",
               "useconv":          run_opts["useconv"],
               "days_extra":       run_opts["days_extra"],
               "day0":             run_opts["day0"]}
    if "incubation_type" in run_setup["mcmcopts"]:
        modelinfo["incubation_type"] = run_setup["mcmcopts"]["incubation_type"]
        print("Incubation type:",modelinfo["incubation_type"])
    
    # if "useconv" in run_setup["mcmcopts"]:
    #     modelinfo["useconv"] = run_setup["mcmcopts"]["useconv"]
    #     if modelinfo["useconv"] == 1:
    #         print("Using fft convolution instead of quadrature")


    nstart   = run_setup["ppopts"]["nstart"]
    nsamples = run_setup["ppopts"]["nsamples"]
    nskip    = (chn.shape[0]-nstart)//nsamples
    pred=[]
    chnPred=chn[nstart::nskip,:]
    print("=======================================")
    print("Mean N/2sigma:",chnPred[:,1].mean(),np.std(chnPred[:,1]))
    print("=======================================")

    # pushed forward pdf
    # XXX
    pred=np.array([modelPred(chnPred[k],modelinfo,n_regions=num_regions,is_cdf=True)[ireg] for k in range(chnPred.shape[0])])
    print('Pushed forward predictions shape: ',pred.shape)
    pred1 = np.zeros(pred.shape)
    for j in range(1,pred.shape[1]):
        pred1[:,j]=pred[:,j]-pred[:,j-1]

    pred = pred1.copy()
    #print(counts_new.shape,day0_new.shape,pred.shape)
    if run_setup["ppopts"]["postpred"] == 1:
        # posterior predictive
        if ireg<3:
            pred = addModelDiscr(pred,chnPred,run_setup,ireg="corr")
        else:
            pred = addModelDiscr(pred,chnPred,run_setup,ireg="uncorr")


    # assemble set of dates 
    f = h5py.File(fh5, 'w')
    dset1 = f.create_dataset("predictions", data=pred,  compression="gzip")
    f.close()
else:
    # retrieve MCMC chain
    file   = h5py.File(fh5, 'r')
    pred   = np.array(file["predictions"])
    file.close()

datesData = np.array([parser.parse(rawdata[0][i,0]) for i in range(rawdata[0].shape[0])])
datesPred = np.concatenate((datesData,
            np.array([datesData[-1]+datetime.timedelta(days=i+1) 
            for i in range(run_setup["ppopts"]["days_extra"])])))

# colormap settings
import matplotlib as mpl
cmap1 = mpl.cm.PuBu
cmap2 = mpl.cm.PuRd

qntList = [0.025]+[0.05*i for i in range(1,20)]+[0.975]
normalize = mpl.colors.Normalize(vmin=0.025,vmax=0.5)
iendData = np.where(datesPred==datesData[-1])[0][0]+1

#--------------------------------------------------------------------------------
# Plot push-forward/posterior predictive PDF
plot_post_pred(datesPred,pred*run_opts["population_data"][ireg]*1e6,datesData,daily_counts[ireg]*run_opts["population_data"][ireg]*1e6,qntList,normalize,iendData,run_setup,nc_new=nc_new)
#quit()

#-------------------------------------------------------
# infection rate

fh5 = run_setup["infopts"]["finfout"]+".h5"
if run_setup["infopts"]["runmodel"]:
    modelinfo={"num_waves": run_setup["modelopts"]["num_waves"],
               "error_model_type":error_model_type,
               "inftype":run_setup["infopts"]["inftype"],
               "day0":run_setup["regioninfo"]["day0"],
               "ndays":run_setup["infopts"]["ndays"]}
    nstart   = run_setup["ppopts"]["nstart"]
    nsamples = run_setup["ppopts"]["nsamples"]
    nskip    = (chn.shape[0]-nstart)//nsamples
    chnPred  = chn[nstart::nskip,:]

    # infection curves
    #print(run_opts["population_data"][ireg])
    infect = np.array([infection(chnPred[k],modelinfo,regionID=ireg) for k in range(chnPred.shape[0])])
    #print(infect.max(),infect.min())
    
    # assemble set of dates 
    f = h5py.File(fh5, 'w')
    dset2 = f.create_dataset("infect", data=infect,  compression="gzip")
    f.close()
else:
    # retrieve MCMC chain
    file   = h5py.File(fh5, 'r')
    infect = np.array(file["infect"])
    file.close()

print("First/End Dates: ",[parser.parse(rawdata[ireg][i,0]) for i in [0,-1]])
#print(parser.parse(rawdata[ireg][0,0]))
ndaysinf=(parser.parse(rawdata[ireg][-1,0])-parser.parse(rawdata[ireg][0,0])).days+run_setup["ppopts"]["days_extra"]+1
datesmean=np.array([parser.parse(rawdata[ireg][0,0])+datetime.timedelta(days=i) for i in range(ndaysinf)])

infall = np.zeros((infect.shape[0],ndaysinf))
for i in range(infect.shape[0]):
	for j in range(infect.shape[2]):
		dtobj = datetime.datetime.fromtimestamp(infect[i,0,j])
		if len(np.where(datesmean==dtobj)[0])>0:
			posID = np.where(datesmean==dtobj)[0][0]
			infall[i,posID] = infect[i,1,j]*run_opts["population_data"][ireg]

iendData = np.where(datesmean==datesData[-1])[0][0]+1
plot_infection_curve(datesmean,infall,qntList,normalize,iendData,run_setup)

#------------------------------------------------------------
# save csv files
from prime_utils import output_epicurves, output_infcurves

nskip = run_setup["csvout"]["nskip"]
#qlist = run_setup["csvout"]["qlist"]
qlist = qntList
fnewc = run_setup["csvout"]["fnewcases"]
finfc = run_setup["csvout"]["finfcurve"]
# if error_model_type=="add":
#     fnewc=fnewc+"_a"
#     finfc=finfc+"_a"
# else:
#     fnewc=fnewc+"_am"
#     finfc=finfc+"_am"

if run_setup["ppopts"]["postpred"] == 1:
    fnewc=fnewc +"_pp"
else:
    fnewc=fnewc+"_pf"

fnewc=fnewc+".csv"
finfc=finfc+".csv"

output_epicurves(pred*run_opts["population_data"][ireg]*1e6,datesPred,new_cases,nskip,qlist,fnewc)
output_infcurves(infall,datesmean,nskip,qlist,finfc)

