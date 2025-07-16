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

def addModelDiscr(pred,chn,setup):
    model_type = setup["mcmcopts"]['model_type']
    error_model_type  = setup["mcmcopts"]['error_model_type']

    # if not specified return
    if "postpred" not in setup["ppopts"]:
        return pred

    # if pp not 1 return
    if setup["ppopts"]["postpred"] != 1:
        return pred

    # posterior-predictive
    if "likl_type" in setup["mcmcopts"] and setup["mcmcopts"]["likl_type"]=="negative_binomial":
        alpha_ind = 4
        if model_type == "twoWave":
            alpha_ind = 8
        elif model_type == "threeWave":
            alpha_ind = 12

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
        else:
            print("WARNING: Error model not recognized, not adding error to prediction")

    return pred


#-------------------------------------------------------
run_setup=json.load(open(sys.argv[1]))
print(run_setup)

#-------------------------------------------------------
# definitions
fdata = run_setup["regioninfo"]["regionname"]+".dat"
fchno = run_setup["regioninfo"]["fchain"]
day0  = run_setup["regioninfo"]["day0"]

model_type        = run_setup["mcmcopts"]['model_type']
error_model_type  = run_setup["mcmcopts"]['error_model_type']


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
    rawdata = np.loadtxt(run_setup["ppopts"]["newdata"],dtype=str)
    nc_new = [np.array([parser.parse(rawdata[i,0]) for i in range(rawdata.shape[0])])]
    nc_new.append(np.array([float(rawdata[i,1]) for i in range(rawdata.shape[0])]))


#-------------------------------------------------------
# extract data from raw data
rawdata = np.loadtxt(fdata,dtype=str)
days_since_day0 = np.array([(parser.parse(rawdata[i,0])-parser.parse(day0)).days for i in range(rawdata.shape[0])])
new_cases = np.array([float(rawdata[i,1]) for i in range(rawdata.shape[0])])

#-------------------------------------------------------
# get sigma for the incubation model
incubation_median = run_setup["incopts"]["incubation_median"]
incubation_sigma  = run_setup["incopts"]["incubation_sigma"]

if "incubation_model" in run_setup["incopts"]:
    inc_model = run_setup["incopts"]["incubation_model"]
else:
    inc_model = "lognormal"

#-------------------------------------------------------
# compute predictions -> push forward or posterior predictive pdf
fh5 = prediction_filename(run_setup)

if run_setup["ppopts"]["runmodel"]:
    modelinfo={"model_type":model_type,
               "error_model_type":error_model_type,
               "error_weight":None,
               "days_since_day0":days_since_day0,
               "new_cases":new_cases,
               "incubation_model": inc_model,
               "incubation_median":incubation_median,
               "incubation_sigma":incubation_sigma,
               "inftype":run_setup["infopts"]["inftype"],
               "days_extra":run_setup["ppopts"]["days_extra"],
               "day0":run_setup["regioninfo"]["day0"]}
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
    pred=np.array([modelPred(chnPred[k],modelinfo,is_cdf=True) for k in range(chnPred.shape[0])])
    pred1 = np.zeros(pred.shape)
    for j in range(1,pred.shape[1]):
        pred1[:,j]=pred[:,j]-pred[:,j-1]

    pred = pred1.copy()
    # posterior predictive
    pred = addModelDiscr(pred,chnPred,run_setup)

    # assemble set of dates 
    f = h5py.File(fh5, 'w')
    dset1 = f.create_dataset("predictions", data=pred,  compression="gzip")
    f.close()
else:
    # retrieve MCMC chain
    file   = h5py.File(fh5, 'r')
    pred   = np.array(file["predictions"])
    file.close()

datesData = np.array([parser.parse(rawdata[i,0]) for i in range(rawdata.shape[0])])
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
plot_post_pred(datesPred,pred,datesData,new_cases,qntList,normalize,iendData,run_setup,nc_new=nc_new)

#-------------------------------------------------------
# infection rate

fh5 = run_setup["infopts"]["finfout"]+".h5"
if run_setup["infopts"]["runmodel"]:
    modelinfo={"model_type":model_type,
               "error_model_type":error_model_type,
               "inftype":run_setup["infopts"]["inftype"],
               "day0":run_setup["regioninfo"]["day0"],
               "ndays":run_setup["infopts"]["ndays"]}
    nstart   = run_setup["ppopts"]["nstart"]
    nsamples = run_setup["ppopts"]["nsamples"]
    nskip    = (chn.shape[0]-nstart)//nsamples
    chnPred  = chn[nstart::nskip,:]
    # infection curves
    infect=np.array([infection(chnPred[k],modelinfo) for k in range(chnPred.shape[0])])
    # assemble set of dates 
    f = h5py.File(fh5, 'w')
    dset2 = f.create_dataset("infect", data=infect,  compression="gzip")
    f.close()
else:
    # retrieve MCMC chain
    file   = h5py.File(fh5, 'r')
    infect = np.array(file["infect"])
    file.close()

ndaysinf=(parser.parse(rawdata[-1,0])-parser.parse(rawdata[0,0])).days+run_setup["ppopts"]["days_extra"]+1
datesmean=np.array([parser.parse(rawdata[0,0])+datetime.timedelta(days=i) for i in range(ndaysinf)])

infall = np.zeros((infect.shape[0],ndaysinf))
for i in range(infect.shape[0]):
	for j in range(infect.shape[2]):
		dtobj = datetime.datetime.fromtimestamp(infect[i,0,j])
		if len(np.where(datesmean==dtobj)[0])>0:
			posID = np.where(datesmean==dtobj)[0][0]
			infall[i,posID] = infect[i,1,j]


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
if error_model_type=="add":
    fnewc=fnewc+"_a"
    finfc=finfc+"_a"
else:
    fnewc=fnewc+"_am"
    finfc=finfc+"_am"

if run_setup["ppopts"]["postpred"] == 1:
    fnewc=fnewc+"_pp"
else:
    fnewc=fnewc+"_pf"

fnewc=fnewc+".csv"
finfc=finfc+".csv"

output_epicurves(pred,datesPred,new_cases,nskip,qlist,fnewc)
output_infcurves(infall,datesmean,nskip,qlist,finfc)

