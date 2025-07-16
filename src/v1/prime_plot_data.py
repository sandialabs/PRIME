import sys
import numpy as np
import h5py
import json
import matplotlib.pyplot as plt
import datetime
from dateutil    import parser
from prime_utils import runningAvg

VERBOSE=0

def main(setupfile):
    r"""
    Plot raw and filtered data for the region specified in the setupfile.


    Parameters
    ----------
    setupfile: string
        json file (.json) including the region name. The "regionname.dat"
        should exist in the path accessible for this script
    """
    #-------------------------------------------------------
    run_setup=json.load(open(setupfile))
    if VERBOSE>0:
    	print(run_setup)
    
    #-------------------------------------------------------
    # definitions
    fdata = run_setup["regioninfo"]["regionname"]+".dat"
    day0  = run_setup["regioninfo"]["day0"]
    ndays_average = run_setup["regioninfo"]["running_avg_obs"]
    
    #-------------------------------------------------------
    # get raw data
    rawdata = np.loadtxt(fdata,dtype=str)
    ndays_data = rawdata.shape[0]
    days_since_day0 = np.array([(parser.parse(rawdata[i,0])-parser.parse(day0)).days for i in range(rawdata.shape[0])])
    new_cases = np.array([float(rawdata[i,1]) for i in range(ndays_data)])
    datesData = np.array([parser.parse(rawdata[i,0]) for i in range(ndays_data)])
    
    #-------------------------------------------------------
    # Filtered Data (with a running average)
    new_cases_filter = runningAvg(np.array([float(rawdata[i,1]) for i in range(ndays_data)]),ndays_average)
    
    #-------------------------------------------------------
    # plot
    
    fig = plt.figure(figsize=(6,6))
    ax=fig.add_axes([0.18,0.26,0.75,0.7])
    
    maxVal = np.max(new_cases)
    
    symtype=run_setup["ppopts"]["newcases"][0]
    ms=run_setup["ppopts"]["newcases"][1]
    
    pl1=ax.plot(datesData,new_cases,symtype,ms=ms)
    
    
    plt.plot(datesData,new_cases,'.k',label='Data')
    
    # TODO should change label depending on input deck options
    filter_label = str(ndays_average)+"-day running average"
    plt.plot(datesData,new_cases_filter,'s-r',label=filter_label)
    
    plt.xticks(rotation=45)
    
    x0=parser.parse(run_setup["ppopts"]["xylim_newcases"][0])
    x1=parser.parse(run_setup["ppopts"]["xylim_newcases"][1])
    
    x1=parser.parse(rawdata[-1,0])+datetime.timedelta(days=run_setup["ppopts"]["days_extra"]) 
    ax.set_xlim([x0,x1])
    
    y0 = 0.0
    y1 = 1.1*maxVal
    
    ax.set_ylim([y0,y1])
    
    xlbs=run_setup["ppopts"]["xylbl_newcases"][0]
    xlbf=run_setup["ppopts"]["xylbl_newcases"][1]
    ax.set_xlabel(xlbs,fontsize=xlbf)
    ylbs=run_setup["ppopts"]["xylbl_newcases"][2]
    ylbf=run_setup["ppopts"]["xylbl_newcases"][3]
    ax.set_ylabel(ylbs,fontsize=ylbf)
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(run_setup["ppopts"]["xyticklbl_newcases"][0])
    
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(run_setup["ppopts"]["xyticklbl_newcases"][1])
    
    plt.legend(fontsize=run_setup["ppopts"]["xyticklbl_newcases"][0])
    
    figname=run_setup["regioninfo"]["regionname"] + "_case_data" 
    plt.savefig(figname+"."+run_setup["ppopts"]["figtype"])
    plt.close() 
    
if __name__ == '__main__':
    main(sys.argv[1])
