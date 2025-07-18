import sys,json
import numpy as np
import h5py
import pickle as pkl
from prime_stats import getKDE
from prime_plot  import plKDE


def main(filename):
    r"""
    Plots 1D and 2D marginal kernel density estimates based on MCMC samples

    Parameters
    ----------
    filename: string
        json file (.json) including run setup information and 
        postprocessing information for an MCMC run. It should 
        specify the name of the file containing the MCMC chain

        or
        
        pickle file (.pkl) with a dictionary containing the 
        KDE distributions.This file is generated by running 
        this script with a json file (see above)    
    """
    if filename[-3:] == "pkl":

        print("Reading KDE data from {}".format(filename))
    
        f = open(filename,"rb")
        kdedict = pkl.load(f)
        f.close()
    
        model_type       = kdedict['model_type']      
        error_model_type = kdedict['error_model_type']
    
    else:

        run_setup=json.load(open(filename))
        print(run_setup)
        
        #-------------------------------------------------------
        # definitions
        fchno = "".join(run_setup["regioninfo"]["region_tag"])+"_mcmc.h5"
        
        print("Computing KDE data from {}".format(fchno))
    
        # retrieve MCMC chain
        file = h5py.File(fchno, 'r')
        chn  = np.array(file["chain"])
        file.close()
    
        # compute KDE
        #chn[:,-4] = np.exp(chn[:,-4])
        #chn[:,-2] = np.exp(chn[:,-2])
        #chn[:,-1] = np.exp(chn[:,-1])
        nstart  = run_setup["ppopts"]["nstart"]
        kdedict = getKDE(chn,nskip=nstart,nthin=50,npts=100,bwfac=1.0)
       
        # Add data on model type to KDE
        error_model_type = run_setup["bayesmod"]["error_model_type"]
        kdedict["model_type"] = run_setup["modelopts"]["num_waves"]
        kdedict["error_model_type"] = error_model_type
        kdedict["regions"] =run_setup["regioninfo"]["region_tag"]
    
        f = open("kdedict.pkl","wb")
        pkl.dump(kdedict,f)
        f.close()

    # Model variable labels
    vnames =["t_0",\
                "N",\
                "k",\
                "\\theta"]*len(kdedict["regions"])
    
    # Error model labels
    if error_model_type == "add":
        vnames += ["\log\,\\tau_{\phi}^2","\lambda","\log\,\sigma_a"]
    elif error_model_type == "addMult":
        vnames += ["\log\,\\tau_{\phi}^2","\lambda","\log\,\sigma_a","\log\,\sigma_m"]
    else:
        sys.exit("Invalid error_model_type, select 'add' or 'addMult'")
    
    # Plot KDE
    plKDE(kdedict,'kde.pdf',vnames=vnames,ncont=31,ds=0.03)


if __name__ == '__main__':
    main(sys.argv[1])

