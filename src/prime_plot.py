import numpy as np
import matplotlib        as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import os

from   dateutil import parser


def plKDE(kdedic,figname,vnames=None,fsizex=12,fsizey=12,ncont=41,fs1=18,fs2=12,
          xs=0.05,ys=0.05,ds=0.01,xe=0.04,ye=0.04):

    mpl.rcParams['xtick.labelsize'] = fs2
    mpl.rcParams['ytick.labelsize'] = fs2

    nvars = len(kdedic['x1D'])
    sx = (1-(nvars-1)*ds-xs-xe)/nvars
    sy = (1-(nvars-1)*ds-ys-ye)/nvars
    majorFormatter = FormatStrFormatter('%6.0e')

    # read distance correlation if availble
    distCorr = False
    if os.path.isfile("distanceCorr.txt"):
        distCorrMat = np.loadtxt("distanceCorr.txt")
        distCorr = True

    # create figure
    fig = plt.figure(figsize=(fsizex,fsizey))
    # place diagonal plots
    subs=[]
    # add diagonal plots
    for i in range(nvars):
        subs.append(fig.add_axes([xs+i*(sx+ds),ys+(nvars-1-i)*(sy+ds),sx,sy]))
    # add lower triangular plots
    for i in range(nvars-1):
        for j in range(i+1):
            subs.append(fig.add_axes([xs+j*(sx+ds),ys+(nvars-2-i)*(sy+ds),sx,sy]))
            if distCorr:
                subs[-1].set_title("{:.2}".format(distCorrMat[i+1,j]),fontweight='bold')
    subsnp=np.array(subs)
    # Plot 1D pdfs
    for i in range(nvars):
        pl1d=subsnp[i].plot(kdedic['x1D'][i],kdedic['p1D'][i])
    # Plot 2D pdfs
    for i in range(nvars*(nvars-1)//2):
        pl2d=subsnp[nvars+i].contour(kdedic['x2D'][i],kdedic['y2D'][i],kdedic['p2D'][i],ncont)
    # no ticks and ticklabels
    for i in range(nvars,len(subsnp)):
        xti=subsnp[i].set_xticks([])
        yti=subsnp[i].set_yticks([])
    # Set y labels on the right for diagonal plots
    if vnames is not None:
        for i in range(nvars):
            subsnp[i].xaxis.tick_top()
            subsnp[i].yaxis.tick_right()
            subsnp[i].yaxis.set_label_position("right")
            subsnp[i].set_ylabel(r'$p('+vnames[i]+')$', fontsize=fs1,rotation=90)
        icnt=0
        for i in range(len(subsnp)-nvars+1,len(subsnp)):
            icnt += 1
            subsnp[i].set_xlabel(r'$'+vnames[icnt-1]+'$', fontsize=fs1)
        icnt=nvars
        for i in range(nvars,2*nvars-1):
            icnt += i-nvars
            subsnp[icnt].set_ylabel(r'$'+vnames[i-nvars+1]+'$',fontsize=fs1, rotation=90)
        subsnp[nvars-1].set_xlabel(r'$'+vnames[nvars-1]+'$', fontsize=fs1)
    # no ticks or ticklabels on the diagonal plots
    for i in range(nvars):
        xti=subsnp[i].set_xticks([])
        yti=subsnp[i].set_yticks([])
    # save
    plt.savefig(figname)
    plt.show()


#--------------------------------------------------------------------------------
# Plot push-forward/posterior predictive PDF
def plot_post_pred(datesPred,pred,datesData,new_cases,
    qntList,normalize,iendData,run_setup,nc_new=None):

    # colormap settings
    import matplotlib as mpl
    cmap1 = mpl.cm.PuBu
    cmap2 = mpl.cm.PuRd

    fig = plt.figure(figsize=(6,6))
    ax=fig.add_axes([0.18,0.26,0.75,0.7])

    maxVal = -1.e100
    for i in range(len(qntList)-1):
        qnt0 = np.quantile(pred,qntList[i],axis=0)
        qnt1 = np.quantile(pred,qntList[i+1],axis=0)
        midPt = 0.5*(qntList[i]+qntList[i+1])
        alph = 0.8
        if qntList[i] >= 0.5:
            pl2 = ax.fill_between(datesPred[:iendData],qnt0[:iendData],qnt1[:iendData],color=cmap1(normalize(1-midPt)),alpha=alph,zorder=1)
            pl2 = ax.fill_between(datesPred[iendData-1:],qnt0[iendData-1:],qnt1[iendData-1:],color=cmap2(normalize(1-midPt)),alpha=alph,zorder=1)
        else:
            pl2 = ax.fill_between(datesPred[:iendData],qnt0[:iendData],qnt1[:iendData],color=cmap1(normalize(midPt)),alpha=alph,zorder=1)
            pl2 = ax.fill_between(datesPred[iendData-1:],qnt0[iendData-1:],qnt1[iendData-1:],color=cmap2(normalize(midPt)),alpha=alph,zorder=1)

    for i in range(len(run_setup["ppopts"]["quantile_newcases"])):
        qnt = run_setup["ppopts"]["quantile_newcases"][i]
        ltp = run_setup["ppopts"]["linetype_newcases"][i]
        lwd = run_setup["ppopts"]["linewidth_newcases"][i]
        qntPred=np.quantile(pred,qnt,axis=0)
        maxVal = max(maxVal,qntPred.max())
        pl1=ax.plot(datesPred,np.quantile(pred,qnt,axis=0),ltp,lw=lwd)

    symtype=run_setup["ppopts"]["newcases"][0]
    ms=run_setup["ppopts"]["newcases"][1]
    if nc_new is not None:
        pl3=ax.plot(nc_new[0],nc_new[1],'ko',ms=ms+1,mfc='w',mec='k',mew=1.5)

    pl3=ax.plot(datesData,new_cases,symtype,ms=ms)

    plt.xticks(rotation=45)

    x0=parser.parse(run_setup["ppopts"]["xylim_newcases"][0])
    x1=parser.parse(run_setup["ppopts"]["xylim_newcases"][1])

    # comment next lines
    # x1=parser.parse(rawdata[-1,0])+datetime.timedelta(days=run_setup["ppopts"]["days_extra"]) 
    x1=datesPred[-1]
    ax.set_xlim([x0,x1])

    if (len(run_setup["ppopts"]["xylim_newcases"])==4):
        y0=run_setup["ppopts"]["xylim_newcases"][2]
        y1=run_setup["ppopts"]["xylim_newcases"][3]
    else:
        y0 = 0.0
        y1 = 1.1*maxVal

    # comment next lines
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

    figname=run_setup["ppopts"]["fout_newcases"]
    plt.savefig(figname+"."+run_setup["ppopts"]["figtype"])
    plt.close() 


#--------------------------------------------------------------------------------
# Plot push-forward/posterior predictive PDF
def plot_infection_curve(datesmean,infall,qntList,normalize,iendData,run_setup):

    # colormap settings
    import matplotlib as mpl
    cmap1 = mpl.cm.PuBu
    cmap2 = mpl.cm.PuRd

    fig = plt.figure(figsize=(6,6))
    ax=fig.add_axes([0.18,0.26,0.75,0.7])

    maxVal=-1.e100
    for i in range(len(qntList)-1):
        qnt0 = np.quantile(infall,qntList[i],axis=0)
        qnt1 = np.quantile(infall,qntList[i+1],axis=0)
        midPt = 0.5*(qntList[i]+qntList[i+1])
        alph = 0.8
        if qntList[i] >= 0.5:
            pl2 = ax.fill_between(datesmean[:iendData],qnt0[:iendData],qnt1[:iendData],color=cmap1(normalize(1-midPt)),alpha=alph,zorder=1)
            pl2 = ax.fill_between(datesmean[iendData-1:],qnt0[iendData-1:],qnt1[iendData-1:],color=cmap2(normalize(1-midPt)),alpha=alph,zorder=1)
        else:
            pl2 = ax.fill_between(datesmean[:iendData],qnt0[:iendData],qnt1[:iendData],color=cmap1(normalize(midPt)),alpha=alph,zorder=1)
            pl2 = ax.fill_between(datesmean[iendData-1:],qnt0[iendData-1:],qnt1[iendData-1:],color=cmap2(normalize(midPt)),alpha=alph,zorder=1)

    for i in range(len(run_setup["infopts"]["quantile_inf"])):
        qnt = run_setup["infopts"]["quantile_inf"][i]
        ltp = run_setup["infopts"]["linetype_inf"][i]
        lwd = run_setup["infopts"]["linewidth_inf"][i]
        qntPred=np.quantile(infall,qnt,axis=0)
        maxVal = max(maxVal,qntPred.max())
        pl1=ax.plot(datesmean,np.quantile(infall,qnt,axis=0),ltp,lw=lwd)

    plt.xticks(rotation=45)

    x0=run_setup["infopts"]["xylim_inf"][0]
    x1=run_setup["infopts"]["xylim_inf"][1]
    # ax.set_xlim([parser.parse(x0),parser.parse(x1)])
    ax.set_xlim([parser.parse(x0),datesmean[-1]])

    y0=run_setup["infopts"]["xylim_inf"][2]
    y1=run_setup["infopts"]["xylim_inf"][3]
    #ax.set_ylim([y0,y1])
    ax.set_ylim([y0,1.1*maxVal])

    if "xy_scale" in run_setup["infopts"]:
        # ax.set_xscale(run_setup["infopts"]["xy_scale"][0])
        ax.set_yscale(run_setup["infopts"]["xy_scale"][1])
    else:
        ax.set_yscale("log")
    
    xlbs=run_setup["infopts"]["xylbl_inf"][0]
    xlbf=run_setup["infopts"]["xylbl_inf"][1]
    ax.set_xlabel(xlbs,fontsize=xlbf)
    ylbs=run_setup["infopts"]["xylbl_inf"][2]
    ylbf=run_setup["infopts"]["xylbl_inf"][3]
    ax.set_ylabel(ylbs,fontsize=ylbf)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(run_setup["infopts"]["xyticklbl_inf"][0])

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(run_setup["infopts"]["xyticklbl_inf"][1])

    figname=run_setup["infopts"]["fout_inf"]
    plt.savefig(figname+"."+run_setup["infopts"]["figtype"])
    plt.close()
    #plt.show()

