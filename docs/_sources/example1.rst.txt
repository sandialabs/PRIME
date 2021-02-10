One-Wave Model
==============

The following example demonstrates the one-wave model for daily
confirmed cases data in the state of New Mexico up to May 13th, 2020,
shown in :numref:`casedataNM0513`. 

.. figure:: ./figures/NM_case_data_0513.pdf 
    :width: 75 %
    :name: casedataNM0513

    Daily confirmed cases of COVID-19 in New Mexico up
    to May 13th, 2020. Daily confirmed cases are shown as black 
    symbols, and the corresponding 7-day averaged data 
    shown with red lines and symbols. 

Problem Setup
~~~~~~~~~~~~~

All the data needed by PRIME is contained in a single json file. 
The first block of the json file, "regioninfo", should contain information about the data: 

.. code-block:: JSON

   {
       "regionname":"NM",
       "fchain":"NM_mcmc.h5",
       "day0":"2020-03-01",
       "running_avg_obs":7
   }

The "regionname" provides the name of the file contain daily confirmed case data. In this case,
it is "NM.dat". This file should contain two columns with dates and daily confirmed cases, respectively: 

.. code-block:: RST

   2020-03-11 4
   2020-03-12 2
   2020-03-13 4
   2020-03-14 3
   2020-03-15 4
   2020-03-16 4
   2020-03-17 2
   2020-03-18 5
   2020-03-19 7
   2020-03-20 8
   2020-03-21 14
   2020-03-22 8
   2020-03-23 18
   2020-03-24 17

"fchain" is the name of an HDF5 file containing the MCMC chain along with other useful metadata. 
The "day0" field specifies the day with the index 0; this is important for setting the prior distribution
for :math:`t_0`. The "running_avg_obs" field sets the number of days to compute a running average over, in 
this case 7, as shown in :numref:`casedataNM0513`. 

The second block of the json file sets options for the model and MCMC:

.. code-block:: JSON

   {
       "model_type": "oneWave",
       "error_model_type": "addMult",
       "logfile":"logmcmcNM.txt",
       "nsteps":1000000,
       "nfinal":10000000,
       "useconv":1,
       "incubation_type":"uncertain",
       "gamma":0.2,
       "spllo":[-10,0.0002,0.1, 0.1, 0.0,-20],
       "splhi":[10,0.500,30.0,400.0,10.0,1.0],
       "cini":[0,0.02,6.0,20.2,3.00,0.1],
       "cvini":[0.04,0.001,0.01,0.01,0.01,0.01]
   }

The first two inputs specify how many waves the model has and the error model. 
In this example, we use a single infection wave ("oneWave") with the additive and multiplicative error 
models ("addMult"). The number of steps in the MCMC chain is set by "nsteps" and is 1000000 in this case.
The "useconv" option determines if the integrals over probability distributions should be used (1=on). 
The incubation model is set by "incubation_type", which is set to "uncertain" in this case to model the 
incubation rate as a random variable instead of a fixed value. 

The lists "spllo" and "splhi" contain minimum and maximum values, respectively, of model parameters that MCMC 
can sample. This overrides the bounds of the prior distributions. The lists "cini" and "cvini" contain
initial guesses for the mean and variance of each model parameter. 

The prior distributions are specified in the "bayesmod" block:

.. code-block:: JSON

   {
       "prior_types":["g","u","u","u","u","u"],
       "prior_info":[[0,1],[0,1],[0,1],[0,1],[0,1],[0,1]]
   }

The list "prior_types" contains the type of distribution used for each prior. In this case, we are using a 
Gaussian distribution for the first model parameter, :math:`t_0` and uniform distributions for all others.
The list "prior_info" contains the mean and standard deviation of each distribution. Note that for uniform
distributions, the entry in this list is ignored; the upper and lower bounds are set using entries in "splhi"
and "spllo", respectively. 

For this case, most of the prior distributions were determined by trial and error, with the exception of the
prior for :math:`t_0`, which was set by observing when the increase in cases started and centering the prior
7-10 days before this to account for incubation time. 

Next, the properties of the incubation model are set in the "incopts" section:

.. code-block:: JSON
   
   {
       "incubation_median":5.1,
       "incubation_sigma":0.418,
       "incubation_025":2.2,
       "incubation_975":11.5
   }


These data are used by PRIME to construct a fixed or uncertain incubation rate model. 

A json file containing these section can be used to run the MCMC and output the chain,
but other sections are needed for postprocessing. To run PRIME for this case, simply call 
the "prime_run.py" script followed by the name of the json file. 

New Case Forecast Results
~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: ./figures/NM_newcases_0513.pdf 
    :width: 75 %
    :name: newcasesNM0513

    One-wave forecast for New Mexico on May 13th, 2020. Data used to calibrate the 
    epidemiological model are indicated by filled black circles.
    The shaded color region illustrates either the posterior-predictive distribution 
    with darker colors near the median and
    lighter colors near the low and high quantile values. The blue colors correspond
    to the hindcast dates and red colors to forecasts. The inter-quartile range is
    marked with green lines and the 95\% confidence interval with dashed lines.
    The plot also shows data collected at a later time, with open circles,
    to check the agreement between the forecast and the observed number of cases
    after the model has been calibrated.

Forecast results can be computed by running the postprocessor script "prime_compute_epi_inf_curves.py" 
in the same directory as the run. This script requires several additional sections in the 
json file. Firstly, the "ppopts" section contains the information needed to plot the new
case forecast presented in :numref:`newcasesNM0513`. 

.. code-block:: JSON
   
   {
        "nstart":100000,
        "nsamples":1000,
        "days_extra":10,
        "runmodel":1,
        "postpred":1,
        "newdata": "NM_future.dat",
        "quantile_newcases":[0.025,0.25,0.5,0.75,0.975],
        "linetype_newcases":["b--","g-","r-","g-","b--"],
        "linewidth_newcases":[3,2,3,2,3],
        "fillbetw_newcases":[[0.25,0.5,"g",0.4],[0.5,0.75,"g",0.4]],
        "xylim_newcases":["2020-03-01","2020-04-15",0,300],
        "xylbl_newcases":["Date",16,"Reported New Cases on Date",16],
        "xyticklbl_newcases":[14,14],
        "newcases":["ko",6],
        "figtype":"pdf",
        "fpredout":"NM_epidemic_curve",
        "fout_newcases":"NM_epidemic_curve"
    }

The portion of the MCMC chain used to generate the plot is set by "nstart" and "nsamples". 
"nstart" sets the starting index for the portion of the chain used for postprocessing. 
"nsamples" sets the number of entries in the chain (after index "nstart") to be 
sampled uniformly for postprocessing.

"days_extra" sets how many days out to compute the forecast, in this case 10 days, or until
May 23rd, 2020. 

"runmodel" determines whether or not to run the model for each chain sample to compute new
cases or to read the new case data from the HDF5 file whose name is specified by the "fpredout"
entry in this block. This option should be set to 1 the first time that "prime_compute_epi_inf_curves.py" is
run, but can be set to 0 for subsequent runs, for example if one wants to regenerate a plot. 

"postpred" is set to 1 to plot the posterior predictive and 0 to plot the push forward PDF. 

"newdata" contains the name of an ascii file with future case data. In this case it contains
case data from May 14th, 2020 onwards. 

The next 8 entires in the "ppopts" block above correspond to plot settings. Many of them 
concern the quantile curve plots, starting with which quantiles to show ("quantile_newcases"), 
and the corresponding line color/style ("linetype_newcases"), line width
("linewidth_newcases"), and the color fill between lines ("fillbetw_newcases"). 

Plot limits, labels, and tick font sizes can be set with "xylim_newcases", "xylbl_newcases",
and "xyticklbl_newcases", respectively. Finally, "newcases" contains the a list with the 
color/symbol followed by the symbol size for the daily new case data used for the forecast. 

Finally, "figtype" sets the file format for the forecast plot to written to, "fpredout" contrains 
the name of an HDF5 file containing the data shown in the forecast plot, and "fout_newcases" is
the name of the forecast plot file. The script "prime_compute_epi_inf_curves.py" adds prefixes to indicate
which error models are used and if the posterior predictive is plotted. This means that the figure
will be written to "NM_newcases_amN_pp.pdf" for our example. 

Infection Rate Prediction Results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: ./figures/NM_infcurve_0513.pdf 
    :width: 75 %
    :name: infcurveNM0513

    One-wave infection rate curve forecast for New Mexico on May 13th, 2020. The shaded color 
    region illustrates either the posterior-predictive distribution with darker colors near 
    the median and lighter colors near the low and high quantile values. The blue colors correspond
    to the hindcast dates and red colors to forecasts. The inter-quartile range is
    marked with green lines and the 95\% confidence interval with dashed lines.


To plot the infection rate curve as presented in :numref:`infcurveNM0513`, an "infopts" section is
needed in the json file:

.. code-block:: JSON

   {
       "inftype":"gamma",
       "ndays":180,
       "runmodel":1,
       "postpred":1,
       "quantile_inf":[0.025,0.25,0.5,0.75,0.975],
       "linetype_inf":["b--","g-","r-","g-","b--"],
       "linewidth_inf":[3,2,3,2,3],
       "fillbetw_inf":[[0.25,0.5,"g",0.4],[0.5,0.75,"g",0.4]],
       "xylim_inf":["2020-03-01","2020-05-01",10,1000],
       "xylbl_inf":["Date",16,"Infection Rate [ppl/day]",16],
       "xyticklbl_inf":[14,14],
       "newcases":["ko",6],
       "figtype":"pdf",
       "finfout":"NM_infection_curve",
       "fout_inf":"NM_infection_curve"
   }

Here, "inftype" sets the infection rate curve type, in this case it is a gamma distribution. 
"runmodel" and "postpred" are the same as in the "ppopts" section. The other entries correspond
to the same plot format and file format/name settings in the "ppopts" section. 

Finally, new case and infection rate data can be written out in CSV format if a "csvout" section
is included in the json file:

.. code-block:: JSON

   {
       "nskip":100,
       "finfcurve":"NM_infection_curve",
       "fnewcases":"NM_epidemic_curve",
       "qlist":[0.025,0.25,0.5,0.75,0.975]
   }

Each row of the CSV files contain data corresponding to each date for which case data is availible 
along with the dates for which a forecast is availible. In this case, data from early March to May 23rd,
2020 is included. The data on each row includes the date, forecast quantile(s), and individual samples 
from the MCMC chain. The new case file also contains the reported daily new cases in the last column for
all dates in which it is availible. In this case, daily new cases data is included up to May 13th.  

"nskip" sets the sampling frequency for the MCMC chain. In this example, every 100th chain sample is 
included. Recall that the number of chain samples is set to 1000 in the "ppopts" section under "nsamples". 
This means that the csv files will include 10 samples in this example. 

"finfcurve" and "fnewcases" specify the file names with infection rate and new cases data, 
respectively. 

Finally, "qlist" specifies the quantiles for which to output data. 

JSON Input File
~~~~~~~~~~~~~~~

The enitre json file is included below for completeness:

.. code-block:: JSON

   {
       "regioninfo":{
           "regionname":"NM",
           "fchain":"NM_mcmc.h5",
           "day0":"2020-03-01",
           "running_avg_obs":7
       },
       "mcmcopts":{
           "model_type": "oneWave",
           "error_model_type": "addMult",
           "logfile":"logmcmcNM.txt",
           "nsteps":1000000,
           "nfinal":10000000,
           "useconv":1,
           "incubation_type":"uncertain",
           "gamma":0.2,
           "spllo":[-10,0.0002,0.1, 0.1, 0.0,-20],
           "splhi":[10,0.500,30.0,400.0,10.0,1.0],
           "cini":[0,0.02,6.0,20.2,3.00,0.1],
           "cvini":[0.04,0.001,0.01,0.01,0.01,0.01]
       },
       "bayesmod":{
           "prior_types":["g","u","u","u","u","u"],
           "prior_info":[[0,1],[0,1],[0,1],[0,1],[0,1],[0,1]]
       },
       "incopts":{
           "incubation_median":5.1,
           "incubation_sigma":0.418,
           "incubation_025":2.2,
           "incubation_975":11.5
       },
       "ppopts":{
           "nstart":100000,
           "nsamples":10000,
           "days_extra":10,
           "runmodel":1,
           "postpred":1,
           "newdata": "NM_future.dat",
           "quantile_newcases":[0.025,0.25,0.5,0.75,0.975],
           "linetype_newcases":["b--","g-","r-","g-","b--"],
           "linewidth_newcases":[3,2,3,2,3],
           "fillbetw_newcases":[[0.25,0.5,"g",0.4],[0.5,0.75,"g",0.4]],
           "xylim_newcases":["2020-03-01","2020-04-15",0,300],
           "xylbl_newcases":["Date",16,"Reported New Cases on Date",16],
           "xyticklbl_newcases":[14,14],
           "newcases":["ko",6],
           "figtype":"pdf",
           "fpredout":"NM_epidemic_curve",
           "fout_newcases":"NM_epidemic_curve"
       },
       "infopts":{
           "inftype":"gamma",
           "ndays":180,
           "runmodel":1,
           "postpred":1,
           "quantile_inf":[0.025,0.25,0.5,0.75,0.975],
           "linetype_inf":["b--","g-","r-","g-","b--"],
           "linewidth_inf":[3,2,3,2,3],
           "fillbetw_inf":[[0.25,0.5,"g",0.4],[0.5,0.75,"g",0.4]],
           "xylim_inf":["2020-03-01","2020-05-01",10,1000],
           "xylbl_inf":["Date",16,"Infection Rate [ppl/day]",16],
           "xyticklbl_inf":[14,14],
           "newcases":["ko",6],
           "figtype":"pdf",
           "finfout":"NM_infection_curve",
           "fout_inf":"NM_infection_curve"
       },
       "csvout":{
           "nskip":100,
           "finfcurve":"NM_infection_curve",
           "fnewcases":"NM_epidemic_curve",
           "qlist":[0.025,0.25,0.5,0.75,0.975]
       }
   }
