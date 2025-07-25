Multi-Wave Model
================

The following examples demonstrate the multi-wave capability of PRIME. 

Two Waves
---------

The first example considers daily confirmed case data in the state of New Mexico 
up to August 26th, 2020, shown in :numref:`casedataNM0826`.

.. figure:: ./figures/NM_case_data_0826.pdf 
    :width: 60 %
    :name: casedataNM0826

    Daily confirmed cases of COVID-19 in New Mexico up
    to August 26th, 2020. Daily confirmed cases are shown as black 
    symbols, and the corresponding 7-day averaged data 
    shown with red lines and symbols. 

We use the two-wave model to fit this since there are two large peaks in new cases
in early May and late July 2020. 

Problem Setup
~~~~~~~~~~~~~

While the two wave model can be run independently it is highly recommended to make 
use of one-wave model results to help compute prior distributions for the first wave
in the two-wave model. This is done to i) make use of all available information and 
ii) limit the parameter range for the early wave parameters to enhance the robustness 
of the multi-wave model.   

The priors for the first wave parameters :math:`t_0`, :math:`N_1`, :math:`k_1`,
:math:`\theta_1` can be then be estimated by the mean and variance of the one-wave 
MCMC chain. This gives us an informed guess for part of the parameters, limiting the 
region of parameter space that an MCMC approach needs to explore. A smaller parameter 
range increases the chance that MCMC is able to find find regions of high likelihood 
:math:`p\left(\boldsymbol{y} | \Theta \right)`, making the model more robust. The 
second wave prior distributions should be chosen to include values that approximate 
the case data well; these distributions may need to be selected by trial and error, 
although :math:`\Delta t_2` can be estimated by visually inspecting raw daily new 
case data. 

This example uses Gaussian prior distributions for :math:`t_0` and :math:`\Delta t_2`
and uniform prior distributions for all other model parameters. For Gaussian priors 
the mean and variance are estimated from the one-wave chain. For uniform priors the 
minimum and maximum parameter are set equal to :math:`\mu_{chain} \pm 3 \sigma_{chain}`, 
where :math:`\mu_{chain}` and :math:`\sigma_{chain}` are the mean and standard 
deviation computed from the one-wave chain. Note that positive parameters 
:math:`N_1,k_1,\theta_1` are restricted to positive values; that is, the lower bound 
is the maximum of zero and :math:`\mu_{chain} - 3 \sigma_{chain}`. 

The json file for a two-wave model can be produced from a one-wave model by running
the script "prime_compute_prior_from_mcmc.py" in the directory containing the one-wave 
model json script and HDF5 file with MCMC chain data. The arguments for this script in
this case are the name of json file used to run the one-wave model and the model type 
for the desired json file, in this case "twoWave". For this example, we use
the results from example 1 to generate the two-wave json file (included in the 
"JSON Input Files" section at the end of this example).  

It is recommended to experiment with the prior distributions for the 2nd wave 
parameters :math:`\Delta t_2`, :math:`N_2`, :math:`k_2`, and :math:`\theta_2`. 
The time between the first and second wave, :math:`\Delta t_2`, can be estimated 
by the number of days between the first large increase in cases and the second 
large increase. For this example, we use a wide prior distribution for 
:math:`\Delta t_2`, centered on 60 days with a standard deviation of 15 days. 

The two-wave model is also run with the "prime_run.py" script. 

New Case Forecast Results
~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: ./figures/NM_newcases_2wave_0826.pdf 
    :width: 60 %
    :name: newcasesNM0826

    Two-wave forecast for New Mexico on August 26th, 2020. Symbols and lines are
    the same as in example 1.

The two-wave forecast results are also computed by running the "prime_compute_epi_inf_curves.py"
in the same directory as the run. The postprocessing sections of the json script are
the same for multiple waves. Please refer to example 1 for more details on postprocessing.

Three Waves
-----------

The first example considers daily confirmed case data in the state of New Mexico 
up to November 10th, 2020, shown in :numref:`casedataNM1110`.

.. figure:: ./figures/NM_case_data_1110.pdf 
    :width: 60 %
    :name: casedataNM1110

    Daily confirmed cases of COVID-19 in New Mexico up
    to November 10th, 2020. Daily confirmed cases are shown as black 
    symbols, and the corresponding 7-day averaged data 
    shown with red lines and symbols. 

We use the three-wave model to fit this since there are three large peaks in new cases
in early May, late July 2020, and cases are rising rapidly in November 2020. 

Problem Setup
~~~~~~~~~~~~~

Like the two-wave model, the three-wave model can be run independently but it is highly 
recommended to make use of two-wave model results to help compute prior distributions 
for the first and second waves in the three-wave model. The priors for the parameters 
:math:`t_0`, :math:`N_1`, :math:`k_1`, :math:`\theta_1`, :math:`\Delta t_2`, :math:`N_2`,
:math:`k_2`, and :math:`\theta_2` are estimated from the two-wave model MCMC chain.

The json file for a three-wave model can be produced from a two-wave model by running
the script "prime_compute_prior_from_mcmc.py" in the directory containing the two-wave 
model json script and HDF5 file with MCMC chain data. The arguments for this script in
this case are the name of json file used to run the two-wave model and the model type 
for the desired json file, in this case "threeWave". For this example, we use
the two-wave results from the previous section to generate the three-wave json file 
(included in the "JSON Input Files" section at the end of this example).  

It is recommended to experiment with the prior distributions for the 3rd wave 
parameters :math:`\Delta t_3`, :math:`N_3`, :math:`k_3`, and :math:`\theta_3`. 
The time between the first and third wave, :math:`\Delta t_3`, can be estimated 
by the number of days between the first large increase in cases and the third
large increase. For this example, we use a wide prior distribution for 
:math:`\Delta t_3`, centered on 180 days with a standard deviation of 15 days. 

The three-wave model is also run with the "prime_run.py" script. 

New Case Forecast Results
~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: ./figures/NM_newcases_3wave_1110.pdf 
    :width: 60 %
    :name: newcasesNM1110

    Three-wave forecast for New Mexico on November 10th, 2020. Symbols and lines are
    the same as in example 1.

The three-wave forecast results are also computed by running the "prime_compute_epi_inf_curves.py"
in the same directory as the run. The postprocessing sections of the json script are
the same for multiple waves. Please refer to example 1 for more details on postprocessing.

JSON Input Files
----------------

Two Wave
~~~~~~~~

.. code-block:: JSON

    {
        "regioninfo":{
            "count_data":["NM.dat"],
            "population_data":[2.113],
            "region_tag":["NM"],
            "day0":"2020-03-01",
            "running_avg_obs":7
        },
        "mcmcopts": {
            "logfile": "logmcmcNM.txt",
            "nsteps": 1000000,
            "nfinal": 10000000,
            "gamma": 0.7,
            "spllo": [-2.87,0.0045,3.49,12.32,20,0.0002,0.1,0.1,-30.0,-20],
            "splhi": [2.969,0.0087,5.36,26.17,120,0.5,30.0,400.0,10.0,1.0],
            "cini": [0.0458,0.0066,4.42,19.24,80,0.006,6.0,12,-14,-2],
            "cvini": [0.9494,2.1586e-06,0.097,4,4,1.0e-6,0.01,0.01,0.01,0.01]
        },
        "bayesmod": {
            "prior_types": ["g","u","u","u","g","u","u","u","u","u"],
            "prior_info": [[0.0458,0.974],[0,1],[0,1],[0,1],[60,15.0],[0,1],[0,1],[0,1],[0,1],[0,1]],
            "error_model_type": "addMult",
            "lpf_type":"gaussian"
        },
        "modelopts":{
            "num_waves":2,
            "useconv":1,
            "incubation_median":5.1,
            "incubation_sigma":0.418,
            "incubation_025":2.2,
            "incubation_975":11.5,
            "incubation_model":"lognormal",
            "incubation_type":"stochastic"
        },
        "ppopts": {
            "nstart": 200000,
            "nsamples": 1000,
            "days_extra": 10,
            "runmodel": 1,
            "postpred": 1,
            "quantile_newcases": [0.025,0.25,0.5,0.75,0.975],
            "linetype_newcases": ["b--","g-","r-","g-","b--"],
            "linewidth_newcases": [3,2,3,2,3],
            "fillbetw_newcases": [[0.25,0.5,"g",0.4],[0.5,0.75,"g",0.4]],
            "xylim_newcases": ["2020-03-01","2020-04-15",0,500],
            "xylbl_newcases": ["Date",16,"Reported New Cases on Date",16],
            "xyticklbl_newcases": [14,14],
            "newcases": ["ko",6],
            "figtype": "pdf",
            "fpredout": "NM_epidemic_curve",
            "fout_newcases": "NM_epidemic_curve"
        },
        "infopts": {
            "inftype": "gamma",
            "ndays": 180,
            "runmodel": 1,
            "postpred": 1,
            "quantile_inf": [0.025,0.25,0.5,0.75,0.975],
            "linetype_inf": ["b--","g-","r-","g-","b--"],
            "linewidth_inf": [3,2,3,2,3],
            "fillbetw_inf": [[0.25,0.5,"g",0.4],[0.5,0.75,"g",0.4]],
            "xylim_inf": ["2020-03-01","2020-05-01",10,1000],
            "xylbl_inf": ["Date",16,"Infection Rate [ppl/day]",16],
            "xyticklbl_inf": [14,14],
            "newcases": ["ko",6],
            "figtype": "pdf",
            "finfout": "NM_infection_curve",
            "fout_inf": "NM_infection_curve"
        },
        "csvout": {
            "nskip": 100,
            "finfcurve": "NM_infection_curve",
            "fnewcases": "NM_epidemic_curve",
            "qlist": [0.025,0.25,0.5,0.75,0.975]
        }
    }
 
Three Wave
~~~~~~~~~~

.. code-block:: JSON

    {
        "regioninfo":{
            "count_data":["NM.dat"],
            "population_data":[2.113],
            "region_tag":["NM"],
            "day0":"2020-03-01",
            "running_avg_obs":7
        },
        "mcmcopts": {
            "logfile": "logmcmcNM.txt",
            "nsteps": 1000000,
            "nfinal": 10000000,
            "gamma": 0.7,
            "spllo": [-1.5199,0.00636,3.71538,17.571,70.818,0.0051,0.0,4.976,140,0.0002,0.1,0.1,-30.0,-20],
            "splhi": [1.8198,0.0077,4.6929,24.5749,112.5741,0.00638,11.2084,15.4422,220,0.5,30.0,400.0,10.0,1.0],
            "cini": [0.1499,0.007,4.204,21.0729,91.69619,0.00574,5.5528,10.2094,180,0.01,6.0,20.2,-10,0.1],
            "cvini": [0.3098,2.1092e-07,0.0265,1.3626,48.432,2.0496e-07,3.5539,3.0425,225,0.001,0.01,0.01,0.01,0.01]
        },
        "bayesmod": {
            "prior_types": ["g","u","u","u","g","u","u","u","g","u","u","u","u","u"],
            "prior_info": [
                [0.1499487243752924,0.5566124720557659],[0,1],[0,1],[0,1],
                [91.69619492513007,6.959308502135749],[0,1],[0,1],[0,1],
                [180,15.0],[0,1],[0,1],[0,1],
                [0,1],[0,1]
            ],
            "error_model_type": "addMult",
            "lpf_type":"gaussian"
        },
        "modelopts":{
            "num_waves":3,
            "useconv":1,
            "incubation_median":5.1,
            "incubation_sigma":0.418,
            "incubation_025":2.2,
            "incubation_975":11.5,
            "incubation_model":"lognormal",
            "incubation_type":"stochastic"
        },
        "ppopts": {
            "nstart": 100000,
            "nsamples": 1000,
            "days_extra": 10,
            "runmodel": 1,
            "postpred": 1,
            "newdata":"NM_1124.dat",
            "quantile_newcases": [0.025,0.25,0.5,0.75,0.975],
            "linetype_newcases": ["b--","g-","r-","g-","b--"],
            "linewidth_newcases": [3,2,3,2,3],
            "fillbetw_newcases": [[0.25,0.5,"g",0.4],[0.5,0.75,"g",0.4]],
            "xylim_newcases": ["2020-03-01","2020-04-15",0,1800],
            "xylbl_newcases": ["Date",16,"Reported New Cases on Date",16],
            "xyticklbl_newcases": [14,14],
            "newcases": ["ko",6],
            "figtype": "pdf",
            "fpredout": "NM_epidemic_curve",
            "fout_newcases": "NM_epidemic_curve"
        },
        "infopts": {
            "inftype": "gamma",
            "ndays": 270,
            "runmodel": 1,
            "postpred": 1,
            "quantile_inf": [0.025,0.25,0.5,0.75,0.975],
            "linetype_inf": ["b--","g-","r-","g-","b--"],
            "linewidth_inf": [3,2,3,2,3],
            "fillbetw_inf": [[0.25,0.5,"g",0.4],[0.5,0.75,"g",0.4]],
            "xylim_inf": ["2020-03-01","2020-05-01",10,1000],
            "xylbl_inf": ["Date",16,"Infection Rate [ppl/day]",16],
            "xyticklbl_inf": [14,14],
            "newcases": ["ko",6],
            "figtype": "pdf",
            "finfout": "NM_infection_curve",
            "fout_inf": "NM_infection_curve"
        },
        "csvout": {
            "nskip": 100,
            "finfcurve": "NM_infection_curve",
            "fnewcases": "NM_epidemic_curve",
            "qlist": [0.025,0.25,0.5,0.75,0.975]
        }
    }
