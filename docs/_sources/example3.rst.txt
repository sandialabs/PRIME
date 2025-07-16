Multi-region Model
==================

The following examples demonstrate the multi-region capability of PRIME. We consider 
two- and three-region cases as illustrated in the figure below

.. figure:: ./figures/counties_BSFV.pdf
    :width: 75 %
    :name: countiesNM

    Adjacent counties in New Mexico used to illustrate the multi-region capabilities in PRIME.

Two Regions Example
-------------------

The first example considers daily confirmed case data in the state of New Mexico 
up to September 15 26th, 2020.  

Problem Setup
~~~~~~~~~~~~~

The setup for running multiple regions is similar with the 1-region model. The snippet below
provides details on the setup for the Bernalillo and Santa Fe counties of New Mexico 
inferred jointly. The entries for "spllo" (parameters' low bounds), "splhi" (parameters' high bounds),
"cini" (initial guess), "cvini" (diagonal entries for the MCMC proposal covariance) now contain
10 entries (4 parameters for each region, 
:math:`\log\sigma_a`, :math:`\log\sigma_m`, :math:`\log\tau_{\phi}^2`, :math:`\lambda`)

.. code-block:: JSON

    "regioninfo":{
        "count_data":["BERNALILLO.dat","SANTA_FE.dat"],
        "population_data":[0.670, 0.145],
        "region_tag":["BERNALILLO","SANTA_FE"],
        "day0":"2020-06-10",
        "running_avg_obs":7
    },
    "mcmcopts":{
        "logfile":"logmcmcNM.txt",
        "nsteps":2000000,
        "nfinal":10000000,
        "gamma":0.2,
        "spllo":[-30, 0.00001,0.1, 0.1,
                 -30, 0.00001,0.1, 0.1,
                 -30.0, 0.0001, -30.0, -20],
        "splhi":[ 40, 0.500,50.0,400.0,
                  40, 0.500,50.0,400.0,
                  10.0, 0.99, 10.0, 1.0],
        "cini":[-4.9, 0.0067, 3.97, 11.74,
                -15.8, 0.006, 4.65, 16.27,
                 3.25, 0.6, -13.02, -6.8],
        "cvini":[1e0, 6e-8, 1e-1, 1e-1,
                 1e0, 6e-8, 1e-1, 1e-1,
                 1e-2, 1e-4, 1e-2, 1e-2]
    }

Two-region Inference Results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The posterior predictive results for the case described above are shown in  :numref:`B2reg` for Bernalillo County 
and in :numref:`SF2reg` for the Santa Fe County.

.. list-table:: 

    * - .. figure:: ./figures/B_forecast_2reg.pdf 
           :width: 100 %
           :name: B2reg

           2-region inference: Bernalillo County

      - .. figure:: ./figures/SF_forecast_2reg.pdf
           :width: 100 %
           :name: SF2reg

           2-region inference: Santa Fe County


Three Regions Example
---------------------

The setup for the three-region example follows the same pattern as the one for two regions. 

.. list-table:: 

    * - .. figure:: ./figures/B_forecast_3reg.pdf 
           :width: 100 %
           :name: B3reg

           3-region inference: Bernalillo County

      - .. figure:: ./figures/SF_forecast_3reg.pdf
           :width: 100 %
           :name: SF3reg

           3-region inference: Santa Fe County

      - .. figure:: ./figures/V_forecast_3reg.pdf
           :width: 100 %
           :name: V3reg

           3-region inference: Valencia County


JSON Input Files
----------------

Two Regions
~~~~~~~~~~~

.. code-block:: JSON

    {
        "regioninfo":{
            "count_data":["BERNALILLO.dat","SANTA_FE.dat"],
            "population_data":[0.670,0.145],
            "region_tag":["BERNALILLO","SANTA_FE"],
            "day0":"2020-06-10",
            "running_avg_obs":7
        },
        "mcmcopts":{
            "logfile":"logmcmcNM.txt",
            "nsteps":2000000,
            "nfinal":10000000,
            "gamma":0.2,
            "spllo":[-30,0.00001,0.1,0.1,-30,0.00001,0.1,0.1,-30.0,0.0001,-30.0,-20],
            "splhi":[ 40,0.500,50.0,400.0,40,0.500,50.0,400.0,10.0,0.99,10.0,1.0],
            "cini":[-4.9,0.0067,3.97,11.74,-15.8,0.006,4.65,16.27,3.25,0.6,-13.02,-6.8],
            "cvini":[1e0,6e-8,1e-1,1e-1,1e0,6e-8,1e-1,1e-1,1e-2,1e-4,1e-2,1e-2]
        },
        "bayesmod":{
            "prior_types":["g","u","u","u","g","u","u","u","lgm","u","u","u"],
            "prior_info":[[0,10],[0,1],[0,1],[0,1],[0,10],[0,1],[0,1],[0,1],[10,2],[0,1],[0,1],[0,1]],
            "error_model_type": "addMult",
            "lpf_type":"gaussian"
        },
        "modelopts":{
            "num_waves":1,
            "useconv":1,
            "incubation_median":5.1,
            "incubation_sigma":0.418,
            "incubation_025":2.2,
            "incubation_975":11.5,
            "incubation_model":"lognormal",
            "incubation_type":"stochastic"
        },
        "ppopts":{
            "nstart":500000,
            "nsamples":5000,
            "days_extra":14,
            "runmodel":1,
            "postpred":1,
            "quantile_newcases":[0.025,0.25,0.5,0.75,0.975],
            "linetype_newcases":["b--","g-","r-","g-","b--"],
            "linewidth_newcases":[3,2,3,2,3],
            "fillbetw_newcases":[[0.25,0.5,"g",0.4],[0.5,0.75,"g",0.4]],
            "xylim_newcases":["2020-06-01","2020-07-15",0,120],
            "xylbl_newcases":["Date",16,"Reported New Cases on Date",16],
            "xyticklbl_newcases":[14,14],
            "newcases":["ko",6],
            "figtype":"pdf",
            "fpredout":"_forecast",
            "fout_newcases":"_forecast"
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
            "xylim_inf":["2020-06-01","2020-05-01",1,200],
            "xylbl_inf":["Date",16,"Infection Rate [ppl/day]",16],
            "xyticklbl_inf":[14,14],
            "newcases":["ko",6],
            "figtype":"pdf",
            "finfout":"_infection",
            "fout_inf":"_infection"
        },
        "csvout":{
            "nskip":50,
            "finfcurve":"_infection",
            "fnewcases":"_forecast",
            "qlist":[0.025,0.25,0.5,0.75,0.975]
        }
    }

Three Regions
~~~~~~~~~~~~~

.. code-block:: JSON

    {
        "regioninfo":{
            "count_data":["BERNALILLO.dat","SANTA_FE.dat","VALENCIA.dat"],
            "population_data":[0.670,0.145,0.077],
            "region_tag":["BERNALILLO","SANTA_FE","VALENCIA"],
            "day0":"2020-06-10",
            "running_avg_obs":7
        },
        "mcmcopts":{
            "logfile":"logmcmcNM.txt",
            "nsteps":2000000,
            "nfinal":10000000,
            "gamma":0.2,
            "spllo":[-30,0.00001,0.1,0.1,-30,0.00001,0.1,0.1,-30,0.00001,0.1,0.1,-30.0,0.0001,-30.0,-20],
            "splhi":[ 40,0.500,50.0,400.0,40,0.500,50.0,400.0,30,0.500,30.0,400.0,10.0,0.99,10.0,2.0],
            "cini":[-4.9,0.0067,3.97,11.74,-5.8,0.006,4.65,16.27,-5.33,0.0066,3.92,12.26,1,0.5,-13.02,-6.8],
            "cvini":[1e0,6e-8,1e-1,1e-1,1e0,6e-8,1e-1,1e-1,1e0,3e-8,1e-2,1e-1,1e-2,1e-4,1e-2,1e-2]
        },
        "bayesmod":{
            "prior_types":["g","u","u","u","g","u","u","u","g","u","u","u","lgm","u","u","u"],
            "prior_info":[[0,10],[0,1],[0,1],[0,1],[0,10],[0,1],[0,1],[0,1],[0,10],[0,1],[0,1],[0,1],[10,2],[0,1],[0,1],[0,1]],
            "error_model_type": "addMult",
            "lpf_type":"gaussian"
        },
        "modelopts":{
            "num_waves":1,
            "useconv":1,
            "incubation_median":5.1,
            "incubation_sigma":0.418,
            "incubation_025":2.2,
            "incubation_975":11.5,
            "incubation_model":"lognormal",
            "incubation_type":"stochastic"
        },
        "ppopts":{
            "nstart":500000,
            "nsamples":5000,
            "days_extra":14,
            "runmodel":1,
            "postpred":1,
            "quantile_newcases":[0.025,0.25,0.5,0.75,0.975],
            "linetype_newcases":["b--","g-","r-","g-","b--"],
            "linewidth_newcases":[3,2,3,2,3],
            "fillbetw_newcases":[[0.25,0.5,"g",0.4],[0.5,0.75,"g",0.4]],
            "xylim_newcases":["2020-06-01","2020-07-15",0,120],
            "xylbl_newcases":["Date",16,"Reported New Cases on Date",16],
            "xyticklbl_newcases":[14,14],
            "newdata":"plus10.dat",
            "newcases":["ko",6],
            "figtype":"pdf",
            "fpredout":"_forecast",
            "fout_newcases":"_forecast"
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
            "xylim_inf":["2020-06-01","2020-05-01",1,200],
            "xylbl_inf":["Date",16,"Infection Rate [ppl/day]",16],
            "xyticklbl_inf":[14,14],
            "newcases":["ko",6],
            "xy_scale":["linear","linear"],
            "figtype":"pdf",
            "finfout":"_infection",
            "fout_inf":"_infection"
        },
        "csvout":{
            "nskip":50,
            "finfcurve":"_infection",
            "fnewcases":"_forecast",
            "qlist":[0.025,0.25,0.5,0.75,0.975]
        }
    }
