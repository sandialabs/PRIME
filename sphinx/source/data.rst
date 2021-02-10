Data
====

The number of people developing symptoms daily are compared to data obtained 
from several sources at the national, state, or regional levels [JHUCOVID19]_ [NYTCOVID19]_. 
We found that, for some states or regions, 
the reported daily counts exhibited a significant amount of noise. This is 
caused by variation in testing capabilities and sometimes by how data is 
aggregated from region to region over the coarse of a week. To filter the 
noise observed in daily case count data, we use of 7-day rolling averages. 
Time series of daily counts, raw data with black symbols and filtered with 
red symbols, are presented in :numref:`casedata`. We employ filtered 
data in the examples shown in this manual.

.. figure:: ./figures/FL_NM_case_data.pdf 
    :width: 90 %
    :name: casedata

    Daily confirmed cases of COVID-19 aggregated at state
    level, shown in black symbols, and the corresponding
    7-day averaged data shown with red lines and symbols. 

.. [JHUCOVID19] `COVID-19 Data Repository by the Center for Systems Science and Engineering (CSSE) at Johns Hopkins University <https://github.com/CSSEGISandData/COVID-19>`_

.. [NYTCOVID19] `Coronavirus (Covid-19) Data in the United States <https://github.com/nytimes/covid-19-data>`_

