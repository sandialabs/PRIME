About
=====

PRIME is a modeling framework designed for the "real-time'" characterization and forecasting of
partially observed epidemics. Characterization is the estimation of infection spread
parameters using daily counts of symptomatic patients. The method is designed to help guide
medical resource allocation in the early epoch of the outbreak. The estimation problem is posed
as one of Bayesian inference and solved using a Markov Chain Monte Carlo technique. The framework
can accommodate multiple epidemic waves and can help identify different disease dynamics at the
regional, state, and country levels. We include examples using publicly available COVID-19 data.

Team
----
PRIME was written and developed by Cosmin Safta, Jaideep Ray, Patrick Blonigan, and Kamaljit Chowdhary. 
We would like to acknowledge helpful suggestions made by several colleagues:
Erin Acquesta, Thomas Catanach, Bert Debusschere, Sean DeRosa, Pat Finley, Edgar Galvan, 
Gianluca Geraci, John D. Jakeman, Mohammad Khalil, Khachik Sargsyan, and Teresa Portone.

Acknowledgments
---------------
Sandia National Laboratories is a multi-mission laboratory managed and operated by
National Technology and Engineering Solutions of Sandia, LLC., a wholly owned subsidiary of
Honeywell International, Inc., for the U.S. Department of Energy's National Nuclear
Security Administration under contract DE-NA-0003525.
This work was funded in part by the Laboratory Directed Research \& Development (LDRD)
program at Sandia National Laboratories and by the US Department of Energy (DOE) Office of
Science through the National Virtual Biotechnology Laboratory, a consortium of national
laboratories (Argonne, Los Alamos, Oak Ridge, and Sandia) focused on responding to COVID-19,
with funding provided by the Coronavirus CARES Act. The views expressed in the
article do not necessarily represent the views of the U.S. Department of Energy
or the United States Government.

Requirements
------------

The following python packages are required for PRIME, in addition to other
default python packages.

* dateutil, h5py, matplotlib, numpy, scipy

We have tested PRIME using with python versions 3.6-3.8.