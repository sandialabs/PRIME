Model description
=================

.. |t0| replace:: :math:`t_0`
.. |thj| replace:: :math:`\theta_j`
.. |kj| replace:: :math:`k_j`
.. |dtj| replace:: :math:`\Delta t_j`
.. |pLN| replace:: :math:`f_{LN}`
.. |cLN| replace:: :math:`F_{LN}`

PRIME implements an epidemiological model to characterize and forecast the rate at
which people turn symptomatic from disease over time. For the purpose of this model,
we assume that once people develop symptoms, they have ready access to medical
services and can be diagnosed readily. From this perspective, these forecasts
represent a lower bound on the actual number of people that are infected with
COVID-19 as the people currently infected, but still incubating, are not accounted
for. A fraction of the population infected might also exhibit minor or no symptoms
at all and might not seek medical advice. Therefore, these cases will not be part
of patient counts released by health officials. The epidemiological model consists
of two canonical elements: an infection rate model and an incubation rate model.
One or more infection rate models are then combined through a convolution with the
incubation rate model to yield the number of cases that turn symptomatic daily.

Infection Rate Model
--------------------

The infection rate component is modeled as a Gamma distribution

.. math::
   :nowrap:

   \begin{equation}
   f_\Gamma(t;k_j,\theta_j)=\theta_j^{-k_j}t^{k_j-1}\exp(-t/\theta_j)\big/\Gamma(k_j)
   \end{equation}

with shape |kj| and scale |thj|. The choice of values for the pair (|kj|,|thj|) can 
accomodate both sharp increases in the number of infections, which would correspond 
to strained medical resources, as well as weaker gradients corresponding to a smaller 
pressure on the available medical resources. :numref:`inf-rate` show example
infection rate curves for several shape and scale parameter values.


.. figure:: ./figures/inf_rate.pdf
    :width: 90 %
    :name: inf-rate

    Infection rate models with fixed scale parameters :math:`\theta=10` (left frame) 
    and fixed shape parameter :math:`k=3` (right frame).

.. _inc-rate-section:

Incubation Rate Model
---------------------

PRIME employs a lognormal incubation distribution for COVID-19 [Lauer2020]_.
The probability density function (PDF), |pLN|, and cumulative distribution function (CDF),
|cLN|, of the lognormal distribution are given by

.. math::
   :nowrap:

   \begin{equation}
   f_{LN}(t;\mu,\sigma)=\frac{1}{t\sigma\sqrt{2}}\exp\left(-\frac{(\log t-\mu)^2}{2\sigma^2}\right),\,
   F_{LN}(t;\mu,\sigma)=\frac{1}{2}\mathrm{erfc}\left(-\frac{\log t-\mu}{\sigma\sqrt{2}}\right)
   \label{eq:inceq}
   \end{equation}

In this toolkit we model the mean :math:`\mu` as Student's t distribution with :math:`n=36` degrees of 
freedom which provided the closest agreement for the 95\% confidence interval with the data 
in [Lauer2020]_. Similarly, the standard deviation :math:`\sigma` is assumed to have a 
chi-square distribution. The resulting 95\% CIs are :math:`[1.48, 1.76]` and :math:`[0.320, 0.515]` 
for :math:`\mu` and :math:`\sigma`, respectively. The left frame in :numref:`inc-rate` shows the family of PDFs
with :math:`\mu` and :math:`\sigma` drawn from Student's t and :math:`\chi^2` distributions,
respectively. The nominal incubation PDF is shown in black in
this frame. The impact of the uncertainty in the incubation model
parameters is displayed in the right frame of this figure. For example,
7 days after infection, there is a large variability (60\%-90\%) in the
fraction of infected people that completed the incubation phase and started
displaying symptoms. This variability decreases at later times, e.g. after
10 days more then 85\% of case completed the incubation process.

.. figure:: ./figures/inc_rate.pdf
    :width: 90 %
    :name: inc-rate

    Probability density functions for the incubation model (left frame) and fraction of people for 
    which incubation ran its course after 7, 10, and 14 days respectively (right frame)

.. _single-wave:

Single Wave Model
-----------------

With these assumptions the number of people infected *and* with completed incubation period in the 
time range :math:`[t_{i-1},t_i]` can be written as a convolution between the infection rate and the 
incubation rate [Safta2020]_

.. math::
   :nowrap:

   \begin{equation}
   n_i\approx N(t_i-t_{i-1})\int_{t_0}^{t_i} f_\Gamma(\tau-t_0;k,\theta)
   f_{LN}(t_i-\tau;\mu,\sigma)d\tau
   \label{eq:symptApprox}
   \end{equation}


Here :math:`N` represents the total number of cases over the course of the epidemic, :math:`(t_i-t_{i-1})`
is typically equal to 1 day, and the time parameter |t0| represents the start of
the epidemic. In the expression above the sub-script *j* was neglected since the model contains one
wave only.


.. _multi-wave:

Multi-Wave Model
-----------------

The multiple wave model is an extension of the single wave model presented above. In the multi-wave model, 
a set of infection curves are superimposed to approximate the evolution of the epidemic that exhibits 
multiple peaks across certain regions. The resulting model is written as

.. math::
   :nowrap:

   \begin{equation}
   n_i\approx (t_i-t_{i-1}) \int_{t_0}^{t_i} \left( \sum_{j=1}^K N_j f_\Gamma(\tau-t_0-\Delta t_j;k_j,\theta_j) \right)
   f_{LN}(t_i-\tau;\mu,\sigma)d\tau  
   \label{eq:symptApproxMultiWave}
   \end{equation}

Here :math:`N_j` represents the total number of cases over the course of the *j-th* wave, and |dtj| represents 
the time shift for the *j-th* infection curve with respect to the start of the epidemic |t0|.

.. [Lauer2020] `Lauer S.A., Grantz K.H., Bi Q., Jones F.K., Zheng Q., Meredith H.R., Azman A.S., Reich N.G., Lessler J., The Incubation Period of Coronavirus Disease 2019 (COVID-19) From Publicly Reported Confirmed Cases: Estimation and Application, Annals of Internal Medicine (2020), <https://doi.org/10.7326/M20-0504>`_

.. [Safta2020] `Safta C., Ray J., and Sargsyan K., Characterization of partially observed epidemics through Bayesian inference: application to COVID-19, Computational Mechanics (2020), <https://doi.org/10.1007/s00466-020-01897-z>`_
