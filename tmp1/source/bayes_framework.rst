Bayesian Framework
==================

.. |by| replace:: :math:`{\boldsymbol{y}}`
.. |bn| replace:: :math:`{\boldsymbol{n}}`
.. |eps| replace:: :math:`{\epsilon}`

Given data, |by|, in the form of time-series of daily counts and the model
predictions |bn| for the number of new symptomatic counts daily for the same time
range, we employ a Bayesian framework to calibrate the epidemiological model parameters.
The discrepancy between the data and the model is written as

.. math::
   :nowrap:

   \begin{equation}
   \boldsymbol{y} = \boldsymbol{n}(\Theta)+\epsilon
   \label{eq:discrepancy}
   \end{equation}

Here, :math:`\Theta` is a vector of model parameters,
and |eps| represents the statistical discrepancy between the model and the data.
The elements of :math:`\Theta` depend on the number of epidemic waves being modeled.

.. math::
   :nowrap:

   \begin{equation}
   \Theta=\Theta^{(1)}\cup\Theta^{(2)}\cup\ldots\cup\Theta^{(K)}\cup\Theta^{(\epsilon)}
   \label{eq:param_vec_def}
   \end{equation}

where :math:`\Theta^{(j)}`  are the parameter for the *j*-th wave of infections,
*K* is the number of waves and :math:`\Theta^{(\epsilon)}`
are parameters for the error model. The parameters for each wave are given by

.. math::
   :nowrap:

   \begin{equation}
   \Theta^{(j)}=\{\Delta t_j,N_j,k_j,\theta_j\}
   \end{equation}

For the first epidemic wave :math:`\Delta t_1=t_0`.

The error model encapsulates, in this context, both errors in the observations as well
as errors due to imperfect modeling choices. The observation errors include variations
due to testing capabilities as well as errors when tests are interpreted. Values for
the vector of parameters :math:`\Theta` is estimated in the form of a multivariate PDF
via Bayes theorem

.. math::
   :nowrap:

   \begin{equation}
   p(\Theta\vert \boldsymbol{y})\propto p(\boldsymbol{y}\vert\Theta) p(\Theta)
   \label{eq:bayes}
   \end{equation}

where :math:`p(\Theta\vert\boldsymbol{y})` is the posterior distribution we are seeking after
observing the data |by|, :math:`p(\boldsymbol{y}\vert\Theta)` is the likelihood of observing
the data |by| for a particular set of values for the model parameters :math:`\Theta`,
and :math:`p(\Theta)` encapsulates any prior information available for the model parameters.
Bayesian methods are well-suited for dealing with heterogeneous sources of uncertainty,
in this case from our modeling assumptions, i.e. model and parametric uncertainties,
as well as the communicated daily counts of COVID-19 new cases, i.e. experimental errors.

Likelihood Construction
-----------------------

The library provides options for  both deterministic and stochastic formulations for 
the incubation model. In the former case the mean and standard deviation of the
incubation model are fixed at their nominal values and the model prediction
:math:`n_i` for day :math:`t_i` is a scalar value that depends on $\Theta$ only. 
In the latter case, the incubation model is stochastic with mean and standard deviation
of its natural logarithm treated as Student's t and $\chi^2$ random variables, 
respectively, as discussed in :ref:`inc-rate-section` Section. Let us denote the underlying
independent random variables by :math:`{\boldsymbol \xi}=\{\xi_\mu,\xi_\sigma,\}`.
The model prediction :math:`n_j({\boldsymbol{\xi}})` for day *j* is now a random 
variable induced by ${\bm \xi}$ plugged into :ref:`single-wave` or :ref:`multi-wave`
and :math:`{\boldsymbol n}({\boldsymbol\xi})` is a random vector.

Deterministic Incubation Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PRIME provides options for both Gaussian and negative binomial formulations for
the statistical discrepancy :math:`\epsilon` between :math:`{\boldsymbol n}` and 
:math:`{\boldsymbol y}`. In the first approach we assume :math:`\epsilon` has a zero-mean
Multivariate Normal (MVN) distribution. Given the sparsity of data, correlations across
time are currently neglected and the liklihood :math:`p({\boldsymbol y}\vert\Theta)` is computed as

.. math::
   :nowrap:

   \begin{equation}
   p({\boldsymbol y}\vert\Theta)=\prod_{i=1}^D\pi_{n_i(\Theta)}(y_i)
   =(2\pi)^{-D/2}\prod_{i=1}^D \sigma_i^{-1}
   \exp\left(-\frac{(y_i-n_i)^2}{2\sigma_i^2}\right)
   \label{eq:deticg}
   \end{equation}


with 

.. math::
   :nowrap:

   \begin{equation}
   \sigma_i=\sigma_a+\sigma_m\, n_i(\Theta)
   \end{equation}

The additive, :math:`\sigma_a`, and multiplicative, :math:`\sigma_a`, components 
of the error model :math:`\Theta^{(\epsilon)}=\{\sigma_a,\sigma_a\}` will be inferred 
jointly with the model parameters. In practive, PRIME infers the logarithm of 
these parameters to ensure they remain positive. 

The second approach assumes a negative-binomial distribution
for the discrepancy between data and model predictions. The negative-binomial distribution 
is commonly used in epidemiology to model overly dispersed data, e.g. in case 
where the standard deviation exceeds the mean [Lloyd2007]_. For this modeling choice, 
the likelihood of observing the data given a choice for the model parameters is given by

.. math::
   :nowrap:

   \begin{equation}
   p({\boldsymbol y}\vert\Theta)=\prod_{i=1}^D\pi_{n_i(\Theta)}(y_i)
   =\prod_{i=1}^D \binom{y_i+\alpha-1}{\alpha-1}
   \left(\frac{\alpha}{\alpha+n_i(\Theta)}\right)^\alpha
   \left(\frac{n_i(\Theta)}{\alpha+n_i(\Theta)}\right)^{y_i}
   \label{eq:deticnb}
   \end{equation}

where :math:`\alpha>0` is the dispersion parameter, and

.. math::
   :nowrap:

   \begin{equation}
   \binom{y_i+\alpha-1}{\alpha-1}=\frac{\Gamma(y_i+\alpha)}{\Gamma(\alpha)\Gamma(y_i+1)}
   \end{equation}

is the binomial coefficient. For simulations employing a negative binomial 
distribution of discrepancies, the logarithm of the dispersion parameter :math:`\alpha` 
will be inferred jointly with the other model parameters.

Stochastic Incubation Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the {\it stochastic incubation model} the likelihood reads as 

.. math::
   :nowrap:

   \begin{equation}
   p({\boldsymbol y}\vert\Theta)=\pi_{{\boldsymbol n}(\Theta),{\boldsymbol\xi}}
   ({\boldsymbol y})\approx\prod_{i=1}^D\pi_{n_i(\Theta),{\boldsymbol\xi}}(y_i)
   \label{eq:detfullpost}
   \end{equation}

The second expression in the right-hand side above assumes independence of the 
discrepancies between different days. Unlike the deterministic incubation model, the 
likelihood components for each day :math:`\pi_{n_i(\Theta),{\boldsymbol\xi}}(y_i)` 
are not analytically tractable anymore since they now incorporate contributions 
from :math:`{\boldsymbol \xi}`, i.e. from the variability of the parameters of the incubation model. 
One can evaluate the likelihood via kernel density estimation (KDE) by sampling 
:math:`{\boldsymbol \xi}` for each sample of :math:`\Theta`, and combining these samples 
with samples of the assumed discrepancy :math:`\epsilon`, in order to arrive at 
an estimate of :math:`\pi_{n_i(\Theta),{\boldsymbol\xi}}(y_i)`. In fact, by sampling a 
*single* value of :math:`{\boldsymbol \xi}` for each sample of :math:`\Theta`, 
one achieves an unbiased estimate of the likelihood :math:`\pi_{n_i(\Theta),{\boldsymbol\xi}}(y_i)`, 
and given the independent-component assumption, it also leads to an unbiased estimate of the 
full likelihood :math:`\pi_{{\boldsymbol n}(\Theta),{\boldsymbol\xi}}({\boldsymbol y})`. 

Posterior Distribution Sampling
-------------------------------

A Markov Chain Monte Carlo (MCMC) algorithm is used to sample from the
posterior density :math:`p(\Theta\vert\boldsymbol{y})`. MCMC is a class of techniques that
allows sampling from a posterior distribution by constructing a Markov Chain
that has the posterior as its stationary distribution. In particular, PRIME uses a
an adaptive Metropolis algorithm [Haario2001]_. Given the construction corresponding to the 
stochastic incubation model presented above, we employ the unbiased estimate of the approximate 
likelihood. This is the essence of the pseudo-marginal MCMC algorithm [Andrieu2009]_ guaranteeing 
that the accepted MCMC samples correspond to the posterior distribution.
At each MCMC step we draw a random sample :math:`\xi` from its distribution, and then 
we estimate the likelihood in a way similar to the deterministic incubation model.

:numref:`my-custom-label` displays 1D and 2D joint marginal
distributions based on two-wave model results. We use the Raftery-Lewis diagnostic [Raftery1992]_ 
to determine the number of MCMC samples required for converged statistics corresponding
to stationary posterior distributions for :math:`\Theta`. The required number of samples
is of the order :math:`o(10^5-10^6)` depending on the geographical region employed in
the inference. The resulting Effective Sample Size [Kass1998]_ varies between *8,000* and
*15,000* samples depending on each parameter which is sufficient to estimate joint
distributions for the model parameters.

.. figure:: ./figures/CA_0819_kde.pdf
    :width: 90 %
    :name: my-custom-label

    1D and 2D joint marginal distributions the components of
    :math:`\Theta=\{t_0,N_1,k_1,\theta_1,\Delta t_2,N_2,k_2,\theta_2,\log\sigma_a, \log\sigma_m\}` 
    for data from California up to 2020-08-19. Distance correlations for each pair of 
    parameters is displayed above each joint marginal distribution.


Posterior Predictive Tests
--------------------------

We will employ Bayesian posterior-predictive distributions [Lynch2004]_ to assess the 
predictive skill of the statistical model. The Bayesian posterior-predictive distribution, 
defined below is computed by marginalization of the likelihood over the posterior distribution 
of model parameters :math:`\Theta`:

.. math::
  :nowrap:

  \begin{equation}
  p_{\mathrm{pp}}\left(\boldsymbol{y}^{\mathrm{(pp)}}\vert\boldsymbol{y}\right)=\int_{\boldsymbol{\Theta}}
  p(\boldsymbol{y}^{\mathrm{(pp)}}\vert\Theta)
  p(\Theta\vert\boldsymbol{y}) d\Theta.
  \label{eq:ppd}
  \end{equation}

The posterior predictive distribution 
:math:`p_{\mathrm{pp}}\left(\boldsymbol{y}^{\mathrm{(pp)}}\vert\boldsymbol{y}\right)`
is estimated through sampling, using the parameter samples readily available 
from the MCMC exploration of the parameter space, i.e. similar to results shown 
in Fig.~\ref{fig:mcmc}. Typically we subsample the MCMC chain to about 
10-15K samples that will be used to generate posterior predictive statistics.
After the model evaluations :math:`\boldsymbol{y}=\boldsymbol{n}(\Theta)` 
are completed, we add random noise consistent
with the likelihood model settings presented in \ref{sec:lk}. The resulting
samples are used to compute summary statistics corresponding to 
:math:`p_{\mathrm{pp}}\left(\boldsymbol{y}^{\mathrm{(pp)}}\vert\boldsymbol{y}\right)`.

The posterior-predictive distribution results can be used in hindcast mode, 
to check how well the model follows the data, and for short-term
forecasts for the spread dynamics of this disease. In the hindcast
regime, the infection rate is convolved with the incubation rate model to generate
statistics for :math:`\boldsymbol{y}^{\mathrm{(pp)}}` that will be 
compared against :math:`\boldsymbol{y}`, the data used to infer the model parameters. 
The same functional form can be used to generate statistics for 
:math:`\boldsymbol{y}^{\mathrm{(pp)}}` beyond the set of dates for which 
data was available. We limit these forecasts to 7--10 days as our infection 
rate model does not count for changes in social dynamics that can significantly 
impact the epidemic over a longer timerange.


Model Selection
---------------

Quantitative comparisons between models can be made with several metrics defined in the following sections.  

AIC
~~~

The Akaike Information Criteria (AIC) [Akaike1974]_ is defined as

.. math::
   :nowrap:

    \begin{equation}
    AIC = 2 m_{\Theta} - 2 \ln (L_{max}), 
    \label{eq:AIC}
    \end{equation}

where :math:`m_{\Theta}` is the number of parameters in :math:`\Theta` and :math:`L_{max}` 
is the maximum value of the likelihood :math:`p(\boldsymbol{y}\vert\Theta)`. This is 
estimated by the maximum likelihood in the MCMC chain. Given a choice of models, the 
model with the smallest AIC value is considered to be the highest quality model. 

BIC
~~~

The Bayesian Information Criteria (BIC) [Schwarz1978]_ is defined as

.. math::
   :nowrap:

   \begin{equation}
   BIC = m_{\Theta} \ln(d) - 2 \ln (L_{max}), 
   \label{eq:BIC}
   \end{equation}

where :math:`d` is the number of observations, equal to the length of the array 
:math:`\boldsymbol{n}`. Given a choice of models, the model with the smallest BIC value 
is considered to be the highest quality model. 

CRPS
~~~~

The Continuous Ranked Probability Score (CRPS) [Gneiting2007]_ measures the difference between the CDF of the provided
data and that of the forecast/predicted data, i.e., data generated based on the posterior 
predictive distribution. It is computed by summing up marginal distributions for each day that data is available

.. math::
   :nowrap:

   \begin{equation}
   CRPS = \frac{1}{d} \sum_{j=1}^{d} \int_{-\infty}^{\infty} \left(\mathcal{F}_{pp,j}
   \left(y^{(pp)}_j| \boldsymbol{y} \right) - \mathcal{H}_{y_j}\left(y^{(pp)}_j \right) \right)^2 dy^{(pp)}_j, 
   \label{eq:crps}
   \end{equation}

where :math:`y^{(pp)}_j \equiv y^{(pp)}(t_j)` is new daily case predictions on day :math:`j` 
obtained via the posterior-predictive distribution, :math:`y_j \equiv y(t_j)` is new daily 
case data on day :math:`j` and :math:`\mathcal{F}_{pp,j}` is the 1-D marginal posterior 
predictive CDF for day :math:`j` computed using 1-D marginal posterior predictive distributions

.. math::
   :nowrap:    

   \begin{equation}
   \mathcal{F}_{pp,j}( y^{(pp)}_j | \boldsymbol{n} ) = \int_{-\infty}^{y^{(pp)}_j} p_{pp,j} 
   \left(y^{(pp)'}_j | \boldsymbol{n} \right) dy^{(pp)'}_j
   \end{equation}

where

.. math::
   :nowrap:

   \begin{equation}
   p_{pp,j} \left(y^{(pp)}_j | \boldsymbol{n} \right) = \int p_{pp}(\boldsymbol{y}^{(pp)} | \boldsymbol{y}) d\boldsymbol{y}^{(pp)}_{\sim j}
   \end{equation}

is the marginal 1-D posterior predictive density corresponding to day :math:`j`, based 
on :math:`p_{pp}\left( \boldsymbol{y}^{(pp)} | \boldsymbol{y} \right)`. Here, 
:math:`d\boldsymbol{y}^{(pp)}_{\sim j} \equiv dy^{(pp)}_1 \cdots dy^{(pp)}_{j-1} dy^{(pp)}_{j+1} \cdots dy^{(pp)}_d`. 
The CDF of the provided case data :math:`\boldsymbol{y}` is approximated as a Heaviside 
function centered at :math:`y_j`, :math:`\mathcal{H}_{y_j}(y^{(pp)}_j) = 1_{y^{(pp)}_j \ge y_j}`. 
Like AIC and BIC, the model with the smallest value of CRPS is considered to be of higher quality than other models. 

.. [Akaike1974] `Akaike, H., A new look at the statistical model identification, Annals of Statistics (2009) <https://dx.doi.org/10.1109/TAC.1974.1100705>`_

.. [Andrieu2009] `Andrieu C., Roberts G.O., The pseudo-marginal approach for efficient Monte Carlo computations, Annals of Statistics (2009) <https://dx.doi.org/10.1214/07-AOS574>`_

.. [Gneiting2007] `Gneiting T., Raftery A., Strictly Proper Scoring Rules, Prediction, and Estimation (2007) <https://sites.stat.washington.edu/raftery/Research/PDF/Gneiting2007jasa.pdf>`_

.. [Haario2001] `Haario H., Saksman E., Tamminen J., An adaptive Metropolis algorithm, Bernoulli (2001) <https://projecteuclid.org/euclid.bj/1080222083>`_

.. [Kass1998] `Raftery A.E., Lewis S., Markov Chain Monte Carlo in Practice: A Roundtable Discussion, The American Statistician (1998) <http://dx.doi.org/10.1080/00031305.1998.10480547>`_

.. [Lloyd2007] `Lloyd-Smith J.O., Maximum Likelihood Estimation of the Negative Binomial Dispersion Parameter for Highly Overdispersed Data, with Applications to Infectious Diseases, Public Library of Science (2007) <http://dx.doi.org/10.1371/journal.pone.0000180>`_

.. [Lynch2004] `Lynch S.M., Western B., Bayesian posterior predictive checks for complex models, Sociological Methods and Research (2004) <http://dx.doi.org/10.1177/0049124103257303>`_

.. [Raftery1992] `Raftery A.E., Lewis S., How Many Iterations in the Gibbs Sampler?, Bayesian Statistics 4 (1992) <http://people.ee.duke.edu/~lcarin/raftery92how.pdf>`_

.. [Schwarz1978] `Schwarz G., Estimating the Dimension of a Model, Bayesian Statistics 4 (1992) <http://dx.doi.org/10.1214/aos/1176344136>`_

