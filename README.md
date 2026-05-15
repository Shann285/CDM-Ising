# CDM-Ising
Bayesian regularized estimation of cognitive diagnosis models incorporating Ising networks

We propose two generalizations of CDMs, where the dependence between latent attributes or item response residuals is modeled by an Ising network. The first framework, termed as latent Ising CDM (LI-CDM), formulates an Ising network between latent attributes. This framework facilitates the exploratory estimation of the independence or dependence among latent attributes. The second framework, termed as residual Ising CDM (RI-CDM), models the residuals of item responses using an Ising network. Within this framework, it is feasible to develop CDMs that depart from the assumption of local item independence. Bayesian regularized procedures are developed for estimating both models, allowing parameter estimation and network structure identification to be conducted within a unified framework.

The LICDM-N500K3L.R file is the Bayesian regularized procedure for the LI-CDMs under the condition of N=500, K=3 and large interactions.

The LICDM-N500K4L.R file is the Bayesian regularized procedure for the LI-CDMs under the condition of N=500, K=4 and large interactions.

The RICDM-N500K3L.R file is the Bayesian regularized procedure for the RI-CDMs under the condition of N=500, K=3 and large interactions.

The RICDM-N500K4L.R file is the Bayesian regularized procedure for the RI-CDMs under the condition of N=500, K=4 and large interactions.

The 2026-data-fraction2.R file includes the estimation procedures for the traditional DINA model and the RI-CDM for the real data.
