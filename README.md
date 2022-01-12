# Accounting for temporal change in multiple biodiversity patterns improves the inference of metacommunity processes

This repository contains code and data associated with the manuscript “Accounting for temporal change in multiple biodiversity patterns improves the inference of metacommunity processes” published in *Ecology*.

Guzman, L. M., Thompson, P. L., Viana, D.S., Vanschoenwinkel, B., Horváth, Z., Ptacnik, R., Jeliazkov, A., Gascón, S., Lemmens, P., Anton-Pardo, M., Langenheder, S., Meester, L.D., Chase, J.M

## In a nutshell

Metacommunity ecology has focused on using observational and analytical approaches to disentangle effects of processes on community assembly. Many methods have been proposed for this purpose, most notably multivariate analyses of species abundance and its association with variation in spatial and environmental conditions. Here, we simulate metacommunities emphasizing three main processes that underlie a continuum of metacommunity dynamics following the paper by [Thompson et al. 2020](https://onlinelibrary.wiley.com/doi/10.1111/ele.13568). We then use random forests to evaluate the performance of classical and novel time-integrated summary statistics of metacommunity dynamics to understand the strength of these three processes. We found that: (i) time series are necessary to disentangle metacommunity processes, (ii) each of the three studied processes is distinguished with different descriptors, (iii) each descriptor is differently sensitive to temporal and spatial sampling effort.

The code used to generate the simulations of the metacommunities can be found [here](https://zenodo.org/record/3708350#.X1Fw2y2z0Us).

The code presented here calcualtes the summary statistics for each simulation run and runs random forests to compare the sets of statistics. 

This project was funded by [iDiv](https://www.idiv.de/en/index.html)
