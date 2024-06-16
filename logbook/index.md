# Simple Difference-in-Differences Estimation in Fixed-T Panels

[Nicholas Brown](https://sites.google.com/msu.edu/nicholasbrown)<sup>1</sup>,
[Kyle Butts](https://www.kylebutts.com/)<sup>2</sup>, and
[Joakim Westerlund](https://www.deakin.edu.au/about-deakin/people/joakim-westerlund)<sup>3,4</sup>
<br>
<sup>1</sup>Queen's University, <sup>2</sup>University of Colorado: Boulder, <sup>3</sup>Lund University, <sup>4</sup>Deakin University

#### [Paper](https://arxiv.org/abs/2301.11358)


## Abstract

The present paper proposes a new treatment effects estimator that is valid when the
number of time periods is small, and the parallel trends condition holds conditional on
covariates and unobserved heterogeneity in the form of interactive fixed effects. The 
estimator also allow the control variables to be affected by treatment and it enables 
estimation of the resulting indirect effect on the outcome variable. The asymptotic 
properties of the estimator are established and their accuracy in small samples is 
investigated using Monte Carlo simulations. The empirical usefulness of the estimator 
is illustrated using as an example the effect of increased trade competition on firm 
markups in China.


## Replication

### Monte Carlo Simulations

1. `code/simulations/simulation-1.jl`

 - Simulations presented in Table 1 and Table 2
 - There are set of helper functions in `code/simulations/helpers.module.jl` containing estimation functions and the DGP.

### Application

1. `code/Trade-Liberalization-and-Markup-Dispersion/analysis.R`

- Takes the data from [Lu and Yu (2015)](https://www.aeaweb.org/articles?id=10.1257/app.20140350) and recreates Figure 2 in their paper.
- Estimates the C<sup>2</sup>ED<sup>2</sup> model to estimate the effect of TWO ascension on markup dispersion in high-tariff industries.
- Decomposes the effect into the mediated effect via decrease in marginal cost dispersion.


## Citation

```
@article{brown2023difference,
  title={Difference-in-Differences via Common Correlated Effects},
  author={Brown, Nicholas and Butts, Kyle and Westerlund, Joakim},
  journal={arXiv preprint arXiv:2301.11358},
  year={2023}
}
```

