# bunsen: Robust Estimation with Covariate Adjustment for Survival Endpoint in Clinical Trials

The bunsen package aims to provide a robust and easy-to-use interface for estimating marginal or unconditional Hazard Ratio (HR) and Restricted Mean Survival Time (RMST) when adjusting prognostic covariates in clinical trials.

## Methodology

bunsen is developed based on three key papers:

-   HR: [Rhian Daniel et al. (2020)](https://onlinelibrary.wiley.com/doi/10.1002/bimj.201900297)
-   RMST: [Theodore Karrison et al. (2018)](https://journals.sagepub.com/doi/10.1177/1740774518759281)
-   Extension of above two methods: [Jiawei Wei et al. (2024)](https://www.tandfonline.com/doi/full/10.1080/19466315.2023.2292774)

## Example

``` r
data('oak')

cox_event <- coxph(Surv(OS, os.status) ~ trt+btmb+pdl1, data=oak)

cox_censor <- coxph(Surv(OS, 1-os.status) ~trt+btmb+pdl1, data=oak)

get_marginal_effect(trt = 'trt',cox_event,cox_censor,M=10000,data=oak)

# Calculating point estimate in local clustermq using multiprocess...
# Starting 4 processes ...
# Running 4 calculations (12 objs/71.2 Kb common; 1 calls/chunk) ...
# Master: [12.2 secs 1.6% CPU]; Worker: [avg 18.6% CPU, max 305.1 Mb]                                                                                                          
# Calculating SE in clustermq using bootstrap N = 1000...
# Submitting 100 worker jobs (ID: cmq6349) ...
# Running 1,000 calculations (23 objs/255.6 Kb common; 1 calls/chunk) ...
# Master: [19.8 secs 7.8% CPU]; Worker: [avg 36.7% CPU, max 308 Mb]                                                                                                            
#         HR         se       2.5%      97.5% 
# -0.4631677  0.1047327 -0.6667394 -0.2504570 
```

## Package authors

-   Xinlei Deng (Maintainer)
-   Dominic Magirr
-   Craig Wang
-   Alexander Przybylski
-   Mark Baillie

## Acknowledgements

-   Jiawei Wei
-   Lukas Andreas Widmer

## References
