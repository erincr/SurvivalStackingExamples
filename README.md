# Survival Stacking in Action
This repository includes code to accompany the paper: "A review of survival stacking: a method to cast survival regression analysis as a classification problem".

Abstract:
> While there are many well-developed data science methods for classification and regression, there are relatively few methods for working with right-censored data. 
> Here, we review survival stacking, a method for casting a survival regression analysis problem as a classification problem, thereby allowing the use of general 
> classification methods and software in a survival setting. Inspired by the Cox partial likelihood, survival stacking collects features and outcomes of survival 
> data in a large data frame with a binary outcome. We show that survival stacking with logistic regression is approximately equivalent to the Cox proportional hazards 
> model. We further illustrate survival stacking on real and simulated data. By reframing survival regression problems as classification problems, survival stacking removes 
> the reliance on specialized tools for survival regression, and makes it straightforward for data scientists to use well-known learning algorithms and software for 
> classification in the survival setting. This in turn lowers the barrier for flexible survival modeling.

This repository includes three examples of survival stacking in action. All are similar: they perform survival regression, and then measure and plot performance. 
Common functions are in the file `Helper.R`. The other files show examples on different datasets:
1. `Rotterdam_GBSG_Example.R` performs survival regression for the Rotterdam and GBSG datasets[^1].
2. `Simulated_Examples.R` does survival regression for right censored data simulated using `simsurv`[^2].
3. `Simulated_Truncated_Examples.R` does survival regression for simulated data as in (2), additionally with left truncation.

[^1]: Patrick Royston and Douglas G Altman. External validation of a Cox prognostic model: principles and methods. BMC
medical research methodology, 13(1):1–15, 2013.

[^2]: Samuel L. Brilleman, Rory Wolfe, Margarita Moreno-Betancur, and Michael J. Crowther. Simulating survival data using
the simsurv R package. Journal of Statistical Software, 97(3):1–27, 2020.

