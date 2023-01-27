# `survival_adjustments`

A repository for adjusted Kaplan-Meier survival estimations funded by the grant OP VVV IGA/A, CZ.02.2.69/0.0/0.0/19_073/0016936 with no. 18/2021, which has been provided by the Internal Grant Agency of the Prague University.

## Installation

```
for(
    my_function in c(
        
        "getMySurvivalTable",
        "getMyKaplanMeierEstimator",
        "getMyAdjustedKaplanMeierEstimator"
        
    )
){
    
    source(
        paste(
            "https://raw.githubusercontent.com/LStepanek/",
            "survival_adjustments/main/",
            my_function,
            ".R",
            sep = ""
        )
    )
    
}

```

<p align="center">
  <img src="https://raw.githubusercontent.com/LStepanek/survival_adjustments/main/_EU_grant_funding_logo_.jpg">
</p>
