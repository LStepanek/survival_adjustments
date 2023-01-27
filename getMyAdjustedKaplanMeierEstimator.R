###############################################################################
###############################################################################
###############################################################################

getMyAdjustedKaplanMeierEstimator <- function(
    
    my_survival_table,
    my_delta = 0.5,
    confidence_bounds = c("greenwood", "log-log"),
    confidence_level = 0.05
    
){
    
    # '''
    # On input, it uses "my_survival_table", which is an output
    # of function getMySurvivalTable(), and a value or vector
    # of delta parameter "my_delta".
    # On output, it returns adjusted Kaplan-Meier estimator and a lower
    # and upper bound of 100 * (1 - "confidence_level") % confidence
    # interval based on eithe Greenwood formula (when confidence_bounds
    # = "greenwood"), or log-log transformation (when confidence_bounds
    # = "log-log").
    
    if(
        length(my_delta) == 1
    ){
        my_delta <- rep(
            my_delta,
            dim(my_survival_table)[1]
        )
    }
    
    my_kaplan_meier <- rep(0, dim(my_survival_table)[1])
    my_summation <- rep(0, dim(my_survival_table)[1])
    my_lower_bound <- rep(0, dim(my_survival_table)[1])
    my_upper_bound <- rep(0, dim(my_survival_table)[1])
    
    my_kaplan_meier[1] <- 1 - my_delta[1] * my_survival_table[
        1,
        "n_of_events"
    ] / my_survival_table[
        1,
        "n_at_risk"
    ]
    
    if(
        dim(my_survival_table)[1] > 1
    ){
        
        for(
            i in 2:dim(my_survival_table)[1]
        ){
            
            my_kaplan_meier[i] <- my_kaplan_meier[i - 1] * (
                1 - my_delta[i] * my_survival_table[
                    i,
                    "n_of_events"
                ] / my_survival_table[
                    i,
                    "n_at_risk"
                ]
            )
            
        }
        
    }
    
    my_summation[1] <- my_delta[1] * my_survival_table[
        1,
        "n_of_events"
    ] / (
        my_survival_table[
            1,
            "n_at_risk"
        ] * (
            my_survival_table[
                1,
                "n_at_risk"
            ] - my_delta[1] * my_survival_table[
                1,
                "n_of_events"
            ]
        )
    )
    
    if(
        dim(my_survival_table)[1] > 1
    ){
        
        for(
            i in 2:dim(my_survival_table)[1]
        ){
            
            my_summation[i] <- my_summation[i - 1] + (
                my_delta[i] * my_survival_table[
                    i,
                    "n_of_events"
                ] / (
                    my_survival_table[
                        i,
                        "n_at_risk"
                    ] * (
                        my_survival_table[
                            i,
                            "n_at_risk"
                        ] - my_delta[i] * my_survival_table[
                            i,
                            "n_of_events"
                        ]
                    )
                )
            )
            
        }
        
    }
    
    if(
        confidence_bounds == "greenwood"
    ){
        
        my_lower_bound <- my_kaplan_meier - qnorm(
            1 - confidence_level / 2
        ) * sqrt(
            (my_kaplan_meier ^ 2) * my_summation
        )
        my_upper_bound <- my_kaplan_meier + qnorm(
            1 - confidence_level / 2
        ) * sqrt(
            (my_kaplan_meier ^ 2) * my_summation
        )
        
    }
    
    if(
        confidence_bounds == "log-log"
    ){
        
        my_lower_bound <- my_kaplan_meier ^ (
            exp(
                + qnorm(
                    1 - confidence_level / 2
                ) * sqrt(
                    1 / (log(my_kaplan_meier) ^ 2) * my_summation
                )
            )
        )
        my_upper_bound <- my_kaplan_meier ^ (
            exp(
                - qnorm(
                    1 - confidence_level / 2
                ) * sqrt(
                    1 / (log(my_kaplan_meier) ^ 2) * my_summation
                )
            )
        )
        
    }
    
    return(
        data.frame(
            
            "kaplan_meier_survival" = my_kaplan_meier,
            "confidence_lower_bound" = my_lower_bound,
            "confidence_upper_bound" = my_upper_bound,
            stringsAsFactors = FALSE
            
        )
    )
    
}


## ----------------------------------------------------------------------------

###############################################################################
###############################################################################
###############################################################################





