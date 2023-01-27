###############################################################################
###############################################################################
###############################################################################

getMySurvivalTable <- function(
    
    my_times,
    my_events
    
){
    
    # '''
    # Returns a survival table, i. e. for a given time point, it outputs
    # -- time points ("time_point") 
    # -- numbers of individuals at risk ("n_at_risk"),
    # -- numbers of events ("n_of_events"),
    # -- numbers of censored individuals ("n_of_censored")
    # using two vectors of the same length -- "my_times" and "my_events"
    # for time and event component of time-to-event variable.
    # '''
    
    my_output <- NULL
    
    for(
        my_time in sort(unique(my_times))
    ){
        
        my_output <- rbind(
            my_output,
            c(
                "time_point" = my_time,
                "n_at_risk" = length(
                    my_events[
                        which(
                            my_times >= my_time
                        )
                    ]
                ),
                "n_of_events" = sum(
                    my_events[
                        which(
                            my_times == my_time
                        )
                    ]
                )
            )
        )
        
    }
    
    my_output <- my_output[
        -which(
            my_output[, "n_of_events"] == 0
        )
        ,
    ]
    
    n_of_censored <- NULL
    
    for(
        i in 1:(dim(my_output)[1] - 1)
    ){
        
        n_of_censored <- c(
            
            n_of_censored,
            my_output[i, "n_at_risk"] - 
            my_output[i, "n_of_events"] - 
            my_output[i + 1, "n_at_risk"]
            
        )
        
    }
    
    my_output <- data.frame(
        
        my_output,
        "n_of_censored" = c(
            n_of_censored,
            NA
        ),
        stringsAsFactors = FALSE
        
    )
    
    return(
        my_output
    )
    
}


## ----------------------------------------------------------------------------

###############################################################################
###############################################################################
###############################################################################





