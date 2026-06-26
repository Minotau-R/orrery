#' Calculate by feature variance from a beta distributions
#'
#' @description See rbeta.
#' @param count_sample A vector of count data.
#' @param log_transformed A boolean, whether to return log-transformed values.
#'
#' @export
#'
.getBetaVars <- function(count_sample, log_transformed = F){

    if(!log_transformed){
        return(sapply(count_sample,
                      FUN = function(x){.getBetaVar(k = x, n = sum(count_sample))}, simplify = T))
    }

    if(log_transformed){
        return(sapply(count_sample,
                      FUN = function(x){getLogBetaVar(k = x, n = sum(count_sample))}, simplify = T))
    }

}


#' Calculate the variance from a beta distribution.
#'
#' @description See rbeta. Called by .getBetaVars.
#' @param k An integer. How many counts were observerd
#' @param n An integer. How many obserations were made
#'
#' @export
#'
.getBetaVar <- function(k, n){
    a = k + 1
    b = n - k + 1
    return(a*b/((a+b*b)*(a + b + 1)))
}


#' Calculate the log variance from a beta distribution.
#'
#' @description See rbeta. Called by .getBetaVars.
#' @param k An integer. How many counts were observerd
#' @param n An integer. How many obserations were made
#'
#' @export
#'
getLogBetaVar <- function(k, n){
    a = k + 1
    b = n - k + 1
    return(trigamma(a) - trigamma(a + b))
}


#' Calculate by feature mean from a beta distributions
#'
#' @description See rbeta.
#' @param count_sample A vector of count data.
#' @param log_transformed A boolean, whether to return log-transformed values.
#'
#' @export
#'
.getBetaMeans <- function(count_sample, log_transformed = F){

    if(!log_transformed){
        return(sapply(count_sample, FUN = function(x){.getBetaMean(k = x, n = sum(count_sample))}, simplify = T))
    }
    if(log_transformed){
        return(sapply(count_sample, FUN = function(x){getLogBetaMean(k = x, n = sum(count_sample))}, simplify = T))
    }

}


#' Calculate the mean from a beta distribution.
#'
#' @description See rbeta. Called by .getBetaMeans.
#' @param k An integer. How many counts were observed
#' @param n An integer. How many observations were made
#'
#' @export
#'
.getBetaMean <- function(k, n){
    a = k + 1
    b = n - k + 1
    return(a/(a+b))
}



#' Calculate the log mean from a beta distribution.
#'
#' @description See rbeta. Called by .getBetaMeans.
#' @param k An integer. How many counts were observerd
#' @param n An integer. How many obserations were made
#'
#' @export
#'
getLogBetaMean <- function(k, n){
    a = k + 1
    b = n - k + 1
    return(digamma(a) - digamma(a + b))
}


#' Calculate by feature mean from a CLR transformed approximation
#'
#' @description See rbeta.
#' @param count_sample A vector of count data.
#'
#' @export
#'
.getCLRMeans <- function(count_sample){
    .getBetaMeans(count_sample = count_sample, log_transformed = T) -
        mean(.getBetaMeans(count_sample = count_sample, log_transformed = T))
}


#' Calculate by feature variance from a CLR transformed approximation
#'
#' @description See rbeta.
#' @param count_sample A vector of count data.
#'
#' @export
#'
.getCLRVars <- function(count_sample){
    (.getBetaVars(count_sample = count_sample, log_transformed = T)) +
        (
            sum(.getBetaVars(count_sample = count_sample, log_transformed = T)) /
                (length(count_sample)* length(count_sample))
        )

}

#Leave these here for now
#' #' Calculate the covariance of two marginals in a Dirichlet.
#' #'
#' #' @description See rbeta.
#' #' @param k_i An integer. How many counts were observed in the first marginal
#' #' @param k_j An integer. How many counts were observed in the first marginal
#' #' @param n An integer. How many observations were made
#' #'
#' #' @export
#' #'
#' .getCovDirichlet <- function(k_i, k_j, n){
#'   a_i = k_i + 1
#'   a_j = k_j + 1
#'   a_0 = n + 1
#'   -(a_i)*(a_j)/(((a_0)*(a_0)) *((a_0) + 1))
#' }
#'
#'
#' #' Calculate the sum of all covariances between marginals in a Dirichlet.
#' #'
#' #' @description See rbeta.
#' #' @param count_sample A vector of count data.
#' #'
#' #' @export
#' #'
#' .getCovSum <- function(count_sample){
#'   n = sum(count_sample)
#'   covInd <- combn(x = count_sample, m = 2)
#'   covSum = 0
#'   for(i in 1:ncol(covInd)){
#'     covSum = covSum + .getCovDir(k_i = covInd[1,i], k_j = covInd[2,i], n = n)
#'   }
#'   return(covSum*2)
#' }
