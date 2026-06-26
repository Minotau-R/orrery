#' Estimate the approximated means from a CLR-transformed count sample.
#'
#' @description See rbeta.
#' @param count_table A table of count data with rows as features and columns as
#'     samples.
#' @param CLR_transformed A boolean, whether to return centered log-ratio (CLR)
#'     transformed values.
#' @param cols_as_features A boolean. Toggles whether rows or columns are
#'     samples
#'
.getTableMeans <- function(
        count_table, CLR_transformed = TRUE, cols_as_features = FALSE
        ){
    margin = 2 - cols_as_features

    if(CLR_transformed){
        table_means = apply(X = count_table,
                            MARGIN = margin,
                            FUN = .getCLRMeans)
    }
    if(!CLR_transformed){
        table_means = apply(X = count_table,
                            MARGIN = margin,
                            FUN = .getBetaMeans,
                            log_transformed = FALSE)
    }

    if(cols_as_features){table_means = t(table_means)}

    return(table_means)
}

#' Estimate the approximated vars from a CLR-transformed count sample.
#'
#' @description See rbeta.
#' @param count_table A table of count data with rows as features and columns as
#'     samples.
#' @param CLR_transformed A boolean, whether to return centered log-ratio (CLR)
#'     values.
#' @param cols_as_features A boolean. Toggles whether rows or columns are
#'     samples.
#'
.getTableVars <- function(
        count_table, CLR_transformed = TRUE, cols_as_features = FALSE
        ){
    margin = 2 - cols_as_features

    if(CLR_transformed){
        table_vars = apply(X = count_table,
                           MARGIN = margin,
                           FUN = .getCLRVars)
    }
    if(!CLR_transformed){
        table_vars = apply(X = count_table,
                           MARGIN = margin,
                           FUN = .getBetaVars,
                           log_transformed = FALSE)
    }

    if(cols_as_features){table_vars = t(table_vars)}

    return(table_vars)
}


#' compute CLR using Aitchison's method
#' @description CLR transform a table.
#' @param x A table of compositional data without zeroes.
#' @param cols_as_features A boolean. Toggles whether rows or columns are samples.
#' @return A table of CLR-transformed data
#'
.clr <- function(x, cols_as_features = FALSE){
    margin = 2 - cols_as_features
    apply(X = x, MARGIN = margin, FUN = function(y){log(y)-mean(log(y))})
}

#' Compute shrunk proportions using the 'qsh' method described by Ionas Erb
#' \href{https://doi.org/10.48550/arXiv.2205.09215}{here}.
#'
#' @description Shrink relative count data using the arithmetic mean.
#' @param x A table of compositional data.
#' @param cols_as_features A boolean. Toggles whether rows or columns are samples.
#' @return A relative count table shrunk towards the arithmetic mean.
#'
.estQshr <- function(x, cols_as_features = FALSE){
    margin = 2 - cols_as_features

    apply(X = x, MARGIN = margin, FUN = function(n){
        #Total counts
        N <- sum(n)

        #Relative counts
        q <- n/N

        #Number of components
        D <- length(n)

        #Define maximum entropy (equidistribution) prior
        Tau <- rep(1/D,D)

        #Define the shrinkage parameter lambda
        l <- (1-sum(q^2)) / ((N-1)*sum((Tau-q)^2))

        #Bound lamda at 0 and 1.
        if(l<0){l=0}
        if(l>1){l=1}

        # Shrink the empirical estimator of the relative counts data
        # by mixing it with the naive equidistant prior.
        return(l/D+(1-l)*n/N)
    })
}

#' Shrink and CLR-transform a relative count table.
#'
#' @description Compute shrunk proportions using the 'qsh' method described by Ionas Erb \href{https://doi.org/10.48550/arXiv.2205.09215}{here}.
#' @param count_table A table of count data with rows as features and columns as samples.
#' @param cols_as_features A boolean. Toggles whether rows or columns are samples.
#' @export
#'
.sCLR <- function(count_table, cols_as_features = FALSE){

    #Shrink the count table
    res <- .estQshr(count_table, cols_as_features = cols_as_features)

    #CLR transform the shrunk relative count estimates. Leave boolean at FALSE.
    res = .clr(res, cols_as_features = FALSE)

    if(cols_as_features){res = t(res)}

    return(res)
}

#' Estimate volatility from count data in Aitchision distance.
#'
#' @description Compute feature-wise variance over time for a group of
#'     measurements.
#' @param count_table A table of count data with rows as features and columns as
#'     samples.
#' @param cols_as_features A boolean. Toggles whether rows or columns are
#'     samples.
#' @param adjust A boolean. Whether to Apply shrinkage. Experimental.
#'
#'
#'
.estVolatility <- function(
        count_table, cols_as_features = FALSE, adjust = FALSE
        ){

    vars <- .estVariance(count_table = count_table,
                        cols_as_features = cols_as_features,
                        adjust = adjust)

    sum(vars)
}

#' Estimate by feature variance from count data in Aitchision distance.
#'
#' @description Compute feature-wise variance over time for a group of measurements.
#' @param count_table A table of count data with rows as features and columns as samples.
#' @param cols_as_features A boolean. Toggles whether rows or columns are samples.
#' @param adjust A boolean. Whether to Apply shrinkage. Experimental.
#'
#' @importFrom stats sd var
#' @export
#'
.estVariance <- function(count_table, cols_as_features = FALSE, adjust = TRUE){
    margin = 1 + cols_as_features
    #nsamples = dim(count_table)[2-cols_as_features]

    if(!adjust){
        estMeans <- .getTableMeans(
            count_table = count_table, cols_as_features = cols_as_features
            )
    }

    if(adjust){
        estMeans <- .sCLR(
            count_table = count_table, cols_as_features = cols_as_features
            )
    }

    apply(X      = estMeans,
          MARGIN = margin,
          FUN    = var)
}
