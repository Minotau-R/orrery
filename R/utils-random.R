#' samples from the probability density function of a feature in count data.
#'
#' @description See rbeta.
#' @param draws An integer. How many samples should be taken
#' @param k An integer. How many counts were observerd
#' @param n An integer. How many obserations were made
#' @importFrom stats rbeta rnorm
#'
#' @export
#'
sampleBetaBinom <- function(draws, k, n){
    #wrap around rbeta for the special case of Binom()
    rbeta(n = draws, shape1 = k + 1, shape2 = n - k + 1)
}

#' samples from the log() of the geometric mean.
#'
#' @param samples An integer. How many samples should be taken
#' @param count_sample A vector of count data.
#' @param log_transformed A boolean, whether to return log-transformed values.
#'
#' @export
#'
sampleGeomMeam <- function(samples, count_sample, log_transformed = T){
    # Make a matrix with resamples for each feature.
    sample_mat = sampleEachBetaBinom(
        samples = samples, count_sample = count_sample
        )
    # Take the average of the log-transformed values, which corresponds to
    # the log-transformed geometric mean.
    geom_mean = apply(X = log(sample_mat), MARGIN = 1, FUN = mean)

    if(!log_transformed){geom_mean = exp(geom_mean)}

    return(geom_mean)
}

#' samples from the probability density function of a feature in count data for each feature.
#'
#' @param samples An integer. How many samples should be taken
#' @param count_sample A vector of count data.
#'
#' @export
#'
sampleEachBetaBinom <- function(samples, count_sample){
    # FIXME refactor to vapply
    return(
        sapply(count_sample, FUN = function(x){
        sampleBetaBinom(draws = samples, k = x, n = sum(count_sample))
            }, simplify = TRUE)
        )
}

#' samples from the probability density function of a CLR-transformed count sample.
#'
#' @param samples An integer. How many samples should be taken
#' @param count_sample A vector of count data.
#'
#' @export
#'
sampleCLR <- function(samples, count_sample){
    log_ratios = log(
        sampleEachBetaBinom(samples = samples, count_sample = count_sample)
        ) - sampleGeomMeam(
            samples = samples, count_sample = count_sample,
            log_transformed = TRUE
            )

    return(log_ratios)
}

#' samples from the approximated probability density function of the
#'     log-transformed geometric mean of a count sample.
#'
#' @param samples An integer. How many samples should be taken
#' @param count_sample A vector of count data.
#' @param log_transformed A boolean, whether to return log-transformed values.
#'
#' @export
#'
sampleGeomMeanApprox <- function(
        samples, count_sample, log_transformed = TRUE
        ){
    # Estimate a mean and standard deviation for a log-transformed normal
    # approximation.
    mu_hat = mean(
        .getBetaMeans(count_sample = count_sample, log_transformed = TRUE)
        )
    sd_hat = sqrt(
        sum(.getBetaVars(count_sample = count_sample, log_transformed = T)) /
            (length(count_sample)* length(count_sample))
        )

    geom_mean = rnorm(n = samples, mean = mu_hat, sd = sd_hat)

    if(!log_transformed){geom_mean = exp(geom_mean)}

    return(geom_mean)
}

#' samples from the approximated probability density function of a feature in
#'     count data for each feature.
#'
#' @param samples An integer. How many samples should be taken
#' @param count_sample A vector of count data.
#' @param log_transformed A boolean, whether to return log-transformed values.
#'
#' @export
#'
sampleEachApprox <- function(samples, count_sample, log_transformed = TRUE){
    mapply(
        FUN = rnorm,
           mean = .getBetaMeans(count_sample = count_sample,
                               log_transformed = log_transformed),
           sd = sqrt(.getBetaVars(count_sample = count_sample,
                                  log_transformed = log_transformed)),
           MoreArgs = list(n = samples)
           )
}

#' samples from the approximated probability density function of a CLR-transformed count sample.
#'
#' @param samples An integer. How many samples should be taken
#' @param count_sample A vector of count data.
#'
#' @export
#'
sampleCLRApprox <- function(samples, count_sample){

    mu_hat = .getBetaMeans(count_sample = count_sample, log_transformed = TRUE) -
        mean(.getBetaMeans(count_sample = count_sample, log_transformed = TRUE))

    sd_hat = sqrt(
        (.getBetaVars(count_sample = count_sample, log_transformed = TRUE)) +
            (sum(.getBetaVars(
                count_sample = count_sample, log_transformed = TRUE)
                ) /
                 (length(count_sample)* length(count_sample)))
    )

    log_ratios = mapply(FUN = rnorm,
                        mean = mu_hat,
                        sd   = sd_hat,
                        MoreArgs = list(n = samples))
    return(log_ratios)
}



#' samples from the approximated probability density function of a
#'     CLR-transformed count sample.
#'
#' @param samples An integer. How many samples should be taken
#' @param count_table A vector of count data.
#'
#' @export
#'
sampleTableCLR <- function(samples, count_table){
    CLR_out <- matrix(0, nrow = nrow(count_table), ncol = ncol(count_table))

    for(i in 1:samples){
        CLR_out = CLR_out + apply(
            X = count_table, MARGIN = 2, FUN = sampleCLRApprox, samples = 1L
        )
    }
    return(CLR_out/samples)

}

#' samples from the approximated probability density function of a longitudinal
#'     set of CLR-transformed count samples.
#'
#' @param samples An integer. How many samples should be taken.
#' @param count_table A vector of count data.
#' @param cols_as_features A boolean. Toggles whether rows or columns are
#'     samples.
#' @param adjust A boolean. Whether to Apply shrinkage. Experimental.
#'
#' @returns a matrix with colums as samples and rows as features.
#'
#' @export
#'
sampleLongitudinal <- function(
        samples, count_table, cols_as_features = F, adjust = FALSE
        ){
    margin = 1 + cols_as_features
    #nsamples = dim(count_table)[2-cols_as_features]

    mus <- .getTableMeans(
        count_table = count_table, cols_as_features = cols_as_features
        )

    if(adjust){
        mus <- .sCLR(
            count_table = count_table, cols_as_features = cols_as_features
            )
    }

    sds <- apply(X = mus, MARGIN = margin, FUN = sd)

    nrows <- length(sds)

    matrix(
        rnorm(samples * nrows, rep(mus, samples), rep(sds, samples)),
        nrow = nrows
        )

}


#' samples from the probability density function of the difference of two
#'     CLR-transformed count samples.
#'
#' @param samples An integer. How many samples should be taken
#' @param x A vector of count data.
#' @param y A second vector of count data.
#' @export
#'
sampleCLRdiff <- function(samples, x, y){
    sampleCLR(samples = samples, count_sample = x) -
        sampleCLR(samples = samples, count_sample = y)
}

#' samples from the probability density function of the aitchison distance of
#'     two count samples.
#'
#' @param samples An integer. How many samples should be taken
#' @param x A vector of count data.
#' @param y A second vector of count data.
#' @export
#'
sampleAitchison <- function(samples, x, y){
    sqrt(rowSums(
        (sampleCLRdiff(samples, x, y)^2)
    ))
}
