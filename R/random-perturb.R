#' Perturb a vector
#'
#' @description Randomly perturb a CLR-transformed vector by a given aitchison distance.
#' @param x A vector of CLR-transformed data.
#' @param by A numeric. How much aitchison distance to perturb the vector x by.
#' @param prop A proportion (0-1). What proportion of the features should be the perturbation be divided over.
#' @importFrom dirmult rdirichlet
#'
perturb_by <- function(x, by, prop = 1){
    #stopifnot("All features in x must be positive." = all(x > 0))
    stopifnot("prop must be between 0-1." = (prop <=1 & prop >=0))

    #Which features to alter
    target = sample(1:length(x), max(2, round(length(x)*prop), digits = 0))

    #Generate container fro perturnbation
    res_perturb = vector(mode = "numeric", length = length(target))

    #Decide which features are going to be positive and negative
    division <- c(T, F, sample(c(T, F), size = length(target)-2, replace = T))

    #Distribute half of the squared distance to the positive division
    res_perturb[ division] =  sqrt(((by^2)/2)*(
        rdirichlet(n = 1, alpha = rep(1, sum(division))))
        )

    #Distribute half of the squared distance to the negative division
    res_perturb[!division] = -sqrt(((by^2)/2)*(
        rdirichlet(n = 1, alpha = rep(1, sum(!division)))))

    #Perturb the input
    x[target] = x[target] + res_perturb

    return(x)
}

#' Perturb a CLR-transformed vector
#'
#' @description Randomly perturb a CLR-transformed vector by a given aitchison
#'     distance.
#' @param x A vector of CLR-transformed data.
#' @param by A numeric. How much aitchison distance to perturb the vector x by.
#' @param prop A proportion (0-1). What proportion of the features should be the
#'     perturbation be divided over.
#'
#' @return A perturbed vector
#'
#' @export
#'
perturb_by_CLR <- function(x, by, prop = 1){
    perturb_by(x = x, by = by, prop = prop)
}

#' Perturb a relative abundance vector
#'
#' @description Randomly perturb a relative abundance vector by a given
#'     aitchison distance.
#' @param x A vector of relative abundance data.
#' @param by A numeric. How much aitchison distance to perturb the vector x by.
#' @param prop A proportion (0-1). What proportion of the features should be the
#'     perturbation be divided over.
#' @param alpha An numeric. The alpha parameter for rdirichlet.
#'
#' @return A perturbed vector
#'
#' @export
#'
perturb_by_relab <- function(x, by, prop = 1, alpha = 1){
    #stopifnot("All features in x must be positive." = all(x > 0))
    stopifnot("prop must be between 0-1." = (prop <=1 & prop >=0))

    x_nonzero = x != 0
    #Which features to alter
    target = sample(
        seq_len(sum(x_nonzero)), max(2, round(sum(x_nonzero)*prop), digits = 0)
        )

    #Generate container for perturbation
    res_perturb = rep(1, length(target))

    #Decide which features are going to be positive and negative
    division = c(T, F, sample(c(T, F), size = length(target)-2, replace = T))

    #Distribute half of the squared distance to the positive division
    res_perturb[ division] =  exp(sqrt(
        ((by^2)/2) *
            rdirichlet(n = 1, alpha = rep(alpha, sum(division)))))


    #Distribute half of the squared distance to the negative division
    res_perturb[!division] = exp(-sqrt(
        ((by^2)/2) *
            rdirichlet(n = 1, alpha = rep(alpha, sum(!division)))
    ))

    #Perturb the input
    x[x_nonzero][target] = x[x_nonzero][target] * res_perturb

    return(x)
}


#' Perturb a count data vector
#'
#' @description Randomly perturb a count data vector by a given aitchison
#'     distance.
#' @param x A vector of count data.
#' @param by A numeric. How much aitchison distance to perturb the vector x by.
#' @param prop A proportion (0-1). What proportion of the features should be
#'     the perturbation be divided over.
#'
#' @return A perturbed vector
#'
#' @export
#'
perturb_by_count <- function(x, by, prop = 1){
    #Record total observations
    s = sum(x)

    #CLR transform
    x = log(x) - mean(log(x))

    #Perturb
    x = perturb_by(x = x, by = by, prop = prop)

    #softmax to un-CLR
    x = exp(x) / sum(exp(x))

    #Return result
    return(x * (s/length(x)))
}

