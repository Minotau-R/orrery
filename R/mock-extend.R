#' Extend a community forward in time with sample-aware naming
#'
#' @description Predict forward simulation of a compositional time series
#' @param init_composition named vector of compositions from previous sample
#' @param A fitted temporal dynamics matrix
#' @param Q fitted process covariance matrix
#' @param alpha fitted overdispersion parameter
#' @param n_time numeric. Number of future timepoints to simulate
#' @param read_depth vector or scalar. Sequencing depth per timepoint
#' @param n_replicates numeric. Number of independent simulation paths
#' @param sample_name character. Prefix for sample naming
#' @param start_timepoint numeric. Last observed timepoint index
#' @param technical_replicate logical. If TRUE, all replicates will have the same observed composition at each timepoint. If FALSE, each replicate will have its own observed composition.
#' @importFrom compositions ilr acomp ilrInv
#' @importFrom MASS mvrnorm
#' @importFrom dirmult rdirichlet
#'

orrey_extend <- function(
    init_composition,
    A,
    Q,
    alpha = 10,
    n_time = 10,
    read_depth = rep(50000, n_time),
    n_replicates = 1,
    sample_name = "sample",
    start_timepoint = 0,
    technical_replicate = FALSE
){

p0 <- init_composition / sum(init_composition)
K <- length(p0)

if(length(read_depth) == 1){
    read_depth <- rep(read_depth, n_time)
}
n_time <- as.integer(n_time)
n_replicates <- as.integer(n_replicates)
stopifnot("read_depth must be a scalar or a vector of length n_time"= length(read_depth) == n_time)
stopifnot("Q must be a symmetric matrix" = isSymmetric(Q))
stopifnot("n_time must be a positive integer" = (n_time > 0 & n_time == round(n_time)))
stopifnot("n_replicates must be a positive integer" = (n_replicates > 0 & n_replicates == round(n_replicates)))
print(eigen(Q)$values)
stopifnot("All eigenvalues of Q must be positive" = all(eigen(Q)$values > 0))
species_names <- names(p0)
if (is.null(species_names)) {
    species_names <- paste0("sp", seq_len(K))
}

sample_ids <- as.vector(outer(seq_len(max(1, n_replicates)), start_timepoint + seq_len(n_time), FUN = function(r, t) if (n_replicates > 1) paste0(sample_name, "_t", t, "_r", r) else paste0(sample_name, "_t", t)))

n_samples <- length(sample_ids)


Z_store <- matrix(
    NA,
    nrow = K - 1,
    ncol = n_samples,
    dimnames = list(
        paste0("ilr", seq_len(K - 1)),
        sample_ids
    )
)

P_store <- matrix(
    NA,
    nrow = K,
    ncol = n_samples,
    dimnames = list(
        species_names,
        sample_ids
    )
)

Y_store <- matrix(
    NA,
    nrow = K,
    ncol = n_samples,
    dimnames = list(
        species_names,
        sample_ids
    )
)

z_prev <- ilr(acomp(p0))

for(i in seq_len(n_time)){

    mu_t <- as.vector(A %*% z_prev)

    z_t <- mvrnorm(
        n = 1,
        mu = mu_t,
        Sigma = Q
    )

    p_t <- as.vector(ilrInv(z_t))
    p_t <- p_t / sum(p_t)
    names(p_t) <- species_names

    z_prev <- z_t
    theta <- p_t * alpha
    if (technical_replicate) {


    p_obs <- rdirichlet(1, theta)
    p_obs <- as.vector(p_obs)
    names(p_obs) <- species_names

    for(r in seq_len(n_replicates)){

        y_t <- as.vector(rmultinom(
            1,
            size = read_depth[i],
            prob = p_obs
        ))
        names(y_t) <- species_names

        idx <- (i - 1) * n_replicates + r

        Z_store[, idx] <- z_t
        P_store[, idx] <- p_t
        Y_store[, idx] <- y_t    }

} else {


    for(r in seq_len(n_replicates)){

        p_obs <- rdirichlet(1, theta)
        p_obs <- as.vector(p_obs)
        names(p_obs) <- species_names

        y_t <- as.vector(rmultinom(
            1,
            size = read_depth[i],
            prob = p_obs
        ))
        names(y_t) <- species_names

        idx <- (i - 1) * n_replicates + r

        Z_store[, idx] <- z_t
        P_store[, idx] <- p_t
        Y_store[, idx] <- y_t    }
}}



return(list(
    sample_ids = sample_ids,
    true_latent_z = Z_store,
    true_composition = P_store,
    observed_counts = Y_store,
    parameters = list(
        A = A,
        Q = Q,
        alpha = alpha,
        read_depth = read_depth,
        sample_name = sample_name,
        start_timepoint = start_timepoint,
        species_names = species_names,
        n_replicates = n_replicates
    )
))
}

#' Generate mock parameters for testing
#'
#' @description Generate mock parameters for testing the orrey_extend function.
#' @param n_features numeric. Number of features in the dataset.
#' @param persistence numeric. Persistence parameter for the dynamics matrix.
#' @param interaction_sd numeric. Standard deviation of the interaction terms in the dynamics matrix.
#' @param individual_sd numeric. Standard deviation of independent species-specific variability.
#' @param shared_sd numeric. Standard deviation of shared environmental variability.
#' @param n_factors numeric. Number of latent ecological factors influencing covariance.
#' @param factor_sd numeric. Standard deviation of latent factor effects.

#'
#' @return A list containing the mock parameters.
#'
mock_parameters <- function(
    n_features = 10,
    persistence = 0.95,
    interaction_sd = 0.01,
    individual_sd = 0.05,
    shared_sd = 0.02,
    n_factors = 3,
    factor_sd = 0.05
){

    init <- setNames(runif(n_features), paste0("sp", seq_len(n_features)))
    init <- init / sum(init)

    D <- n_features - 1

    A <- diag(persistence, D)
    A <- A + matrix(rnorm(D^2, sd = interaction_sd), D, D)
    diag(A) <- persistence
    Q_individual <- diag(individual_sd^2, D)

    u <- rep(1, D)
    Q_shared <- shared_sd^2 * outer(u, u)
    L <- matrix(rnorm(D * n_factors, sd = factor_sd), D, n_factors)
    Q_factor <- L %*% t(L)

    Q <- Q_individual + Q_shared + Q_factor

    Q <- (Q + t(Q)) / 2

    list(init = init, A = A, Q = Q)
}