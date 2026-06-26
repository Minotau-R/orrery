#' Compute the geometric mean of a vector.
#'
#' @description Simply shorthand for exp(mean(log(x))).
#' @param x a numeric vector that does not contain zeroes.
#'
#' @export
#'
g_mean <- function(x){exp(mean(log(x)))}

#' Center by the geometric mean.
#'
#' @description Simply divide each value by its corresponding geometric mean.
#' @param x a numeric table that does not contain zeroes.
#'
#' @export
#'
g_center_by <- function(x, MARGIN = MARGIN){
    return(x / apply(X = x, MARGIN = MARGIN, FUN = mean))
}

#' Undo CLR transformation with softmax
#' @description softmax function to transform vector to relative abundance/proportions.
#' @param x A clr-transformed vector
#' @return A relative abundance/proportions vector.
#'
softmax <- function(x) {exp(x) / sum(exp(x))}


#' Center by margin
#' @param x a table
#' @param MARGIN see apply
#'
g_center_by <- function(x, MARGIN = MARGIN){
    return(x / apply(X = x, MARGIN = MARGIN, FUN = mean))
}


#' Double center a table.
#'
#' @description To center: Divide each value by its corresponding row and column geometric mean.
#' @param x a numeric table.
#' @param type One of c("auto", "numeric", "integer"). What type of numbers are represented in the table.
#' @param m Defaults to 0.65. The scaling factor to impute zeroes in the non-integer case.
#' @param log.transform A boolean. Whether to return log-transformed values.
#'
#' @export
#'
double_center <- function(x, type = "auto", m = 0.65, log.transform = FALSE){
    x <- impute_and_scale(x, type = type, m = m)
    x <- log(x)
    x <- sweep(x, MARGIN = 1, rowMeans(x), `-`, check.margin = TRUE)
    x <- sweep(x, MARGIN = 2, colMeans(x), `-`, check.margin = TRUE)

    if(!log.transform) {return(exp(x))}

    return(x)
}

#' Perform a double-centered log-ratio transformation
#' @rdname double_center
#' @export
#'
dclr <- function(x, type = "auto", m = 0.65){
    double_center(x = x, type = type, m = 0.65, log.transform = TRUE)
}


#' Perform a double-centered log-ratio transformation
#' @rdname double_center
#' @export
#'
dclr_c <- function(x, m = 0.65){
    double_center(x = x, type = "numeric", m = 0.65, log.transform = TRUE)
}

#' Perform a double-centered log-ratio transformation
#' @rdname double_center
#' @export
#'
dclr_p <- function(x){
    double_center(x = x, type = "int", log.transform = TRUE)
}

#' Impute zeroes and scale to sum to 1.
#'
#' @description Scales both count data and generic compositional data to proportions
#' @param x a numeric table.
#' @param type One of c("auto", "numeric", "integer"). What type of numbers are represented
#' @param m Defaults to 0.65. The scaling factor to impute zeroes in the non-integer case.
#'
#'
impute_and_scale <- function(x, type = "auto", m = 0.65){
    stopifnot(type %in% c("auto", "numeric", "integer"))

    rn <- row.names(x)

    if(type == "auto"){
        ifelse(all(is.wholenumber(x)),
               yes = {x = apply(x, 2, as.integer); `<-`(type, "integer")},
               no  = {`<-`(type, "numeric")}
        )
    }

    if(type == "integer"){x = impute_and_scale.int(x)}
    if(type == "numeric"){x = impute_and_scale.num(x, m = m)}

    x = sweep(x = x, MARGIN = 2, STATS = colSums(x), FUN = "/", check.margin = FALSE)

    return(`row.names<-`(x, rn))
}

#' Check whether input is (clooe to) counts
#' @noRd
#'
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

#' Replace zeroes with NA
#' @noRd
#'
zero_to_NA <- function(x) {x[x==0] <- NA;x}

#' Find column-wise non-zero minimum
#' @noRd
#'
find_colmins <- function(x){
    apply(x, MARGIN = 2, FUN = min, na.rm = TRUE)
}

#'@noRd
#'
impute_and_scale.int <- function(x){
    .getTableMeans(x, CLR_transformed = FALSE)
}

#'@noRd
#'
impute_and_scale.num <- function(x, m){
    x <- zero_to_NA(x)
    cM <- find_colmins(x) * m
    indx <- which(is.na(x), arr.ind = TRUE)
    x[indx] <- cM[indx[,2]]
    x
}
