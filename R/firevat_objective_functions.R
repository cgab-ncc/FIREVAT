# FIREVAT Objective Functions
#
# Last revised date:
#   March 7, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title Default.Obj.Fn
#' @description Calculates the default objective value for FIREVAT GA optimization.
#'
#' @param C.refined A numeric value between 0 and 1.
#' @param A.refined A numeric value between 0 and 1.
#' @param C.artifactual A numeric value between 0 and 1.
#' @param A.artifactual A numeric value between 0 and 1.
#'
#' @return A numeric value between 0 and 1.
#'
#' @export
Default.Obj.Fn <- function(C.refined, A.refined, C.artifactual, A.artifactual) {
    obj.val <- C.refined * (1 - A.refined) * C.artifactual * A.artifactual
    return(obj.val)
}


#' @title Euc.Obj.Fn
#' @description Calculates the Euclidean-distance based objective value for FIREVAT GA optimization.
#'
#' @param C.refined A numeric value between 0 and 1.
#' @param A.refined A numeric value between 0 and 1.
#' @param C.artifactual A numeric value between 0 and 1.
#' @param A.artifactual A numeric value between 0 and 1.
#'
#' @return A numeric value between 0 and 1.
#'
#' @export
Euc.Obj.Fn <- function(C.refined, A.refined, C.artifactual, A.artifactual) {
    # Calculates the Euclidean distance from the reference point
    reference.point <- c(0,0,0,0)
    curr.point <- c(C.refined, (1 - A.refined), C.artifactual, A.artifactual)
    obj.val <- as.numeric(dist(rbind(reference.point, curr.point)))
    return(obj.val)
}


#' @title Exp.Weighted.Obj.Fn.1
#' @description Calculates the exponentially weighted objective value
#' for FIREVAT GA optimization.
#'
#' @param C.refined A numeric value between 0 and 1.
#' @param A.refined A numeric value between 0 and 1.
#' @param C.artifactual A numeric value between 0 and 1.
#' @param A.artifactual A numeric value between 0 and 1.
#'
#' @return A numeric value between 0 and 1.
#'
#' @export
Exp.Weighted.Obj.Fn.1 <- function(C.refined, A.refined, C.artifactual, A.artifactual) {
    # Logarithmic weight
    # y = log10(0.9x + 0.1) + 1
    #
    # Exponential weight
    # y = (10^(x-1) - 0.1) / 0.9

    C.refined <- log(((0.9* C.refined) + 0.1), base = 10) + 1
    A.refined <- (10**((1 - A.refined) - 1) - 0.1) / 0.9

    C.artifactual <- log(((0.9* C.artifactual) + 0.1), base = 10) + 1
    A.artifactual <- (10**(A.artifactual - 1) - 0.1) / 0.9

    obj.val <- C.refined * A.refined * C.artifactual * A.artifactual
    return(obj.val)
}


#' @title Exp.Weighted.Obj.Fn.2
#' @description Calculates the exponentially weighted objective value
#' for FIREVAT GA optimization.
#'
#' @param C.refined A numeric value between 0 and 1.
#' @param A.refined A numeric value between 0 and 1.
#' @param C.artifactual A numeric value between 0 and 1.
#' @param A.artifactual A numeric value between 0 and 1.
#'
#' @return A numeric value between 0 and 1.
#'
#' @export
Exp.Weighted.Obj.Fn.2 <- function(C.refined, A.refined, C.artifactual, A.artifactual) {
    # Logarithmic weight
    # y = 0.5 * log10(0.99x + 0.01) + 1
    #
    # Exponential weight
    # y = (10^(2*(x-1)) - 0.01) / 0.99

    C.refined <- (0.5 * log(((0.99 * C.refined) + 0.01), base = 10)) + 1
    A.refined <- (10**(2*((1 - A.refined) - 1)) - 0.01) / 0.99

    C.artifactual <- (0.5 * log(((0.99 * C.artifactual) + 0.01), base = 10)) + 1
    A.artifactual <- (10**(2*(A.artifactual - 1)) - 0.01) / 0.99

    obj.val <- C.refined * A.refined * C.artifactual * A.artifactual
    return(obj.val)
}


#' @title Euc.Exp.Weighted.Obj.Fn
#' @description
#' Calculates the Euclidean-distance of logarithmically weighted
#' objective value for FIREVAT GA optimization.
#'
#' @param C.refined A numeric value between 0 and 1.
#' @param A.refined A numeric value between 0 and 1.
#' @param C.artifactual A numeric value between 0 and 1.
#' @param A.artifactual A numeric value between 0 and 1.
#'
#' @return A numeric value between 0 and 1.
#'
#' @export
Euc.Exp.Weighted.Obj.Fn <- function(C.refined, A.refined, C.artifactual, A.artifactual) {
    # Logarithmic weight
    # y = log10(0.9x + 0.1) + 1
    #
    # Exponential weight
    # y = (10^(x-1) - 0.1) / 0.9
    #
    # Calculates the Euclidean distance from the reference point

    C.refined <- log(((0.9* C.refined) + 0.1), base = 10) + 1
    A.refined <- (10**((1 - A.refined) - 1) - 0.1) / 0.9

    C.artifactual <- log(((0.9* C.artifactual) + 0.1), base = 10) + 1
    A.artifactual <- (10**(A.artifactual - 1) - 0.1) / 0.9

    reference.point <- c(0,0,0,0)
    curr.point <- c(C.refined, A.refined, C.artifactual, A.artifactual)
    obj.val <- as.numeric(dist(rbind(reference.point, curr.point)))
    return(obj.val)
}


#' @title Euc.Exp.Weighted.Seq.Art.Only.Obj.Fn.1
#' @description
#' Calculates the Euclidean-distance of logarithmically weighted
#' objective value for FIREVAT GA optimization.
#'
#' @param C.refined A numeric value between 0 and 1.
#' @param A.refined A numeric value between 0 and 1.
#' @param C.artifactual A numeric value between 0 and 1.
#' @param A.artifactual A numeric value between 0 and 1.
#'
#' @return A numeric value between 0 and 1.
#'
#' @export
Euc.Exp.Weighted.Seq.Art.Only.Obj.Fn.1 <- function(C.refined, A.refined, C.artifactual, A.artifactual) {
    # Logarithmic weight
    # y = log10(0.9x + 0.1) + 1
    #
    # Exponential weight
    # y = (10^(x-1) - 0.1) / 0.9
    #
    # Calculates the Euclidean distance from the reference point

    A.refined <- (10**((1 - A.refined) - 1) - 0.1) / 0.9
    A.artifactual <- (10**(A.artifactual - 1) - 0.1) / 0.9

    reference.point <- c(0, 0)
    curr.point <- c(A.refined, A.artifactual)
    obj.val <- as.numeric(dist(rbind(reference.point, curr.point)))
    return(obj.val)
}


#' @title Euc.Exp.Weighted.Seq.Art.Only.Obj.Fn.2
#' @description
#' Calculates the Euclidean-distance of logarithmically weighted
#' objective value for FIREVAT GA optimization.
#'
#' @param C.refined A numeric value between 0 and 1.
#' @param A.refined A numeric value between 0 and 1.
#' @param C.artifactual A numeric value between 0 and 1.
#' @param A.artifactual A numeric value between 0 and 1.
#'
#' @return A numeric value between 0 and 1.
#'
#' @export
Euc.Exp.Weighted.Seq.Art.Only.Obj.Fn.2 <- function(C.refined, A.refined, C.artifactual, A.artifactual) {
    # Logarithmic weight
    # y = log10(0.9x + 0.1) + 1
    #
    # Exponential weight
    # y = (10^(x-1) - 0.1) / 0.9
    #
    # Calculates the Euclidean distance from the reference point

    C.artifactual <- log(((0.9 * C.artifactual) + 0.1), base = 10) + 1
    A.artifactual <- (10**(A.artifactual - 1) - 0.1) / 0.9

    reference.point <- c(0,0)
    curr.point <- c(C.artifactual, A.artifactual)
    obj.val <- as.numeric(dist(rbind(reference.point, curr.point)))
    return(obj.val)
}


#' @title Exp.Weighted.Refined.Seq.Art.Only.Obj.Fn
#' @description
#' Calculates the Euclidean-distance of logarithmically weighted
#' objective value for FIREVAT GA optimization.
#'
#' @param C.refined A numeric value between 0 and 1.
#' @param A.refined A numeric value between 0 and 1.
#' @param C.artifactual A numeric value between 0 and 1.
#' @param A.artifactual A numeric value between 0 and 1.
#'
#' @return A numeric value between 0 and 1.
#'
#' @export
Exp.Weighted.Refined.Seq.Art.Only.Obj.Fn <- function(C.refined, A.refined, C.artifactual, A.artifactual) {
    # Logarithmic weight
    # y = log10(0.9x + 0.1) + 1
    #
    # Exponential weight
    # y = (10^(x-1) - 0.1) / 0.9
    #
    # Calculates the Euclidean distance from the reference point

    A.refined <- (10**((1 - A.refined) - 1) - 0.1) / 0.9
    obj.val <- A.refined
    return(obj.val)
}


