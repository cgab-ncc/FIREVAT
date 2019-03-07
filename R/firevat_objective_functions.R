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
    # Calculates the Euclidean distance from the ideal point
    ideal.point <- c(1,1,1,1)
    curr.point <- c(C.refined, (1 - A.refined), C.artifactual, A.artifactual)
    obj.val <- dist(rbind(ideal.point, curr.point))
    return(obj.val)
}


#' @title Exp.Weighted.Obj.Fn
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
Exp.Weighted.Obj.Fn <- function(C.refined, A.refined, C.artifactual, A.artifactual) {
    # Logarithmic weight
    # y = log10(0.9x + 0.1) + 1
    #
    # Exponential weight
    # y = (10^(x-1) - 0.1) / 0.9

    C.refined <- log(((0.9* C.refined) + 0.1), base = 10) + 1
    A.refined <- (10**((1 - A.refined) - 1) - 0.1) / 0.9

    C.artifactual <- log(((0.9* C.artifactual) + 0.1), base = 10) + 1
    A.artifactual <- (10**(A.artifactual - 1) - 0.1) / 0.9

    obj.val <- C.refined * (1 - A.refined) * C.artifactual * A.artifactual
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
    # Calculates the Euclidean distance from the ideal point

    C.refined <- log(((0.9* C.refined) + 0.1), base = 10) + 1
    A.refined <- (10**((1 - A.refined) - 1) - 0.1) / 0.9

    C.artifactual <- log(((0.9* C.artifactual) + 0.1), base = 10) + 1
    A.artifactual <- (10**(A.artifactual - 1) - 0.1) / 0.9

    ideal.point <- c(1,1,1,1)
    curr.point <- c(C.refined, (1 - A.refined), C.artifactual, A.artifactual)
    obj.val <- dist(rbind(ideal.point, curr.point))
    return(obj.val)
}



