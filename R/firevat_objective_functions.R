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

    obj.val <- sqrt(C.refined**2 + (1 - A.refined)**2 + C.artifactual**2 + A.artifactual**2)
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

    obj.val <- sqrt(C.refined**2 + A.refined**2 + C.artifactual**2 + A.artifactual**2)
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

    obj.val <- sqrt(A.refined**2 + A.artifactual**2)
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

    obj.val <- sqrt(C.artifactual**2 + A.artifactual**2)
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


#' @title Test.Obj.Fn.1
#' @description
#' Test objective function 1
#'
#' @param C.refined A numeric value between 0 and 1.
#' @param A.refined A numeric value between 0 and 1.
#' @param C.artifactual A numeric value between 0 and 1.
#' @param A.artifactual A numeric value between 0 and 1.
#'
#' @return A numeric value between 0 and 1.
#'
#' @export
Test.Obj.Fn.1 <- function(C.refined, A.refined, C.artifactual, A.artifactual) {
    # Logarithmic weight
    # y = log10(0.9x + 0.1) + 1
    #
    # Exponential weight
    # y = (10^(x-1) - 0.1) / 0.9
    #
    # Calculates the Euclidean distance from the reference point

    C.refined <- log(((0.9* C.refined) + 0.1), base = 10) + 1
    A.refined <- 1 - A.refined
    obj.val <- C.refined * A.refined * C.artifactual
    return(obj.val)
}


#' @title Test.Obj.Fn.2
#' @description
#' Test objective function 2
#'
#' @param C.refined A numeric value between 0 and 1.
#' @param A.refined A numeric value between 0 and 1.
#' @param C.artifactual A numeric value between 0 and 1.
#' @param A.artifactual A numeric value between 0 and 1.
#'
#' @return A numeric value between 0 and 1.
#'
#' @export
Test.Obj.Fn.2 <- function(C.refined, A.refined, C.artifactual, A.artifactual) {
    # Logarithmic weight
    # y = log10(0.9x + 0.1) + 1
    #
    # Exponential weight
    # y = (10^(x-1) - 0.1) / 0.9
    #
    # Calculates the Euclidean distance from the reference point


    C.refined <- log(((0.9* C.refined) + 0.1), base = 10) + 1
    A.refined <- (10**((1 - A.refined) - 1) - 0.1) / 0.9
    C.artifactual <- (10**(C.artifactual - 1) - 0.1) / 0.9

    obj.val <- C.refined * A.refined * C.artifactual
    return(obj.val)
}


#' @title Sigmoid.Obj.Fn
#' @description
#' Sigmoid objective function
#'
#' @param C.refined A numeric value between 0 and 1.
#' @param A.refined A numeric value between 0 and 1.
#' @param C.artifactual A numeric value between 0 and 1.
#' @param A.artifactual A numeric value between 0 and 1.
#'
#' @return A numeric value between 0 and 1.
#'
#' @export
Sigmoid.Obj.Fn <- function(C.refined, A.refined, C.artifactual, A.artifactual) {
    # y = 1 / (1 + e^(-k(b+ax^d)))
    A.refined <- (1 - A.refined)

    # Relaxed
    a <- 10
    b <- -6
    d <- 5
    k <- 3
    C.refined <- (1 / (1 + exp(-1* k * (b + a * (C.refined ** d)))))
    C.artifactual <- (1 / (1 + exp(-1* k * (b + a * (C.artifactual ** d)))))

    # Stringent
    # a <- 10
    # b <- -8.5
    # d <- 5
    # k <- 10
    A.refined <- (1 / (1 + exp(-1 * k * (b + a * (A.refined ** d)))))
    A.artifactual <- (1 / (1 + exp(-1 * k * (b + a * (A.artifactual ** d)))))

    obj.val <- C.refined * A.refined * C.artifactual * A.artifactual
    return(obj.val)
}


#' @title Leaky.ReLU.Obj.Fn
#' @description
#' Lkeay ReLU objective function
#'
#' @param C.refined A numeric value between 0 and 1.
#' @param A.refined A numeric value between 0 and 1.
#' @param C.artifactual A numeric value between 0 and 1.
#' @param A.artifactual A numeric value between 0 and 1.
#'
#' @return A numeric value between 0 and 1.
#'
#' @export
Leaky.ReLU.Obj.Fn <- function(C.refined, A.refined, C.artifactual, A.artifactual) {
    # if less than 0.9:
    #   y = 0.05x
    # else:
    #   y = 10(x - 0.9)
    A.refined <- (1 - A.refined)

    if (C.refined < 0.9) {
        C.refined <- 0.111111 * C.refined
    } else {
        C.refined <- 9 * (C.refined + (-0.888888888))
    }

    if (A.refined < 0.9) {
        A.refined <- 0.111111 * A.refined
    } else {
        A.refined <- 9 * (A.refined + (-0.888888888))
    }

    if (C.artifactual < 0.9) {
        C.artifactual <- 0.111111 * C.artifactual
    } else {
        C.artifactual <- 9 * (C.artifactual + (-0.888888888))
    }

    if (A.artifactual < 0.9) {
        A.artifactual <- 0.111111 * A.artifactual
    } else {
        A.artifactual <- 9 * (A.artifactual + (-0.888888888))
    }

    obj.val <- C.refined * A.refined * C.artifactual * A.artifactual
    return(obj.val)
}


#' @title Leaky.ReLU.A.Ref.Obj.Fn
#' @description
#' Leaky ReLU objective function
#'
#' @param C.refined A numeric value between 0 and 1.
#' @param A.refined A numeric value between 0 and 1.
#' @param C.artifactual A numeric value between 0 and 1.
#' @param A.artifactual A numeric value between 0 and 1.
#'
#' @return A numeric value between 0 and 1.
#'
#' @export
Leaky.ReLU.A.Ref.Obj.Fn <- function(C.refined, A.refined, C.artifactual, A.artifactual) {
    # if greater than 0.5:
    #   y = 0.01(1-x)
    # else:
    #   y = 2((1-x) - 0.5) + 0.005
    if (A.refined > 0.5) {
        A.refined <- 0.01 * (1 - A.refined)
    } else {
        A.refined <- 2 * ((1 - A.refined) - 0.5) + 0.005
    }

    obj.val <- C.refined * A.refined * C.artifactual * A.artifactual
    return(obj.val)
}


#' @title Leaky.ReLU.A.Art.Obj.Fn
#' @description
#' Leaky ReLU objective function
#'
#' @param C.refined A numeric value between 0 and 1.
#' @param A.refined A numeric value between 0 and 1.
#' @param C.artifactual A numeric value between 0 and 1.
#' @param A.artifactual A numeric value between 0 and 1.
#'
#' @return A numeric value between 0 and 1.
#'
#' @export
Leaky.ReLU.A.Art.Obj.Fn <- function(C.refined, A.refined, C.artifactual, A.artifactual) {
    # if less than 0.5:
    #   y = 0.01(1-x)
    # else:
    #   y = 2((1-x) - 0.5) + 0.005
    if (A.artifactual < 0.5) {
        A.artifactual <- 0.01 * A.artifactual
    } else {
        A.artifactual <- 2 * (A.artifactual - 0.5) + 0.005
    }

    A.refined <- (1 - A.refined)

    obj.val <- C.refined * A.refined * C.artifactual * A.artifactual
    return(obj.val)
}


#' @title Exp.Weighted.A.Ref.Obj.Fn
#' @description
#' Exponentially weighted objective function
#'
#' @param C.refined A numeric value between 0 and 1.
#' @param A.refined A numeric value between 0 and 1.
#' @param C.artifactual A numeric value between 0 and 1.
#' @param A.artifactual A numeric value between 0 and 1.
#'
#' @return A numeric value between 0 and 1.
#'
#' @export
Exp.Weighted.A.Ref.Obj.Fn <- function(C.refined, A.refined, C.artifactual, A.artifactual) {
    # Exponentially weight A.refined
    # y = 10^e((1-x) - 1)
    A.refined <- 10 ** (exp(1)*((1 - A.refined) - 1))

    obj.val <- C.refined * A.refined * C.artifactual * A.artifactual
    return(obj.val)
}


#' @title Exp.Weighted.A.Art.Obj.Fn
#' @description
#' Exponentially weighted objective function
#'
#' @param C.refined A numeric value between 0 and 1.
#' @param A.refined A numeric value between 0 and 1.
#' @param C.artifactual A numeric value between 0 and 1.
#' @param A.artifactual A numeric value between 0 and 1.
#'
#' @return A numeric value between 0 and 1.
#'
#' @export
Exp.Weighted.A.Art.Obj.Fn <- function(C.refined, A.refined, C.artifactual, A.artifactual) {
    # Exponentially weight A.artifactual
    # y = 10^e((1-x) - 1)
    A.artifactual <- 10 ** (exp(1)*(A.artifactual - 1))
    A.refined <- (1 - A.refined)

    obj.val <- C.refined * A.refined * C.artifactual * A.artifactual
    return(obj.val)
}


#' @title Sig.Extraction.Obj.Fn
#' @description Calculates the signature extraction objective value for FIREVAT GA optimization.
#'
#' @param C.refined A numeric value between 0 and 1.
#' @param A.refined A numeric value between 0 and 1.
#' @param C.artifactual A numeric value between 0 and 1.
#' @param A.artifactual A numeric value between 0 and 1.
#'
#' @return A numeric value between 0 and 1.
#'
#' @export
Sig.Extraction.Obj.Fn <- function(C.refined, A.refined, C.artifactual, A.artifactual) {
    obj.val <- C.refined * A.refined * C.artifactual * (1 - A.artifactual)
    return(obj.val)
}
