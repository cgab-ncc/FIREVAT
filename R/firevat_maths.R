# FIREVAT Math Functions
#
# Last revised date:
#   February 19, 2019
#
# Authors: 
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title ComputeZScore
#'
#' @description Returns a z-score of x given a distribution of values
#'
#' @param values a numeric vector
#' @param x a numeric value
#'
#' @return a numeric value corresponding to the z-score of x
#'
#' @export
#' @importFrom stats sd
ComputeZScore <- function(values, x)  {
    population.standard.dev <- sd(values)*sqrt((length(values)-1)/(length(values)))
    population.mean <- mean(values)
    z.score <- (x - population.mean) / population.standard.dev
    return(z.score)
}


#' @title ComputeZScoreEquiValue
#'
#' @description Returns a numeric value that is equivalent to the specified z.score
#'   in the distribution of 'values'
#'
#' @param z.score numeric value
#' @param values numeric vector
#'
#' @return a numeric value corresponding to the specified z.score in the 'values' distribution
#'
#' @export
#' @importFrom stats sd
ComputeZScoreEquiValue <- function(z.score, values)  {    
    population.standard.dev <- sd(values)*sqrt((length(values)-1)/(length(values)))
    population.mean <- mean(values)
    
    z.score.equivalent.val <- z.score * population.standard.dev + population.mean
    return(z.score.equivalent.val)
}


#' @title DecimalCeiling
#'
#' @description Returns the ceiling of a decimal value
#'   e.g. value = 0.15, decimal = 0.1 returns 0.2
#'
#' @param value numeric value (decimal)
#' @param decimal numeric value (e.g. 0.1, 0.001)
#'
#' @return a numeric value
#'
#' @export
DecimalCeiling <- function(value, decimal)  {    
    ceil.value <- ceiling(value / decimal) * decimal
    return(ceil.value)
}