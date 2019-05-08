# FIREVAT Log Functions
#
# Last revised date:
#   May 08, 2019
#
# Authors:
#   Andy Jinseok Lee (jinseok.lee@ncc.re.kr)
#   Hyunbin Kim (khb7840@ncc.re.kr)
#   Bioinformatics Analysis Team, National Cancer Center Korea


#' @title PrintLog
#' @description
#' Prints log message
#'
#' @param msg String value message to print along with log type and date
#' @param type String value that represents type of this message. 'INFO' by default.
#'
#' @export
PrintLog <- function(msg, type = "INFO") {
    print(paste0(type, " [", Sys.time(), "] ", msg))
}
