
#' @title PrintLog
#' @description
#' Prints log message
#'
#' @param msg String value message to print along with log type and date
#'
#' @export
PrintLog <- function(msg, type = "INFO") {
    print(paste0(type, " [", Sys.time(), "] ", msg))
}
