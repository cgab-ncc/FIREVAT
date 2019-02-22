FIREVATStartupMessage <- function()  {
    msg <- c(paste0(
    "
______ _____ _____  _________      __  ________
|  ____|_   _|  __ \\|  ____\\ \\    / /\\|__   __|
| |__    | | | |__) | |__   \\ \\  / /  \\  | |
|  __|   | | |  _  /|  __|   \\ \\/ / /\\ \\ | |
| |     _| |_| | \\ \\| |____   \\  / ____ \\| |
|_|    |_____|_|  \\_\\______|   \\/_/    \\_\\_|
    "))
    return(msg)
}

.onAttach <- function(lib, pkg) {
    packageStartupMessage("Initializing FIREVAT")
    packageStartupMessage(FIREVATStartupMessage())
    packageStartupMessage(paste0("version ", packageVersion("FIREVAT")))
    packageStartupMessage("\"Need a light? You got it.\"")
}
