#' @title Write netrc file
#' @description Write a netrc file that is valid for accessing urs.earthdata.nasa.gov
#' @details
#' The database is accessed with the user's credentials.
#' A netrc file storing login and password information is required.
#' See \href{https://urs.earthdata.nasa.gov/}{here}
#' for instruction on how to register and set DataSpace credential.
#' @param login A character. Email address used for logging in on DataSpace.
#' @param password A character. Password associated with the login.
#' @param netrcFile A character. A path to where the netrc file should be written.
#' By default will go to your home directory, which is advised
#' @param overwrite A logical. overwrite the existing netrc file?
#'
#' @return A character vector containing the netrc file path
#' @seealso \code{\link{connectDS}} \code{\link{checkNetrc}}
#' @examples
#' # First, create an account in the DataSpace App and read the terms of use
#' # Next, create a netrc file using writeNetrc()
#' writeNetrc(
#'   login = "XXX@email.com",
#'   password = "yourSecretPassword"
#' )
#' # Specify `netrcFile = getNetrcPath()` to write netrc in the default path
#' @export
writeNetrc <- function(login,
                       password,
                       machine = 'urs.earthdata.nasa.gov',
                       netrcFile =  getNetrcPath(),
                       overwrite = FALSE) {
  
  if (checkNetrc() && !overwrite) {
    stop("'", netrcFile, "' already exists. Set `overwrite=TRUE`
         if you'd like to overwrite.",
         call. = FALSE
    )
  }
  
  string <- paste(
    "\nmachine ", machine,
    "login", login,
    "password", password
  )
  
  # create a netrc file
  write(string,  path.expand(netrcFile), append=TRUE)
  
  # set the owner-only permission
  Sys.chmod(netrcFile, mode = "600")
  
  netrcFile
}

#' @title Check netrc file
#' @description Check that there is a netrc file with a valid
#' entry for urs.earthdata.nasa.gov.
#' @param netrcFile A character. File path to netrc file to check.
#' @param onStaging A logical. Whether to check the staging server instead
#' of the production server.
#' @return logical
#' @seealso \code{\link{writeNetrc}}
#' @export
#'
checkNetrc <- function(netrcFile = getNetrcPath()) {
  
  if (!file.exists(netrcFile)) { return(FALSE) }
  
  lines <- gsub("http.*//", "", readLines(netrcFile))
  
  return(any(grepl("urs.earthdata.nasa.gov", lines)))
}


#' @title Get a default netrc file path
#'
#' @description Get a default netrc file path
#'
#' @return A character vector containing the default netrc file path
#'
#' @examples
#' getNetrcPath()
#' @export
getNetrcPath <- function() {
  home <- Sys.getenv("HOME")
  if (whatOS() == "windows") {
    file.path(home, "_netrc")
  } else {
    file.path(home, ".netrc")
  }
}

whatOS = function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}
