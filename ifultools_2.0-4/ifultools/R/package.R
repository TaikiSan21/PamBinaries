".onLoad" <- function(libname, pkgname)
{   
}

".onUnload" <- function(libpath)
{
}

".onAttach" <- function(libname, pkgname)
{
  if (interactive())
  {
    packageStartupMessage('Loading ', pkgname, ' version ', as.character(packageVersion(pkgname)))
  }
}