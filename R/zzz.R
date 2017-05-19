.onAttach <- function(libname, pkgname) {
    citation <- paste0("If you use ", pkgname,
                       " in published research, please cite:")

    packageStartupMessage(citation)
    invisible()
}
