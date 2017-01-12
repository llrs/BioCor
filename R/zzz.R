##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {

    pkgVersion <- packageDescription(pkgname, fields = "Version")
    msg <- paste0(pkgname, " v", pkgVersion, "  ",
                  "For issues: https://github.com/llrs/", pkgname, "\n\n")

    citation <- paste0("If you use ", pkgname,
                       " in published research, please cite:\n")

    packageStartupMessage(msg, citation)
    invisible()
}
