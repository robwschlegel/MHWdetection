#' A placeholder function so pkgdown can compile.
#'
#' This isn't sriously meant to be used for anything, It is just here to trick
#' pkgdown into functioning more like blgdown.
#'
#' @param x A single integer or numeric value.
#'
#' @return The function will return double whatever integer is put into it.
#'
#' @author Robert W. Schlegel
#'
#' @export
#'
#' @examples
#' res <- placeholder(24)
placeholder <- function(x){
  x2 = x *2
  return(x2)
}
