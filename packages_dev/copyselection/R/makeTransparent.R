#' makes a colour transparent
#'
#' This function makes colors transparent.
#' @param someColor the colour that you want to make transparent
#' @param alpha the amount of transparency. Defaults to 140
#' @keywords colors
#' @export
#' @examples
#' makeTransparent("red", 100)


makeTransparent<-function(someColor, alpha=140)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}