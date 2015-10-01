#' convert factors in a dataframe to strings
#'
#' This function converts all the columns containing factors to strings
#' @param dataframe The dataframe containing factors that you want to convert
#' @keywords factors
#' @export
#' @examples
#' convert.factors.to.strings.in.dataframe(dataframe)


convert.factors.to.strings.in.dataframe <- function(dataframe)
{
  class.data  <- sapply(dataframe, class)
  factor.vars <- class.data[class.data == "factor"]
  for (colname in names(factor.vars))
  {
    dataframe[,colname] <- as.character(dataframe[,colname])
  }
  return (dataframe)
}