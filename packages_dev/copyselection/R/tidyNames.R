#' this is an internal function that converts a list of names
#'
#' This function converts lists of names for various plotting functions used in various plotting functions.
#' @param x a list of names
#' @param fula should FULAI and FULAII (and MANDINKAI) be kept as is? Defaults to FALSE
#' @param pel Should PELII be kept as is? Defaults to FALSE
#' @keywords Busby_bespoke
#' @export
#' @examples
#' tidyNames("ESOMALI", fula=T)

tidyNames <- function(x,fula=F,pel=F,khoesan=F,tig=F){
  x <- gsub("IBD","",x)
  x <- gsub("CULTIVATOR","",x)
  x <- gsub("BLACKSMITH","",x)
  x <- gsub("ESOMALI","SOMALI",x)
  x <- gsub("JUHOANSI","JUHOAN",x)
  x <- gsub("JUHOAN","JU/HOANSI",x)
  x <- gsub("XUNV","XUN",x)
  x <- gsub("SEMI.BANTU","SEMI-BANTU",x) 
  x <- gsub("SWBANTU","HERERO",x)
  x <- gsub("KHOMANI","=KHOMANI",x)
  x <- gsub("GUIGHANAKGAL","/GUI//GHANA_KGAL",x)
  x <- gsub("LWK","LUHYA",x)
  x <- gsub("MKK","MAASAI",x)
  x <- gsub("YRI","YORUBA",x)
  if(pel==F)
  {
    x <- gsub("PELII","PEL",x)
  }
  if(fula==F)
  {
    x <- gsub("FULAII","FULA",x)
    x <- gsub("FULAI","FULA",x)
    x <- gsub("MANDINKAII","MANDINKA",x)
    x <- gsub("MANDINKAI","MANDINKA",x)
  }
  
  if(khoesan == T)
  {
      x <- gsub("/GUI//GHANA_KGAL","/GUI//GANA",x)
      x <- gsub("JU/HOANSI","JU/\'HOANSI",x)
      x <- gsub("XUN","!XUN",x)
      x <- gsub("==KHOMANI","=KHOMANI",x)
  }
  if(tig == T)
  {
     x <- gsub("TYGRAY","TIGRAY",x)
  }
  
  return(x)
}

# tidyNamesUnique <- function(x){
#   x <- unique(gsub("IBD","",x))
#   x <- unique(gsub("CULTIVATOR","",x))
#   x <- unique(gsub("BLACKSMITH","",x))
#   x <- unique(gsub("ESOMALI","SOMALI",x))
#   x <- unique(gsub("JUHOANSI","JUHOAN",x))
#   x <- unique(gsub("JUHOAN","JUHOANSI",x))
#   x <- unique(gsub("XUNV","XUN",x))
#   x <- unique(gsub("SWBANTU","HERERO",x))
#   x <- unique(gsub("KHOMANI","=KHOMANI",x))
#   x <- unique(gsub("GUIGHANAKGAL","/GUI//GHANA_KGAL",x))
#   x <- unique(gsub("JUHOANSI","JU/HOANSI",x))
#   x <- unique(gsub("LWK","LUHYA",x))
#   x <- unique(gsub("MKK","MAASAI",x))
#   x <- unique(gsub("YRI","YORUBA",x))
#   
#   return(x)
# }