#' SpecDist(x, datatype=¡¨abundance¡¨) is a function of estimating species rank distribution for abundance/incidence based data
#' @param x a vector of species abundance or incidence frequency. If datatype = "incidence", then the input format of first entry should be total number of sampling units, and followed by species incidence frequency.
#' @param datatype the data type (abundance or incidence) of input.
#' @return a data.frame object of RAD/RID
SpecDist <- function(x, datatype="abundance"){
  TYPE <- c("abundance", "incidence")
  if(is.na(pmatch(datatype, TYPE)))
    stop("invalid datatype")
  if(pmatch(datatype, TYPE) == -1)
    stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)
  
  if(datatype=="abundance"){  
    out <- rbind(data.frame("probability"=DetAbu(x, zero=TRUE), 
                            "method"="detected"),
                 data.frame("probability"=UndAbu(x), 
                            "method"="undetected"))
  }else if(datatype=="incidence"){
    out <- rbind(data.frame("probability"=DetInc(x, zero=TRUE), 
                            "method"="detected"),
                 data.frame("probability"=UndInc(x), 
                            "method"="undetected"))
  }
  out[order(-out[,1]),]
}

