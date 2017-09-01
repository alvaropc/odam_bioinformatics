\name{ConvertToList}
\alias{ConvertToList}
\docType{data}
\title{
Example of a list}
\description{
Dataset created from a list which needs to be converted as a list again to test AlvaroPonce package functions. See example for more details.
}
\usage{data("ConvertToList")}
\format{
  A data set which needs to be converted to a list with index equal to 5.
  \describe{
   When converted, \code{List} in which you can find characters and numeric objects in order to test the functions 
  }
}


\examples{
data(ConvertToList)
ExampleList <- apply(ConvertToList,2,list)
FinalList=list()
for (i in (1:length(ExampleList))){
  l=unlist(ExampleList[[i]])
  FinalList=lappend(FinalList,l)
  }

lappend(FinalList,1:2)
lcount(FinalList,1)
lcount(FinalList,1,2)
lremove(FinalList,1)
lremove(FinalList,1,2)
}
