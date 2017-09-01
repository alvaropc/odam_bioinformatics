\name{lcount}
\alias{lcount}
\title{
Count an object into a list
}
\description{
Having a list, you can count an object in the whole list or by index.
}
\usage{
lcount(lst, obj,index=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lst}{
A \code{list} in which you want to count an object.
}
 \item{obj}{
 A \code{object} you want to count into your list.
 }
 \item{index}{
 An \code{index} in which you want to count into your list.
 }
}


\value{
The function return the sum of objects found.
}

\author{
Alvaro Ponce Cabrera}


\examples{
l <- list(1,1:3,1:4)
lcount(l,1) #3
lcount(l,1,2) #1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~list }
\keyword{ ~count }
