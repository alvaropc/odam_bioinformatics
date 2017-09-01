\name{lremove}
\alias{lremove}
\title{
Remove an object in a list
}
\description{
Having a list, you can remove an object in the whole list or by index.
}
\usage{
lremove(lst, obj,index=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lst}{
A \code{list} in which you want to remove an object.
}
 \item{obj}{
 A \code{object} you want to remove into your list.
 }
 \item{index}{
 An \code{index} in which you want to remove an objet into your list.
 }
}


\value{
The function return the list without the objects removed.
}

\author{
Alvaro Ponce Cabrera}


\examples{
l <- list(1,1:3,1:4)
lremove(l,1) 
lremove(l,1,2) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~list }
\keyword{ ~remove }
