\name{lappend}
\alias{lappend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add an object to a list
}
\description{
Having a list, you can add a new object in the last part of the list as a new index.
}
\usage{
lappend(lst, obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lst}{
A \code{list} in which you want to add a new object. It will be added as a new index in the final part of the list
}
 \item{obj}{
 A \code{object} you want to add to your list, it could be everything you want to add (vectors, data frames, characters...).
 }
 
}


\value{
The function return the same list you enter with one more index in the final part of the function with your object
}

\author{
Alvaro Ponce Cabrera}


\examples{
l <- list(1,1:3,1:4)
l<-lappend(l,1:5)
l

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~list }
