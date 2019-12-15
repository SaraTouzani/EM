\name{cl}
\alias{cl}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("cl")}
\format{
  The format is:
List of 6
 $ cluster   : int [1:150] 1 1 1 1 1 1 1 1 1 1 ...
 $ size      : 'table' int [1:2(1d)] 88 62
  ..- attr(*, "dimnames")=List of 1
  .. ..$ cluster: chr [1:2] "1" "2"
 $ modes     :'data.frame':	2 obs. of  2 variables:
  ..$ Species         : Factor w/ 3 levels "setosa","versicolor",..: 1 3
  ..$ Sepal.Length_cat: chr [1:2] "petit" "moyen"
 $ withindiff: num [1:2] 69 33
 $ iterations: int 2
 $ weighted  : logi FALSE
 - attr(*, "class")= chr "kmodes"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(cl)
## maybe str(cl) ; plot(cl) ...
}
\keyword{datasets}
