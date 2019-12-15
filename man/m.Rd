\name{m}
\alias{m}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("m")}
\format{
  The format is:
List of 9
 $ cluster     : int [1:150] 2 2 2 2 2 2 2 2 2 2 ...
 $ centers     : num [1:2, 1:4] 6.3 5.01 2.89 3.37 4.96 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:2] "1" "2"
  .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
 $ totss       : num 681
 $ withinss    : num [1:2] 123.8 28.6
 $ tot.withinss: num 152
 $ betweenss   : num 529
 $ size        : int [1:2] 97 53
 $ iter        : int 1
 $ ifault      : int 0
 - attr(*, "class")= chr "kmeans"
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
data(m)
## maybe str(m) ; plot(m) ...
}
\keyword{datasets}
