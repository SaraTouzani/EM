\name{EM}
\alias{EM}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("EM")}
\format{
  The format is:
Formal class 'Strategy' [package "Rmixmod"] with 11 slots
  ..@ algo             : chr "EM"
  ..@ nbTry            : num 1
  ..@ initMethod       : chr "smallEM"
  ..@ nbTryInInit      : num 10
  ..@ nbIterationInInit: num 5
  ..@ nbIterationInAlgo: num 200
  ..@ epsilonInInit    : num 0.001
  ..@ epsilonInAlgo    : num 0.001
  ..@ seed             : num -1
  ..@ parameter        :Formal class 'GaussianParameter' [package "Rmixmod"] with 4 slots
  .. .. ..@ mean       : logi[0 , 0 ] 
  .. .. ..@ variance   :List of 1
  .. .. .. ..$ : num 0
  .. .. ..@ proportions: num(0) 
  .. .. ..@ nbFreeParam: num(0) 
  ..@ labels           : num(0) 
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
data(EM)
## maybe str(EM) ; plot(EM) ...
}
\keyword{datasets}
