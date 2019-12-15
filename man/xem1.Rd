\name{xem1}
\alias{xem1}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("xem1")}
\format{
  The format is:
Formal class 'MixmodCluster' [package "Rmixmod"] with 19 slots
  ..@ strategy   :Formal class 'Strategy' [package "Rmixmod"] with 11 slots
  .. .. ..@ algo             : chr "EM"
  .. .. ..@ nbTry            : num 1
  .. .. ..@ initMethod       : chr "smallEM"
  .. .. ..@ nbTryInInit      : num 10
  .. .. ..@ nbIterationInInit: num 5
  .. .. ..@ nbIterationInAlgo: num 200
  .. .. ..@ epsilonInInit    : num 0.001
  .. .. ..@ epsilonInAlgo    : num 0.001
  .. .. ..@ seed             : num -1
  .. .. ..@ parameter        :Formal class 'GaussianParameter' [package "Rmixmod"] with 4 slots
  .. .. .. .. ..@ mean       : logi[0 , 0 ] 
  .. .. .. .. ..@ variance   :List of 1
  .. .. .. .. .. ..$ : num 0
  .. .. .. .. ..@ proportions: num(0) 
  .. .. .. .. ..@ nbFreeParam: num(0) 
  .. .. ..@ labels           : num(0) 
  ..@ bestResult :Formal class 'MixmodResults' [package "Rmixmod"] with 9 slots
  .. .. ..@ nbCluster     : int 3
  .. .. ..@ model         : chr "Gaussian_pk_Lk_C"
  .. .. ..@ criterion     : chr "BIC"
  .. .. ..@ criterionValue: num 2322
  .. .. ..@ parameters    :Formal class 'GaussianParameter' [package "Rmixmod"] with 4 slots
  .. .. .. .. ..@ mean       : num [1:3, 1:2] 3.91 2.04 4.51 78.38 54.48 ...
  .. .. .. .. ..@ variance   :List of 3
  .. .. .. .. .. ..$ : num [1:2, 1:2] 0.105 0.689 0.689 45.672
  .. .. .. .. .. ..$ : num [1:2, 1:2] 0.0742 0.4847 0.4847 32.1356
  .. .. .. .. .. ..$ : num [1:2, 1:2] 0.0677 0.442 0.442 29.3043
  .. .. .. .. ..@ proportions: num [1:3] 0.24 0.356 0.404
  .. .. .. .. ..@ nbFreeParam: num 13
  .. .. ..@ likelihood    : num -1125
  .. .. ..@ partition     : num [1:272] 1 2 1 2 3 2 3 1 2 3 ...
  .. .. ..@ proba         : num [1:272, 1:3] 9.93e-01 1.20e-10 1.00 3.05e-06 7.05e-02 ...
  .. .. ..@ error         : chr "No error"
  ..@ data       : num [1:272, 1:2] 3.6 1.8 3.33 2.28 4.53 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:2] "Duration" "Waiting.Time"
  ..@ dataType   : chr "quantitative"
  ..@ factor     : num(0) 
  ..@ nbCluster  : num 3
  ..@ knownLabels: int(0) 
  ..@ weight     : num(0) 
  ..@ nbVariable : int 2
  ..@ nbSample   : int 272
  ..@ criterion  : chr "BIC"
  ..@ models     :Formal class 'GaussianModel' [package "Rmixmod"] with 4 slots
  .. .. ..@ family           : chr "general"
  .. .. ..@ listModels       : chr "Gaussian_pk_Lk_C"
  .. .. ..@ free.proportions : logi TRUE
  .. .. ..@ equal.proportions: logi FALSE
  ..@ error      : logi FALSE
  ..@ results    :List of 1
  .. ..$ :Formal class 'MixmodResults' [package "Rmixmod"] with 9 slots
  .. .. .. ..@ nbCluster     : int 3
  .. .. .. ..@ model         : chr "Gaussian_pk_Lk_C"
  .. .. .. ..@ criterion     : chr "BIC"
  .. .. .. ..@ criterionValue: num 2322
  .. .. .. ..@ parameters    :Formal class 'GaussianParameter' [package "Rmixmod"] with 4 slots
  .. .. .. .. .. ..@ mean       : num [1:3, 1:2] 3.91 2.04 4.51 78.38 54.48 ...
  .. .. .. .. .. ..@ variance   :List of 3
  .. .. .. .. .. .. ..$ : num [1:2, 1:2] 0.105 0.689 0.689 45.672
  .. .. .. .. .. .. ..$ : num [1:2, 1:2] 0.0742 0.4847 0.4847 32.1356
  .. .. .. .. .. .. ..$ : num [1:2, 1:2] 0.0677 0.442 0.442 29.3043
  .. .. .. .. .. ..@ proportions: num [1:3] 0.24 0.356 0.404
  .. .. .. .. .. ..@ nbFreeParam: num 13
  .. .. .. ..@ likelihood    : num -1125
  .. .. .. ..@ partition     : num [1:272] 1 2 1 2 3 2 3 1 2 3 ...
  .. .. .. ..@ proba         : num [1:272, 1:3] 9.93e-01 1.20e-10 1.00 3.05e-06 7.05e-02 ...
  .. .. .. ..@ error         : chr "No error"
  ..@ xmlIn      :Formal class 'MixmodXmlInput' [package "Rmixmod"] with 6 slots
  .. .. ..@ file            : chr(0) 
  .. .. ..@ numFormat       : chr "humanReadable"
  .. .. ..@ conversionOnly  : logi FALSE
  .. .. ..@ gaussianModel   :Formal class 'S4' [package ""] with 0 slots
 list()
  .. .. ..@ multinomialModel:Formal class 'S4' [package ""] with 0 slots
 list()
  .. .. ..@ compositeModel  :Formal class 'S4' [package ""] with 0 slots
 list()
  ..@ xmlOut     : chr ""
  ..@ seed       : num -1
  ..@ trace      : num 0
  ..@ massiccc   : num 0
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
data(xem1)
## maybe str(xem1) ; plot(xem1) ...
}
\keyword{datasets}
