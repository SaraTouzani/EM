\name{iris.res}
\alias{iris.res}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("iris.res")}
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
  .. .. ..@ criterionValue: num 605
  .. .. ..@ parameters    :Formal class 'GaussianParameter' [package "Rmixmod"] with 4 slots
  .. .. .. .. ..@ mean       : num [1:3, 1:4] 5.93 6.55 5.01 2.76 2.97 ...
  .. .. .. .. ..@ variance   :List of 3
  .. .. .. .. .. ..$ : num [1:4, 1:4] 0.2152 0.0889 0.1233 0.0331 0.0889 ...
  .. .. .. .. .. ..$ : num [1:4, 1:4] 0.4058 0.1675 0.2324 0.0624 0.1675 ...
  .. .. .. .. .. ..$ : num [1:4, 1:4] 0.155 0.064 0.0888 0.0238 0.064 ...
  .. .. .. .. ..@ proportions: num [1:3] 0.313 0.354 0.333
  .. .. .. .. ..@ nbFreeParam: num 26
  .. .. ..@ likelihood    : num -238
  .. .. ..@ partition     : num [1:150] 3 3 3 3 3 3 3 3 3 3 ...
  .. .. ..@ proba         : num [1:150, 1:3] 1.10e-24 1.22e-20 1.99e-22 3.52e-19 5.83e-25 ...
  .. .. ..@ error         : chr "No error"
  ..@ data       : num [1:150, 1:4] 5.1 4.9 4.7 4.6 5 5.4 4.6 5 4.4 4.9 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..@ dataType   : chr "quantitative"
  ..@ factor     : num(0) 
  ..@ nbCluster  : num 3
  ..@ knownLabels: int(0) 
  ..@ weight     : num(0) 
  ..@ nbVariable : int 4
  ..@ nbSample   : int 150
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
  .. .. .. ..@ criterionValue: num 605
  .. .. .. ..@ parameters    :Formal class 'GaussianParameter' [package "Rmixmod"] with 4 slots
  .. .. .. .. .. ..@ mean       : num [1:3, 1:4] 5.93 6.55 5.01 2.76 2.97 ...
  .. .. .. .. .. ..@ variance   :List of 3
  .. .. .. .. .. .. ..$ : num [1:4, 1:4] 0.2152 0.0889 0.1233 0.0331 0.0889 ...
  .. .. .. .. .. .. ..$ : num [1:4, 1:4] 0.4058 0.1675 0.2324 0.0624 0.1675 ...
  .. .. .. .. .. .. ..$ : num [1:4, 1:4] 0.155 0.064 0.0888 0.0238 0.064 ...
  .. .. .. .. .. ..@ proportions: num [1:3] 0.313 0.354 0.333
  .. .. .. .. .. ..@ nbFreeParam: num 26
  .. .. .. ..@ likelihood    : num -238
  .. .. .. ..@ partition     : num [1:150] 3 3 3 3 3 3 3 3 3 3 ...
  .. .. .. ..@ proba         : num [1:150, 1:3] 1.10e-24 1.22e-20 1.99e-22 3.52e-19 5.83e-25 ...
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
data(iris.res)
## maybe str(iris.res) ; plot(iris.res) ...
}
\keyword{datasets}
