\name{xem}
\alias{xem}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("xem")}
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
  .. .. ..@ nbCluster     : int 2
  .. .. ..@ model         : chr "Binary_pk_Ekjh"
  .. .. ..@ criterion     : chr "BIC"
  .. .. ..@ criterionValue: num 519
  .. .. ..@ parameters    :Formal class 'MultinomialParameter' [package "Rmixmod"] with 5 slots
  .. .. .. .. ..@ center     : num [1:2, 1:5] 2 1 2 3 2 1 2 1 1 1
  .. .. .. .. ..@ scatter    :List of 2
  .. .. .. .. .. ..$ : num [1:5, 1:5] 0.428 0.1203 0.0509 0.3641 0.1074 ...
  .. .. .. .. .. ..$ : num [1:5, 1:5] 0.4937 0.0761 0.1521 0.039 0.0577 ...
  .. .. .. .. ..@ factor     : num [1:5] 2 4 5 5 3
  .. .. .. .. ..@ proportions: num [1:2] 0.346 0.654
  .. .. .. .. ..@ nbFreeParam: num 29
  .. .. ..@ likelihood    : num -198
  .. .. ..@ partition     : num [1:69] 1 1 2 2 2 2 2 1 1 1 ...
  .. .. ..@ proba         : num [1:69, 1:2] 0.998001 0.998377 0.000245 0.019384 0.019384 ...
  .. .. ..@ error         : chr "No error"
  ..@ data       : int [1:69, 1:5] 1 2 2 1 1 1 2 2 2 2 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:69] "1" "2" "3" "4" ...
  .. .. ..$ : chr [1:5] "gender" "eyebrow" "collar" "sub-caudal" ...
  ..@ dataType   : chr "qualitative"
  ..@ factor     : Named int [1:5] 2 4 5 5 3
  .. ..- attr(*, "names")= chr [1:5] "gender" "eyebrow" "collar" "sub-caudal" ...
  ..@ nbCluster  : num 2
  ..@ knownLabels: int(0) 
  ..@ weight     : num(0) 
  ..@ nbVariable : int 5
  ..@ nbSample   : int 69
  ..@ criterion  : chr "BIC"
  ..@ models     :Formal class 'MultinomialModel' [package "Rmixmod"] with 5 slots
  .. .. ..@ variable.independency : logi(0) 
  .. .. ..@ component.independency: logi(0) 
  .. .. ..@ listModels            : chr "Binary_pk_Ekjh"
  .. .. ..@ free.proportions      : logi TRUE
  .. .. ..@ equal.proportions     : logi FALSE
  ..@ error      : logi FALSE
  ..@ results    :List of 1
  .. ..$ :Formal class 'MixmodResults' [package "Rmixmod"] with 9 slots
  .. .. .. ..@ nbCluster     : int 2
  .. .. .. ..@ model         : chr "Binary_pk_Ekjh"
  .. .. .. ..@ criterion     : chr "BIC"
  .. .. .. ..@ criterionValue: num 519
  .. .. .. ..@ parameters    :Formal class 'MultinomialParameter' [package "Rmixmod"] with 5 slots
  .. .. .. .. .. ..@ center     : num [1:2, 1:5] 2 1 2 3 2 1 2 1 1 1
  .. .. .. .. .. ..@ scatter    :List of 2
  .. .. .. .. .. .. ..$ : num [1:5, 1:5] 0.428 0.1203 0.0509 0.3641 0.1074 ...
  .. .. .. .. .. .. ..$ : num [1:5, 1:5] 0.4937 0.0761 0.1521 0.039 0.0577 ...
  .. .. .. .. .. ..@ factor     : num [1:5] 2 4 5 5 3
  .. .. .. .. .. ..@ proportions: num [1:2] 0.346 0.654
  .. .. .. .. .. ..@ nbFreeParam: num 29
  .. .. .. ..@ likelihood    : num -198
  .. .. .. ..@ partition     : num [1:69] 1 1 2 2 2 2 2 1 1 1 ...
  .. .. .. ..@ proba         : num [1:69, 1:2] 0.998001 0.998377 0.000245 0.019384 0.019384 ...
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
data(xem)
## maybe str(xem) ; plot(xem) ...
}
\keyword{datasets}
