\name{mod1}
\alias{mod1}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("mod1")}
\format{
  The format is:
List of 15
 $ call          : language Mclust(data = iris[, 1:4], G = 5)
 $ data          : num [1:150, 1:4] 5.1 4.9 4.7 4.6 5 5.4 4.6 5 4.4 4.9 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
 $ modelName     : chr "EEE"
 $ n             : int 150
 $ d             : int 4
 $ G             : int 5
 $ BIC           : 'mclustBIC' num [1, 1:14] -783 -743 -742 -688 -767 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr "5"
  .. ..$ : chr [1:14] "EII" "VII" "EEI" "VEI" ...
  ..- attr(*, "G")= num 5
  ..- attr(*, "modelNames")= chr [1:14] "EII" "VII" "EEI" "VEI" ...
  ..- attr(*, "control")=List of 4
  .. ..$ eps     : num 2.22e-16
  .. ..$ tol     : num [1:2] 1.00e-05 1.49e-08
  .. ..$ itmax   : int [1:2] 2147483647 2147483647
  .. ..$ equalPro: logi FALSE
  ..- attr(*, "initialization")=List of 3
  .. ..$ hcPairs: 'hc' num [1:2, 1:149] 102 143 128 139 18 41 1 18 1 28 ...
  .. .. ..- attr(*, "initialPartition")= num [1:150] 1 2 3 4 5 6 7 8 9 10 ...
  .. .. .. ..- attr(*, "unique")= int 150
  .. .. ..- attr(*, "dimensions")= int [1:2] 150 4
  .. .. ..- attr(*, "modelName")= chr "VVV"
  .. .. ..- attr(*, "call")= language hc(data = data, modelName = mclust.options("hcModelName")[1])
  .. .. ..- attr(*, "use")= chr "SVD"
  .. .. ..- attr(*, "data")= num [1:150, 1:4] -0.494 -0.454 -0.516 -0.502 -0.522 ...
  .. ..$ subset : NULL
  .. ..$ noise  : NULL
  ..- attr(*, "warn")= logi FALSE
  ..- attr(*, "n")= int 150
  ..- attr(*, "d")= int 4
  ..- attr(*, "oneD")= logi FALSE
  ..- attr(*, "criterion")= chr "BIC"
  ..- attr(*, "returnCodes")= num [1, 1:14] 0 0 0 0 0 0 0 0 -1 0 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr "5"
  .. .. ..$ : chr [1:14] "EII" "VII" "EEI" "VEI" ...
 $ bic           : num -605
 $ loglik        : num -217
 $ df            : num 34
 $ hypvol        : num NA
 $ parameters    :List of 3
  ..$ pro     : num [1:5] 0.0578 0.2755 0.3279 0.2301 0.1087
  ..$ mean    : num [1:4, 1:5] 5.41 3.92 1.432 0.256 4.921 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. ..$ : NULL
  ..$ variance:List of 6
  .. ..$ modelName: chr "EEE"
  .. ..$ d        : int 4
  .. ..$ G        : int 5
  .. ..$ sigma    : num [1:4, 1:4, 1:5] 0.2242 0.0793 0.1406 0.0515 0.0793 ...
  .. .. ..- attr(*, "dimnames")=List of 3
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : NULL
  .. ..$ Sigma    : num [1:4, 1:4] 0.2242 0.0793 0.1406 0.0515 0.0793 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ cholSigma: num [1:4, 1:4] 0.474 0 0 0 0.168 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
 $ z             : num [1:150, 1:5] 0.095949 0.001418 0.00498 0.000645 0.13493 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : NULL
 $ classification: num [1:150] 2 2 2 2 2 2 2 2 2 2 ...
 $ uncertainty   : num [1:150] 0.095949 0.001418 0.00498 0.000645 0.13493 ...
 - attr(*, "class")= chr "Mclust"
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
data(mod1)
## maybe str(mod1) ; plot(mod1) ...
}
\keyword{datasets}
