\name{mod}
\alias{mod}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("mod")}
\format{
  The format is:
List of 15
 $ call          : language Mclust(data = iris[, 1:4])
 $ data          : num [1:150, 1:4] 5.1 4.9 4.7 4.6 5 5.4 4.6 5 4.4 4.9 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
 $ modelName     : chr "VEV"
 $ n             : int 150
 $ d             : int 4
 $ G             : int 2
 $ BIC           : 'mclustBIC' num [1:9, 1:14] -1804 -1123 -879 -894 -783 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:9] "1" "2" "3" "4" ...
  .. ..$ : chr [1:14] "EII" "VII" "EEI" "VEI" ...
  ..- attr(*, "G")= num [1:9] 1 2 3 4 5 6 7 8 9
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
  ..- attr(*, "returnCodes")= num [1:9, 1:14] 0 0 0 0 0 0 0 0 0 0 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:9] "1" "2" "3" "4" ...
  .. .. ..$ : chr [1:14] "EII" "VII" "EEI" "VEI" ...
 $ bic           : num -562
 $ loglik        : num -216
 $ df            : num 26
 $ hypvol        : num NA
 $ parameters    :List of 4
  ..$ pro     : num [1:2] 0.333 0.667
  ..$ mean    : num [1:4, 1:2] 5.006 3.428 1.462 0.246 6.262 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. ..$ : NULL
  ..$ variance:List of 7
  .. ..$ modelName  : chr "VEV"
  .. ..$ d          : int 4
  .. ..$ G          : int 2
  .. ..$ sigma      : num [1:4, 1:4, 1:2] 0.1507 0.1308 0.0208 0.0131 0.1308 ...
  .. .. ..- attr(*, "dimnames")=List of 3
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : NULL
  .. ..$ scale      : num [1:2] 0.0377 0.1331
  .. ..$ shape      : num [1:4] 7.911 0.923 0.63 0.217
  .. ..$ orientation: num [1:4, 1:4, 1:2] -0.6691 -0.7341 -0.0965 -0.0636 0.5979 ...
  .. .. ..- attr(*, "dimnames")=List of 3
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : NULL
  ..$ Vinv    : NULL
 $ z             : num [1:150, 1:2] 1 1 1 1 1 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : NULL
 $ classification: num [1:150] 1 1 1 1 1 1 1 1 1 1 ...
 $ uncertainty   : num [1:150] 2.51e-11 5.56e-08 3.64e-09 8.61e-08 8.50e-12 ...
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
data(mod)
## maybe str(mod) ; plot(mod) ...
}
\keyword{datasets}
