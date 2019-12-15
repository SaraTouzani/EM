\name{teta}
\alias{teta}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("teta")}
\format{
  The format is:
List of 3
 $ :List of 11
  ..$ iteration     : num 1
  ..$ Mean          : num [1:3, 1:5] 3 3 3 5.84 5.84 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.703 -0.052 1.306 0.53 -0.052 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.6977 -0.0463 1.2883 0.5219 -0.0463 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.6811 -0.0422 1.266 0.5129 -0.0422 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 1 1.11 1.24
  ..$ Alphak        :List of 1
  .. ..$ :List of 2
  .. .. ..$ :'data.frame':	3 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:3] 0.333 0.333 0.333
  .. .. ..$ :'data.frame':	4 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:4] 0.08 0.327 0.38 0.213
  ..$ log_like      : num -735
  ..$ Q             : num [1:3] 2.19 2.43 2.71
  ..$ Bic           : num 1564
  ..$ Icl           : num 1187
  ..$ proba_individu: num [1:150, 1] 0.02499 0.00959 0.01442 0.01224 0.01231 ...
  ..$ part_MAP      : num [1:150, 1:3] 1 1 1 1 1 1 1 1 1 1 ...
 $ :List of 11
  ..$ iteration     : num 2
  ..$ Mean          : num [1:3, 1:5] 3 3 3 5.84 5.84 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.6811 -0.0422 1.2658 0.5128 -0.0422 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.6811 -0.0422 1.2658 0.5128 -0.0422 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.6811 -0.0422 1.2658 0.5128 -0.0422 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 1 1.11 1.24
  ..$ Alphak        :List of 1
  .. ..$ :List of 2
  .. .. ..$ :'data.frame':	3 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:3] 0.333 0.333 0.333
  .. .. ..$ :'data.frame':	4 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:4] 0.08 0.327 0.38 0.213
  ..$ log_like      : num -734
  ..$ Q             : num [1:3] 2.23 2.47 2.75
  ..$ Bic           : num 1564
  ..$ Icl           : num 1187
  ..$ proba_individu: num [1:150, 1] 0.02539 0.00995 0.01464 0.01215 0.01209 ...
  ..$ part_MAP      : num [1:150, 1:3] 1 1 1 1 1 1 1 1 1 1 ...
 $ :List of 11
  ..$ iteration     : num 3
  ..$ Mean          : num [1:3, 1:5] 3 3 3 5.84 5.84 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.6811 -0.0422 1.2658 0.5128 -0.0422 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.6811 -0.0422 1.2658 0.5128 -0.0422 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.6811 -0.0422 1.2658 0.5128 -0.0422 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 1 1.11 1.24
  ..$ Alphak        :List of 1
  .. ..$ :List of 2
  .. .. ..$ :'data.frame':	3 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:3] 0.333 0.333 0.333
  .. .. ..$ :'data.frame':	4 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:4] 0.08 0.327 0.38 0.213
  ..$ log_like      : num -734
  ..$ Q             : num [1:3] 2.23 2.47 2.75
  ..$ Bic           : num 1564
  ..$ Icl           : num 1187
  ..$ proba_individu: num [1:150, 1] 0.02539 0.00995 0.01464 0.01215 0.01209 ...
  ..$ part_MAP      : num [1:150, 1:3] 1 1 1 1 1 1 1 1 1 1 ...
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
data(teta)
## maybe str(teta) ; plot(teta) ...
}
\keyword{datasets}
