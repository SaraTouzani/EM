\name{teta_sans_init}
\alias{teta_sans_init}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("teta_sans_init")}
\format{
  The format is:
List of 3
 $ :List of 11
  ..$ iteration     : num 1
  ..$ Mean          : num [1:3, 1:5] 1 2 3 5.01 5.62 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.924 -0.511 2.642 1.074 -0.511 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.476 -0.414 1.827 0.755 -0.414 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.2566 0.0646 0.3074 0.1071 0.0646 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.333 0.26 0.407
  ..$ Alphak        :List of 3
  .. ..$ :List of 2
  .. .. ..$ :'data.frame':	3 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:3] 1 0 0
  .. .. ..$ :'data.frame':	4 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:4] 0 0 0.44 0.56
  .. ..$ :List of 2
  .. .. ..$ :'data.frame':	3 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:3] 0 0.769 0.231
  .. .. ..$ :'data.frame':	4 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:4] 0.00 1.96e-06 8.97e-01 1.03e-01
  .. ..$ :List of 2
  .. .. ..$ :'data.frame':	3 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:3] 0 0.328 0.672
  .. .. ..$ :'data.frame':	4 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:4] 0.197 0.803 0 0
  ..$ log_like      : num -544
  ..$ Q             : num [1, 1] 64.7
  ..$ Bic           : num 1384
  ..$ Icl           : num 1280
  ..$ proba_individu: num [1:150, 1:3] 0.376 0.183 0.372 0.227 0.444 ...
  ..$ part_MAP      : num [1:150, 1:3] 1 1 1 1 1 1 1 1 1 1 ...
 $ :List of 11
  ..$ iteration     : num 2
  ..$ Mean          : num [1:3, 1:5] 1 2 3 5.01 5.62 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.0951 0.042 0.1047 0.0385 0.042 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.2244 0.0493 0.2322 0.067 0.0493 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.333 0.26 0.407
  ..$ Alphak        :List of 3
  .. ..$ :List of 2
  .. .. ..$ :'data.frame':	3 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:3] 1 0 0
  .. .. ..$ :'data.frame':	4 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:4] 0 0 0.44 0.56
  .. ..$ :List of 2
  .. .. ..$ :'data.frame':	3 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:3] 0 0.769 0.231
  .. .. ..$ :'data.frame':	4 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:4] 0.00 1.43e-07 8.97e-01 1.03e-01
  .. ..$ :List of 2
  .. .. ..$ :'data.frame':	3 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:3] 0 0.328 0.672
  .. .. ..$ :'data.frame':	4 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:4] 0.197 0.803 0 0
  ..$ log_like      : num -362
  ..$ Q             : num [1, 1] 125
  ..$ Bic           : num 1019
  ..$ Icl           : num 915
  ..$ proba_individu: num [1:150, 1:3] 6.35 3.51 5.28 4.25 6.89 ...
  ..$ part_MAP      : num [1:150, 1:3] 1 1 1 1 1 1 1 1 1 1 ...
 $ :List of 11
  ..$ iteration     : num 3
  ..$ Mean          : num [1:3, 1:5] 1 2 3 5.01 5.62 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.0951 0.042 0.1047 0.0385 0.042 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.2244 0.0493 0.2322 0.067 0.0493 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.333 0.26 0.407
  ..$ Alphak        :List of 3
  .. ..$ :List of 2
  .. .. ..$ :'data.frame':	3 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:3] 1 0 0
  .. .. ..$ :'data.frame':	4 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:4] 0 0 0.44 0.56
  .. ..$ :List of 2
  .. .. ..$ :'data.frame':	3 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:3] 0 0.769 0.231
  .. .. ..$ :'data.frame':	4 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:4] 0.00 5.62e-08 8.97e-01 1.03e-01
  .. ..$ :List of 2
  .. .. ..$ :'data.frame':	3 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:3] 0 0.328 0.672
  .. .. ..$ :'data.frame':	4 obs. of  1 variable:
  .. .. .. ..$ proba: num [1:4] 0.197 0.803 0 0
  ..$ log_like      : num -362
  ..$ Q             : num [1, 1] 125
  ..$ Bic           : num 1019
  ..$ Icl           : num 915
  ..$ proba_individu: num [1:150, 1:3] 6.35 3.51 5.28 4.25 6.89 ...
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
data(teta_sans_init)
## maybe str(teta_sans_init) ; plot(teta_sans_init) ...
}
\keyword{datasets}
