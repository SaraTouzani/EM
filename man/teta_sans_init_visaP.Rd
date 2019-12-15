\name{teta_sans_init_visaP}
\alias{teta_sans_init_visaP}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("teta_sans_init_visaP")}
\format{
  The format is:
List of 2
 $ :List of 10
  ..$ iteration     : num 1
  ..$ Mean          : num [1:2, 1:8] 1 2 -0.0874 11.6251 0.037 ...
  ..$ Sigma         :List of 2
  .. ..$ : num [1:7, 1:7] 5.26596 0.44638 -1.16809 0.16202 -0.00443 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:7] "V1" "V2" "V3" "V4" ...
  .. .. .. ..$ : chr [1:7] "V1" "V2" "V3" "V4" ...
  .. ..$ : num [1:7, 1:7] 26.33 -1.71 5.87 -8.88 2.61 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:7] "V1" "V2" "V3" "V4" ...
  .. .. .. ..$ : chr [1:7] "V1" "V2" "V3" "V4" ...
  ..$ Proba         : num [1:2] 0.99254 0.00746
  ..$ log_like      : num -14973
  ..$ Q             : num [1, 1] 734
  ..$ Bic           : num 30442
  ..$ Icl           : num 29699
  ..$ proba_individu: num [1:1072, 1:2] 1.11e-09 2.03e-06 2.46e-07 1.25e-11 3.43e-06 ...
  ..$ part_MAP      : num [1:1072, 1:2] 1 1 1 1 1 ...
 $ :List of 10
  ..$ iteration     : num 2
  ..$ Mean          : num [1:2, 1:8] 1 2 -0.0874 11.6253 0.037 ...
  ..$ Sigma         :List of 2
  .. ..$ : num [1:7, 1:7] 5.26596 0.44639 -1.16809 0.16202 -0.00443 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:7] "V2" "V3" "V4" "V5" ...
  .. .. .. ..$ : chr [1:7] "V2" "V3" "V4" "V5" ...
  .. ..$ : num [1:7, 1:7] 26.33 -1.71 5.87 -8.88 2.61 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:7] "V2" "V3" "V4" "V5" ...
  .. .. .. ..$ : chr [1:7] "V2" "V3" "V4" "V5" ...
  ..$ Proba         : num [1:2] 0.99254 0.00746
  ..$ log_like      : num -14973
  ..$ Q             : num [1, 1] 734
  ..$ Bic           : num 30442
  ..$ Icl           : num 29699
  ..$ proba_individu: num [1:1072, 1:2] 1.11e-09 2.03e-06 2.46e-07 1.25e-11 3.43e-06 ...
  ..$ part_MAP      : num [1:1072, 1:2] 1 1 1 1 1 ...
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
data(teta_sans_init_visaP)
## maybe str(teta_sans_init_visaP) ; plot(teta_sans_init_visaP) ...
}
\keyword{datasets}
