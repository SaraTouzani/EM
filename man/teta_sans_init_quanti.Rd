\name{teta_sans_init_quanti}
\alias{teta_sans_init_quanti}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("teta_sans_init_quanti")}
\format{
  The format is:
List of 23
 $ :List of 10
  ..$ iteration     : num 1
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.88 6.13 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.1956 0.0765 0.1094 0.015 0.0765 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.439 0.116 0.458 0.164 0.116 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.117 0.55 0.333
  ..$ log_like      : num -211
  ..$ Q             : num [1, 1] 163
  ..$ Bic           : num 643
  ..$ Icl           : num 547
  ..$ proba_individu: num [1:150, 1:3] 2.07e-117 1.32e-104 7.63e-110 2.76e-98 3.68e-118 ...
  ..$ part_MAP      : num [1:150, 1:3] 5.02e-119 7.38e-106 2.83e-111 1.27e-99 1.05e-119 ...
 $ :List of 10
  ..$ iteration     : num 2
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.84 6.15 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.28997 0.10402 0.18917 -0.00296 0.10402 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.386 0.0928 0.3865 0.1368 0.0928 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.111 0.556 0.333
  ..$ log_like      : num -204
  ..$ Q             : num [1, 1] 164
  ..$ Bic           : num 628
  ..$ Icl           : num 533
  ..$ proba_individu: num [1:150, 1:3] 3.72e-114 3.56e-99 5.15e-105 7.34e-94 2.01e-115 ...
  ..$ part_MAP      : num [1:150, 1:3] 8.56e-116 1.89e-100 1.82e-106 3.21e-95 5.42e-117 ...
 $ :List of 10
  ..$ iteration     : num 3
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.79 6.14 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.33585 0.11895 0.22567 -0.00879 0.11895 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.3784 0.0849 0.379 0.1348 0.0849 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.126 0.541 0.333
  ..$ log_like      : num -201
  ..$ Q             : num [1, 1] 165
  ..$ Bic           : num 623
  ..$ Icl           : num 527
  ..$ proba_individu: num [1:150, 1:3] 3.88e-107 6.77e-92 2.92e-98 3.10e-88 5.87e-109 ...
  ..$ part_MAP      : num [1:150, 1:3] 1.02e-108 4.07e-93 1.17e-99 1.54e-89 1.80e-110 ...
 $ :List of 10
  ..$ iteration     : num 4
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.75 6.13 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.33436 0.11954 0.22039 -0.00237 0.11954 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.3813 0.0819 0.382 0.1323 0.0819 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.142 0.525 0.333
  ..$ log_like      : num -199
  ..$ Q             : num [1, 1] 165
  ..$ Bic           : num 618
  ..$ Icl           : num 522
  ..$ proba_individu: num [1:150, 1:3] 5.06e-102 1.62e-86 2.52e-93 6.33e-84 2.99e-104 ...
  ..$ part_MAP      : num [1:150, 1:3] 1.49e-103 1.09e-87 1.13e-94 3.54e-85 1.03e-105 ...
 $ :List of 10
  ..$ iteration     : num 5
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.71 6.12 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.3286 0.1198 0.21 0.0062 0.1198 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.3876 0.0806 0.3884 0.1312 0.0806 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.157 0.51 0.333
  ..$ log_like      : num -198
  ..$ Q             : num [1, 1] 165
  ..$ Bic           : num 616
  ..$ Icl           : num 519
  ..$ proba_individu: num [1:150, 1:3] 2.38e-100 3.49e-84 1.85e-91 4.39e-82 7.06e-103 ...
  ..$ part_MAP      : num [1:150, 1:3] 7.77e-102 2.62e-85 9.24e-93 2.72e-83 2.70e-104 ...
 $ :List of 10
  ..$ iteration     : num 6
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.67 6.12 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.3254 0.1204 0.2037 0.0116 0.1204 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.3947 0.0799 0.3957 0.1314 0.0799 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.171 0.496 0.333
  ..$ log_like      : num -196
  ..$ Q             : num [1, 1] 165
  ..$ Bic           : num 613
  ..$ Icl           : num 517
  ..$ proba_individu: num [1:150, 1:3] 1.28e-98 5.72e-82 1.27e-89 2.59e-80 2.13e-101 ...
  ..$ part_MAP      : num [1:150, 1:3] 4.53e-100 4.67e-83 6.88e-91 1.75e-81 8.87e-103 ...
 $ :List of 10
  ..$ iteration     : num 7
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.64 6.12 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.3238 0.1216 0.2006 0.0158 0.1216 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.4022 0.0795 0.4039 0.1324 0.0795 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.182 0.484 0.333
  ..$ log_like      : num -196
  ..$ Q             : num [1, 1] 166
  ..$ Bic           : num 612
  ..$ Icl           : num 514
  ..$ proba_individu: num [1:150, 1:3] 2.04e-96 1.60e-79 2.09e-87 3.08e-78 2.27e-99 ...
  ..$ part_MAP      : num [1:150, 1:3] 7.75e-98 1.40e-80 1.21e-88 2.21e-79 1.01e-100 ...
 $ :List of 10
  ..$ iteration     : num 8
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.61 6.12 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.3243 0.1239 0.2001 0.0207 0.1239 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.4095 0.0796 0.4128 0.1338 0.0796 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.192 0.475 0.333
  ..$ log_like      : num -195
  ..$ Q             : num [1, 1] 166
  ..$ Bic           : num 610
  ..$ Icl           : num 513
  ..$ proba_individu: num [1:150, 1:3] 2.30e-94 2.24e-77 2.48e-85 2.85e-76 2.12e-97 ...
  ..$ part_MAP      : num [1:150, 1:3] 9.20e-96 2.06e-78 1.51e-86 2.16e-77 9.95e-99 ...
 $ :List of 10
  ..$ iteration     : num 9
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.59 6.12 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.3283 0.1269 0.2024 0.0266 0.1269 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.4153 0.0798 0.4215 0.1354 0.0798 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.201 0.466 0.333
  ..$ log_like      : num -194
  ..$ Q             : num [1, 1] 167
  ..$ Bic           : num 609
  ..$ Icl           : num 511
  ..$ proba_individu: num [1:150, 1:3] 3.82e-92 3.12e-75 4.36e-83 4.10e-74 3.71e-95 ...
  ..$ part_MAP      : num [1:150, 1:3] 1.59e-93 3.00e-76 2.78e-84 3.25e-75 1.81e-96 ...
 $ :List of 10
  ..$ iteration     : num 10
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.57 6.12 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.3371 0.1306 0.2078 0.0336 0.1306 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.4187 0.0799 0.4293 0.1369 0.0799 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.208 0.459 0.333
  ..$ log_like      : num -194
  ..$ Q             : num [1, 1] 167
  ..$ Bic           : num 608
  ..$ Icl           : num 510
  ..$ proba_individu: num [1:150, 1:3] 2.35e-89 1.03e-72 2.69e-80 1.95e-71 2.98e-92 ...
  ..$ part_MAP      : num [1:150, 1:3] 1.02e-90 1.03e-73 1.78e-81 1.60e-72 1.51e-93 ...
 $ :List of 10
  ..$ iteration     : num 11
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.54 6.13 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.3502 0.1345 0.2158 0.0413 0.1345 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.4197 0.0799 0.4362 0.1382 0.0799 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.216 0.451 0.333
  ..$ log_like      : num -193
  ..$ Q             : num [1, 1] 167
  ..$ Bic           : num 606
  ..$ Icl           : num 508
  ..$ proba_individu: num [1:150, 1:3] 3.96e-86 6.56e-70 4.19e-77 2.16e-68 7.54e-89 ...
  ..$ part_MAP      : num [1:150, 1:3] 1.77e-87 6.76e-71 2.88e-78 1.84e-69 3.96e-90 ...
 $ :List of 10
  ..$ iteration     : num 12
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.52 6.13 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.361 0.1369 0.2225 0.0482 0.1369 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.4206 0.0803 0.4427 0.1394 0.0803 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.222 0.444 0.333
  ..$ log_like      : num -192
  ..$ Q             : num [1, 1] 167
  ..$ Bic           : num 605
  ..$ Icl           : num 507
  ..$ proba_individu: num [1:150, 1:3] 6.30e-83 3.47e-67 5.41e-74 1.82e-65 1.79e-85 ...
  ..$ part_MAP      : num [1:150, 1:3] 2.91e-84 3.69e-68 3.83e-75 1.60e-66 9.69e-87 ...
 $ :List of 10
  ..$ iteration     : num 13
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.51 6.13 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.3632 0.1365 0.224 0.0525 0.1365 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.4238 0.0815 0.4495 0.1406 0.0815 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.228 0.438 0.333
  ..$ log_like      : num -192
  ..$ Q             : num [1, 1] 168
  ..$ Bic           : num 604
  ..$ Icl           : num 506
  ..$ proba_individu: num [1:150, 1:3] 5.96e-80 1.10e-64 3.52e-71 7.03e-63 2.26e-82 ...
  ..$ part_MAP      : num [1:150, 1:3] 2.83e-81 1.20e-65 2.56e-72 6.34e-64 1.26e-83 ...
 $ :List of 10
  ..$ iteration     : num 14
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.5 6.13 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.359 0.134 0.222 0.055 0.134 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.4289 0.0835 0.4565 0.1419 0.0835 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.233 0.434 0.333
  ..$ log_like      : num -192
  ..$ Q             : num [1, 1] 168
  ..$ Bic           : num 604
  ..$ Icl           : num 505
  ..$ proba_individu: num [1:150, 1:3] 3.33e-77 2.24e-62 1.26e-68 1.43e-60 1.55e-79 ...
  ..$ part_MAP      : num [1:150, 1:3] 1.61e-78 2.50e-63 9.31e-70 1.32e-61 8.79e-81 ...
 $ :List of 10
  ..$ iteration     : num 15
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.49 6.14 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.3535 0.1311 0.2186 0.0567 0.1311 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.435 0.086 0.464 0.143 0.086 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.236 0.431 0.333
  ..$ log_like      : num -192
  ..$ Q             : num [1, 1] 168
  ..$ Bic           : num 604
  ..$ Icl           : num 505
  ..$ proba_individu: num [1:150, 1:3] 9.91e-75 2.77e-60 2.46e-66 1.65e-58 5.35e-77 ...
  ..$ part_MAP      : num [1:150, 1:3] 4.87e-76 3.13e-61 1.85e-67 1.53e-59 3.08e-78 ...
 $ :List of 10
  ..$ iteration     : num 16
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.48 6.14 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.3466 0.1274 0.2145 0.0582 0.1274 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.4421 0.0892 0.4721 0.1454 0.0892 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.239 0.428 0.333
  ..$ log_like      : num -191
  ..$ Q             : num [1, 1] 169
  ..$ Bic           : num 603
  ..$ Icl           : num 504
  ..$ proba_individu: num [1:150, 1:3] 1.49e-72 1.99e-58 2.63e-64 1.09e-56 9.12e-75 ...
  ..$ part_MAP      : num [1:150, 1:3] 7.41e-74 2.27e-59 2.00e-65 1.03e-57 5.30e-76 ...
 $ :List of 10
  ..$ iteration     : num 17
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.47 6.14 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.3382 0.1227 0.209 0.0597 0.1227 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.4506 0.0936 0.482 0.1481 0.0936 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.24 0.426 0.333
  ..$ log_like      : num -191
  ..$ Q             : num [1, 1] 169
  ..$ Bic           : num 603
  ..$ Icl           : num 504
  ..$ proba_individu: num [1:150, 1:3] 1.24e-70 8.96e-57 1.69e-62 4.58e-55 8.38e-73 ...
  ..$ part_MAP      : num [1:150, 1:3] 6.18e-72 1.03e-57 1.29e-63 4.34e-56 4.91e-74 ...
 $ :List of 10
  ..$ iteration     : num 18
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.46 6.15 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.3265 0.116 0.2008 0.0615 0.116 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.462 0.1 0.495 0.152 0.1 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.241 0.426 0.333
  ..$ log_like      : num -191
  ..$ Q             : num [1, 1] 169
  ..$ Bic           : num 602
  ..$ Icl           : num 503
  ..$ proba_individu: num [1:150, 1:3] 6.64e-69 3.00e-55 7.84e-61 1.45e-53 4.90e-71 ...
  ..$ part_MAP      : num [1:150, 1:3] 3.32e-70 3.46e-56 6.02e-62 1.38e-54 2.88e-72 ...
 $ :List of 10
  ..$ iteration     : num 19
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.44 6.16 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.3089 0.1056 0.1873 0.0636 0.1056 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.478 0.11 0.515 0.158 0.11 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.24 0.426 0.333
  ..$ log_like      : num -190
  ..$ Q             : num [1, 1] 170
  ..$ Bic           : num 601
  ..$ Icl           : num 501
  ..$ proba_individu: num [1:150, 1:3] 3.06e-67 1.00e-53 3.63e-59 4.58e-52 2.41e-69 ...
  ..$ part_MAP      : num [1:150, 1:3] 1.53e-68 1.15e-54 2.78e-60 4.35e-53 1.41e-70 ...
 $ :List of 10
  ..$ iteration     : num 20
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.42 6.18 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.2853 0.0914 0.1677 0.0659 0.0914 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.499 0.123 0.54 0.166 0.123 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.239 0.428 0.333
  ..$ log_like      : num -189
  ..$ Q             : num [1, 1] 170
  ..$ Bic           : num 598
  ..$ Icl           : num 498
  ..$ proba_individu: num [1:150, 1:3] 2.45e-65 5.66e-52 3.02e-57 2.43e-50 2.03e-67 ...
  ..$ part_MAP      : num [1:150, 1:3] 1.21e-66 6.46e-53 2.29e-58 2.29e-51 1.18e-68 ...
 $ :List of 10
  ..$ iteration     : num 21
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.4 6.19 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.2688 0.0817 0.1534 0.068 0.0817 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.512 0.131 0.558 0.172 0.131 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.238 0.429 0.333
  ..$ log_like      : num -188
  ..$ Q             : num [1, 1] 171
  ..$ Bic           : num 597
  ..$ Icl           : num 497
  ..$ proba_individu: num [1:150, 1:3] 4.74e-63 5.14e-50 4.46e-55 2.10e-48 4.39e-65 ...
  ..$ part_MAP      : num [1:150, 1:3] 2.35e-64 5.86e-51 3.38e-56 1.98e-49 2.55e-66 ...
 $ :List of 10
  ..$ iteration     : num 22
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.39 6.19 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.2646 0.0796 0.1506 0.0695 0.0796 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.516 0.133 0.564 0.174 0.133 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.239 0.428 0.333
  ..$ log_like      : num -188
  ..$ Q             : num [1, 1] 171
  ..$ Bic           : num 597
  ..$ Icl           : num 496
  ..$ proba_individu: num [1:150, 1:3] 5.98e-61 2.51e-48 3.46e-53 1.03e-46 6.58e-63 ...
  ..$ part_MAP      : num [1:150, 1:3] 2.97e-62 2.87e-49 2.64e-54 9.69e-48 3.83e-64 ...
 $ :List of 10
  ..$ iteration     : num 23
  ..$ Mean          : num [1:3, 1:5] 1 2 3 6.39 6.19 ...
  ..$ Sigma         :List of 3
  .. ..$ : num [1:4, 1:4] 0.2637 0.0789 0.1506 0.0701 0.0789 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.518 0.133 0.567 0.175 0.133 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. ..$ : num [1:4, 1:4] 0.1218 0.0972 0.016 0.0101 0.0972 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  .. .. .. ..$ : chr [1:4] "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"
  ..$ Proba         : num [1:3] 0.24 0.427 0.333
  ..$ log_like      : num -188
  ..$ Q             : num [1, 1] 172
  ..$ Bic           : num 596
  ..$ Icl           : num 496
  ..$ proba_individu: num [1:150, 1:3] 2.54e-59 4.68e-47 9.46e-52 2.01e-45 3.30e-61 ...
  ..$ part_MAP      : num [1:150, 1:3] 1.27e-60 5.38e-48 7.23e-53 1.91e-46 1.94e-62 ...
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
data(teta_sans_init_quanti)
## maybe str(teta_sans_init_quanti) ; plot(teta_sans_init_quanti) ...
}
\keyword{datasets}
