\name{meank}
\alias{meank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
meank(Nk, tkx, Xquanti, Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Nk}{
%%     ~~Describe \code{Nk} here~~
}
  \item{tkx}{
%%     ~~Describe \code{tkx} here~~
}
  \item{Xquanti}{
%%     ~~Describe \code{Xquanti} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Nk, tkx, Xquanti, Y) 
{
    K <- length(Nk)
    p <- ncol(Xquanti)
    moyenne <- matrix(0, ncol = p, nrow = K)
    moyenne[, 1] <- unique(Y)
    for (k in 1:K) {
        for (j in 1:p) {
            moyenne[k, j] <- sum((tkx[, k] * Xquanti)[, j])/Nk[k]
        }
    }
    moyenne <- cbind(sort(unique(Y)), moyenne)
    return(moyenne)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
