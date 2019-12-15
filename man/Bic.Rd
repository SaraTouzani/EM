\name{Bic}
\alias{Bic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Bic(f, Y = NULL, nbr_par = NULL, Xquali = NULL, Xquanti = NULL, K = length(unique(Y)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{
%%     ~~Describe \code{f} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{nbr_par}{
%%     ~~Describe \code{nbr_par} here~~
}
  \item{Xquali}{
%%     ~~Describe \code{Xquali} here~~
}
  \item{Xquanti}{
%%     ~~Describe \code{Xquanti} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
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
function (f, Y = NULL, nbr_par = NULL, Xquali = NULL, Xquanti = NULL, 
    K = length(unique(Y))) 
{
    nbr_par <- 0
    if (!is.null(Xquali)) {
        p <- ncol(Xquali)
        for (j in 1:p) {
            nbr_par <- nbr_par + length(unique(Xquali[, j])) - 
                1
        }
        nbr_par <- nbr_par * K
    }
    if (!is.null(Xquanti)) {
        p <- ncol(Xquanti)
        nbr_par <- K * (p * (p + 1)/2 + p + 1) - 1 + nbr_par
    }
    return(-2 * sum(log(f)) + nbr_par * log(length(Y)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
