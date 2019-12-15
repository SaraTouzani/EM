\name{Qt}
\alias{Qt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Qt(tkx, fkx, Pk)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tkx}{
%%     ~~Describe \code{tkx} here~~
}
  \item{fkx}{
%%     ~~Describe \code{fkx} here~~
}
  \item{Pk}{
%%     ~~Describe \code{Pk} here~~
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
function (tkx, fkx, Pk) 
{
    n <- nrow(tkx)
    K <- ncol(tkx)
    somme <- 0
    for (i in 1:n) {
        somme <- somme + ((log1p(Pk) + log1p(unlist(fkx[[i]]))) \%*\% 
            tkx[i, ])
    }
    return(somme)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
