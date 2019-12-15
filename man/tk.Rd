\name{tk}
\alias{tk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tk(fkx, pk)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fkx}{
%%     ~~Describe \code{fkx} here~~
}
  \item{pk}{
%%     ~~Describe \code{pk} here~~
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
function (fkx, pk) 
{
    K <- length(pk)
    n <- length(fkx)
    tk <- matrix(0, ncol = K, nrow = n)
    for (i in 1:n) {
        tk[i, ] <- (as.numeric(unlist(fkx[[i]])) * pk)/(as.numeric(unlist(fkx[[i]])) \%*\% 
            pk)[1, 1]
    }
    return(tk)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
