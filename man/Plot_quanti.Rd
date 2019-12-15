\name{Plot_quanti}
\alias{Plot_quanti}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Plot_quanti(Dataquanti, Em_results, all = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Dataquanti}{
%%     ~~Describe \code{Dataquanti} here~~
}
  \item{Em_results}{
%%     ~~Describe \code{Em_results} here~~
}
  \item{all}{
%%     ~~Describe \code{all} here~~
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
function (Dataquanti, Em_results, all = TRUE) 
{
    tkx <- Em_results[[length(Em_results)]][["part_MAP"]]
    Y <- apply(tkx, MARGIN = 1, which.max)
    if (all) {
        pairs.panels(Dataquanti, gap = 0, bg = Y, pch = 21)
    }
    else {
        i <- 1
        while (i < ((ncol(Dataquanti)) - 6)) {
            pairs.panels(Dataquanti[, i:(i + 5)], gap = 0, bg = Y, 
                pch = 21)
            i <- i + 6
        }
        pairs.panels(Dataquanti[, i:ncol(Dataquanti)], gap = 0, 
            bg = Y, pch = 21)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
