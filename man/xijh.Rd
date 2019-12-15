\name{xijh}
\alias{xijh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xijh(Xquali)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xquali}{
%%     ~~Describe \code{Xquali} here~~
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
function (Xquali) 
{
    if (is.null(dim(Xquali))) {
        p <- 1
        Xijh <- vector("list", 1)
        Xquali <- as.matrix(Xquali)
    }
    else {
        p <- length(colnames(Xquali))
        Xijh <- as.list(colnames(Xquali))
    }
    for (j in 1:p) {
        Xih <- matrix(0, ncol = length(levels(Xquali[, j])), 
            nrow = nrow(Xquali))
        modalite_ij <- matrix("", nrow = nrow(Xquali))
        h_index <- 1
        for (h in sort(levels(Xquali[, j]))) {
            for (i in 1:nrow(Xquali)) {
                if ((Xquali[i, j] == h)) {
                  Xih[i, h_index] <- 1
                  modalite_ij[i] <- h
                }
                else {
                  Xih[i, h_index] <- 0
                }
            }
            h_index <- h_index + 1
        }
        Xijh[[j]] <- data.frame(proba = Xih, modalite = modalite_ij)
    }
    return(Xijh)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
