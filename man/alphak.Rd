\name{alphak}
\alias{alphak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
alphak(tkx, Xijh, Y, Xquali)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tkx}{
%%     ~~Describe \code{tkx} here~~
}
  \item{Xijh}{
%%     ~~Describe \code{Xijh} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
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
function (tkx, Xijh, Y, Xquali) 
{
    alpha_k <- as.list(unique(Y))
    for (k in 1:length(alpha_k)) {
        variables <- as.list(seq_len(length(Xijh)))
        for (j in 1:length(variables)) {
            X <- vector("list", length(levels(Xquali[, j])))
            h_index <- 1
            for (h in sort(levels(Xquali[, j]))) {
                s <- 0
                for (i in 1:nrow(Xquali)) {
                  s <- s + tkx[i, k] * Xijh[[j]][i, -ncol(Xijh[[j]])][h_index]
                }
                X[[h_index]] <- s
                h_index <- h_index + 1
            }
            variables[[j]] <- data.frame(proba = unlist(X)/sum(unlist(X)))
        }
        alpha_k[[k]] <- variables
    }
    return(alpha_k)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
