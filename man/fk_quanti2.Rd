\name{fk_quanti2}
\alias{fk_quanti2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fk_quanti2(Xquanti, mean, sigma, Y, k = length(sigma), n = nrow(Xquanti), p = length(colnames(Xquanti)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xquanti}{
%%     ~~Describe \code{Xquanti} here~~
}
  \item{mean}{
%%     ~~Describe \code{mean} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
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
function (Xquanti, mean, sigma, Y, k = length(sigma), n = nrow(Xquanti), 
    p = length(colnames(Xquanti))) 
{
    K = length(unique(Y))
    fk_quanti_ik <- matrix(1, nrow = n, ncol = length(unique(Y)))
    fkx_q <- as.list(seq_len(length(Y)))
    for (i in 1:n) {
        for (k in 1:K) {
            fk_quanti_ik[i, k] <- (1/(((2 * pi)^(p/2)) * (det(as.matrix(sigma[[k]]))^(0.5)))) * 
                exp((Xquanti[i, ] - as.numeric(mean[k, -1])) \%*\% 
                  ginv(as.matrix(sigma[[k]])) \%*\% as.matrix(Xquanti[i, 
                  ] - as.numeric(mean[k, -1])) * (-0.5))
        }
        fkx_q[[i]] <- fk_quanti_ik[i, ]
    }
    return(fkx_q)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
