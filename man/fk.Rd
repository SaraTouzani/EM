\name{fk}
\alias{fk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fk(Xquanti, mean, sigma, Y, Xijh, alpha_k, K = length(sigma), n = nrow(Xquanti), p = length(colnames(Xquali)))
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
  \item{Xijh}{
%%     ~~Describe \code{Xijh} here~~
}
  \item{alpha_k}{
%%     ~~Describe \code{alpha_k} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
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
function (Xquanti, mean, sigma, Y, Xijh, alpha_k, K = length(sigma), 
    n = nrow(Xquanti), p = length(colnames(Xquali))) 
{
    if (p == 0) {
        p <- 1
    }
    fk_quali_ik <- matrix(1, nrow = n, ncol = length(unique(Y)))
    K = length(unique(Y))
    for (i in 1:n) {
        for (k in 1:K) {
            for (j in 1:p) {
                fk_quali_ik[i, k] <- fk_quali_ik[i, k] * (as.numeric(alpha_k[[k]][[j]][["proba"]] \%*\% 
                  as.numeric(Xijh[[j]][i, -length(Xijh[[j]][i, 
                    ])])))
            }
        }
    }
    fk_quati_ik <- matrix(1, nrow = n, ncol = length(unique(Y)))
    p = length(colnames(Xquanti))
    if (p == 0) {
        p <- 1
    }
    for (i in 1:n) {
        for (k in 1:K) {
            fk_quati_ik[i, k] <- (1/(((2 * pi)^(p/2)) * (det(as.matrix(sigma[[k]]))^(0.5)))) * 
                exp((-0.5) * as.matrix(Xquanti[i, ] - as.numeric(mean[k, 
                  -1])) \%*\% ginv(as.matrix(sigma[[k]])) \%*\% t(as.matrix(Xquanti[i, 
                  ] - as.numeric(mean[k, -1]))))
        }
    }
    fkxi <- as.list(seq_len(nrow(Xquanti)))
    for (i in 1:n) {
        fkxi[[i]] <- data.frame(densite = fk_quali_ik[i, ] * 
            fk_quati_ik[i, ])
    }
    return(fkxi)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
