\name{fk_quali}
\alias{fk_quali}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fk_quali(Y, Xijh, alpha_ki, k = length(alpha_ki), n = nrow(Xquali), p = length(colnames(Xquali)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{Xijh}{
%%     ~~Describe \code{Xijh} here~~
}
  \item{alpha_ki}{
%%     ~~Describe \code{alpha_ki} here~~
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
function (Y, Xijh, alpha_ki, k = length(alpha_ki), n = nrow(Xquali), 
    p = length(colnames(Xquali))) 
{
    if (p == 0) {
        p <- 1
        n <- length(Xquali)
    }
    fk_quali_ik <- matrix(1, nrow = n, ncol = length(unique(Y)))
    K = length(unique(Y))
    fkx_q <- as.list(seq_len(length(Y)))
    for (i in 1:n) {
        for (k in 1:K) {
            for (j in 1:p) {
                fk_quali_ik[i, k] <- fk_quali_ik[i, k] * (as.numeric(alpha_ki[[k]][[j]][["proba"]] \%*\% 
                  as.numeric(Xijh[[j]][i, -length(Xijh[[j]][i, 
                    ])])))
            }
        }
        fkx_q[[i]] <- fk_quali_ik[i, ]
    }
    return(fkx_q)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
