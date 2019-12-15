\name{Em_quali}
\alias{Em_quali}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do function to do de clustering on categorial futures.... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Em_quali(Xquali, alpha_ki = NULL, Xijh_i = NULL, pki = NULL, Y = NULL, methode = "init", e = 0.1, nbr_iteration = NULL, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xquali}{
%%     ~~Describe \code{Xquali} here~~
}
  \item{alpha_ki}{
%%     ~~Describe \code{alpha_ki} here~~
}
  \item{Xijh_i}{
%%     ~~Describe \code{Xijh_i} here~~
}
  \item{pki}{
%%     ~~Describe \code{pki} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{methode}{
%%     ~~Describe \code{methode} here~~
}
  \item{e}{
%%     ~~Describe \code{e} here~~
}
  \item{nbr_iteration}{
%%     ~~Describe \code{nbr_iteration} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
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
function (Xquali, alpha_ki = NULL, Xijh_i = NULL, pki = NULL, 
    Y = NULL, methode = "init", e = 0.1, nbr_iteration = NULL, 
    k) 
{
    if (methode == "kmodes") {
        Y <- kmodes(Xquali, k)$cluster
    }
    else {
        Y <- init_class(iris, k)
    }
    if (is.null(alpha_ki)) {
        alpha_ki <- alphai(Y, Xquali)
    }
    if (is.null(Xijh_i)) {
        Xijh_i <- xijh(Xquali)
    }
    if (is.null(pki)) {
        pki <- as.numeric(table(Y)/length(Y))
    }
    Pki <- pki
    Alphaki <- alpha_ki
    p = length(colnames(Xquali))
    n <- nrow(Xquali)
    fkx <- fk_quali(Y, Xijh_i, alpha_ki, p = p, n = n)
    f <- fx(Pki, fkx, Y)
    tkx <- tk(fkx, Pki)
    li <- l(f)
    Qi <- Qt(tkx, fkx, Pki)
    BIC <- Bic(f, Y, Xquali = Xquali)
    ICL <- BIC - sum(tkx * log(tkx))
    i <- 1
    teta <- list()
    K <- k
    repeat {
        Y <- apply(tkx, MARGIN = 1, which.max)
        Nk <- nk(tkx)
        Pk <- pk(Nk, Xquali)
        Alphak <- alphak(tkx, Xijh_i, Y, Xquali)
        fkx <- fk_quali(Y, Xijh_i, Alphak, p = p, n = n)
        f <- fx(Pk, fkx, Y)
        lf <- l(f)
        tkx <- tk(fkx, Pk)
        Qf <- Qt(tkx, fkx, Pk)
        BIC <- Bic(f, Y, Xquali = Xquali)
        ICL <- BIC - sum(tkx * log(tkx))
        teta[[i]] <- list(iteration = i, Proba = Pk, Alphak = Alphak, 
            log_like = lf, Q = Qf, Bic = BIC, ICL = ICL, proba_individu = fkx, 
            part_MAP = tkx)
        if (!is.null(nbr_iteration)) {
            if (i >= nbr_iteration) 
                break
        }
        if (abs(li - lf) < e) {
            break
        }
        else {
            Pki <- Pk
            Alphaki <- Alphak
            li <- lf
            i <- i + 1
        }
    }
    return(teta)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
