\name{Em_quanti}
\alias{Em_quanti}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do de clustering on guaussian futures. ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Em_quanti(Xquanti, meani = NULL, sigmai = NULL, pki = NULL, Y = NULL, methode = "init", e = 0.1, nbr_iteration = NULL, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xquanti}{
%%     ~~Describe \code{Xquanti} here~~
}
  \item{meani}{
%%     ~~Describe \code{meani} here~~
}
  \item{sigmai}{
%%     ~~Describe \code{sigmai} here~~
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
function (Xquanti, meani = NULL, sigmai = NULL, pki = NULL, Y = NULL, 
    methode = "init", e = 0.1, nbr_iteration = NULL, k) 
{
    if (methode == "kmeans") {
        Y <- kmeans(Xquanti, k)$cluster
    }
    else {
        Y <- init_class(iris, k)
    }
    if (is.null(pki)) {
        pki <- as.numeric(table(Y)/length(Y))
    }
    if (is.null(meani)) {
        meani <- aggregate(Xquanti, by = list(Y), FUN = mean)
    }
    if (is.null(sigmai)) {
        sigmai <- sigma_init(Xquanti, Y)
        Test <- TRUE
        for (i in 1:length(sigmai)) {
            if ((det(as.matrix(sigmai[[i]]))^(0.5)) < 0.001 || 
                is.na(det(as.matrix(sigmai[[2]]))^0.5)) 
                Test <- FALSE
        }
        if (!Test) {
            res <- PCA(Xquanti, scale.unit = TRUE, ncp = Inf, 
                graph = FALSE)
            X <- (res$ind[["coord"]])
            sigmai <- sigma_init(as.data.frame(X), Y)
            o <- 0
            vars <- colnames(Xquanti)
            VAR <- NULL
            XQQ <- NULL
            for (j in 1:ncol(Xquanti)) {
                t <- TRUE
                for (i in 1:k) {
                  if (sum(abs(sigmai[[i]][[j]])) < 0.1) {
                    o <- o + 1
                    t <- FALSE
                  }
                }
                if (t) {
                  XQQ <- cbind(XQQ, X[, j])
                }
            }
            for (j in 2:ncol(XQQ)) {
                sigmaii <- sigma_init(as.data.frame(XQQ[, 1:j]), 
                  Y)
                t <- TRUE
                for (i in 1:k) {
                  if (!is.na(det(as.matrix(sigmaii[[k]]))^(0.5))) {
                    if ((det(as.matrix(sigmaii[[k]]))^(0.5)) < 
                      0.1) {
                      t <- FALSE
                    }
                  }
                  else {
                    t <- FALSE
                  }
                  if (t) {
                    XX <- XQQ[, 1:j]
                  }
                }
                if (!t) {
                  break
                }
            }
            sigmaiii <- sigma_init(as.data.frame(XX), Y)
            sigmai <- sigma_init(as.data.frame(XQQ), Y)
            sigmai <- sigmaiii
            meani <- aggregate(XX, by = list(Y), FUN = mean)
            if (is.null(pki)) {
                pki <- as.numeric(table(Y)/length(Y))
            }
            Meanki <- meani
            Sigki <- sigmai
            Pki <- pki
            fkx <- fk_quanti2(XX, meani, sigmai, Y)
            f <- fx(Pki, fkx, Y)
            tkx <- tk(fkx, Pki)
            li <- l(f)
            Qi <- Qt(tkx, fkx, Pki)
            BIC <- Bic(f, Y, Xquanti = XX)
            ICL <- BIC - sum((matrix(unlist(fkx), ncol = k, byrow = TRUE)) * 
                log(matrix(unlist(fkx), ncol = k, byrow = TRUE)))
            i <- 1
            teta <- list()
            teta[[1]] <- list(iteration = 1, Mean = meani, Sigma = Sigki, 
                Proba = pki, log_like = li, Q = Qi, Bic = BIC, 
                Icl = ICL, proba_individu = matrix(unlist(fkx), 
                  ncol = k, byrow = TRUE), part_MAP = tkx)
            x11()
            repeat {
                pairs.panels(XX, gap = 0, bg = Y, pch = 21)
                Y <- apply(tkx, MARGIN = 1, which.max)
                Nk <- nk(tkx)
                Meank <- meank(Nk, tkx, XX, Y)
                Sigk <- sigmak(Nk, tkx, XX, Meanki)
                Pk <- pk(Nk, XX)
                fkx <- fk_quanti2(XX, Meank, Sigk, Y)
                f <- fx(Pk, fkx, Y)
                lf <- l(f)
                tkx <- tk(fkx, Pk)
                Qf <- Qt(tkx, fkx, Pk)
                BIC <- Bic(f, Y, Xquanti = XX)
                ICL <- BIC - sum(tkx * log1p(tkx))
                teta[[i + 1]] <- list(iteration = i + 1, Mean = Meank, 
                  Sigma = Sigk, Proba = Pk, log_like = lf, Q = Qf, 
                  Bic = BIC, Icl = ICL, proba_individu = matrix(unlist(fkx), 
                    ncol = k, byrow = TRUE), part_MAP = tkx)
                if (!is.null(nbr_iteration)) {
                  if (i >= nbr_iteration) 
                    break
                }
                if (abs(li - lf) < e) {
                  break
                }
                else {
                  Meanki <- as.data.frame(Meank)
                  Sigki <- Sigk
                  Pki <- Pk
                  Alphaki <- Alphak
                  li <- lf
                  i <- i + 1
                }
            }
        }
        else {
            Meanki <- meani
            Sigki <- sigmai
            Pki <- pki
            fkx <- fk_quanti(Xquanti, meani, sigmai, Y)
            f <- fx(Pki, fkx, Y)
            tkx <- tk(fkx, Pki)
            li <- l(f)
            Qi <- Qt(tkx, fkx, Pki)
            BIC <- Bic(f, Y, Xquanti = Xquanti)
            ICL <- BIC - sum((matrix(unlist(fkx), ncol = k, byrow = TRUE)) * 
                log(matrix(unlist(fkx), ncol = k, byrow = TRUE)))
            i <- 1
            teta <- list()
            repeat {
                pairs.panels(Xquanti, gap = 0, bg = Y, pch = 21)
                Y <- apply(tkx, MARGIN = 1, which.max)
                Nk <- nk(tkx)
                Meank <- meank(Nk, tkx, Xquanti, Y)
                Sigk <- sigmak(Nk, tkx, Xquanti, Meanki)
                Pk <- pk(Nk, Xquanti)
                fkx <- fk_quanti(Xquanti, Meank, Sigk, Y)
                f <- fx(Pk, fkx, Y)
                lf <- l(f)
                tkx <- tk(fkx, Pk)
                Qf <- Qt(tkx, fkx, Pk)
                BIC <- Bic(f, Y, Xquanti = Xquanti)
                ICL <- BIC - sum(tkx * log1p(tkx))
                teta[[i]] <- list(iteration = i, Mean = Meank, 
                  Sigma = Sigk, Proba = Pk, log_like = lf, Q = Qf, 
                  Bic = BIC, Icl = ICL, proba_individu = matrix(unlist(fkx), 
                    ncol = k, byrow = TRUE), part_MAP = tkx)
                if (!is.null(nbr_iteration)) {
                  if (i >= nbr_iteration) 
                    break
                }
                if (abs(li - lf) < e) {
                  break
                }
                else {
                  Meanki <- Meank
                  Sigki <- Sigk
                  Pki <- Pk
                  li <- lf
                  i <- i + 1
                }
            }
        }
    }
    else {
        Test <- TRUE
        for (i in 1:length(sigmai)) {
            if ((det(as.matrix(sigmai[[i]]))^(0.5)) < 0.001 || 
                is.na(det(as.matrix(sigmai[[2]]))^0.5)) 
                Test <- FALSE
        }
        if (!Test) {
            res <- PCA(Xquanti, scale.unit = TRUE, ncp = Inf, 
                graph = FALSE)
            X <- (res$ind[["coord"]])
            sigmai <- sigma_init(as.data.frame(X), Y)
            o <- 0
            vars <- colnames(Xquanti)
            VAR <- NULL
            XQQ <- NULL
            for (j in 1:ncol(Xquanti)) {
                t <- TRUE
                for (i in 1:k) {
                  if (sum(abs(sigmai[[i]][[j]])) < 0.1) {
                    o <- o + 1
                    t <- FALSE
                  }
                }
                if (t) {
                  XQQ <- cbind(XQQ, X[, j])
                }
            }
            for (j in 2:ncol(XQQ)) {
                sigmaii <- sigma_init(as.data.frame(XQQ[, 1:j]), 
                  Y)
                t <- TRUE
                for (i in 1:k) {
                  if (!is.na(det(as.matrix(sigmaii[[k]]))^(0.5))) {
                    if ((det(as.matrix(sigmaii[[k]]))^(0.5)) < 
                      0.1) {
                      t <- FALSE
                    }
                  }
                  else {
                    t <- FALSE
                  }
                  if (t) {
                    XX <- XQQ[, 1:j]
                  }
                }
                if (!t) {
                  break
                }
            }
            sigmaiii <- sigma_init(as.data.frame(XX), Y)
            sigmai <- sigma_init(as.data.frame(XQQ), Y)
            sigmai <- sigmaiii
            meani <- aggregate(XX, by = list(Y), FUN = mean)
            if (is.null(pki)) {
                pki <- as.numeric(table(Y)/length(Y))
            }
            Meanki <- meani
            Sigki <- sigmai
            Pki <- pki
            fkx <- fk_quanti2(XX, meani, sigmai, Y)
            f <- fx(Pki, fkx, Y)
            tkx <- tk(fkx, Pki)
            li <- l(f)
            Qi <- Qt(tkx, fkx, Pki)
            BIC <- Bic(f, Y, Xquanti = XX)
            ICL <- BIC - sum((matrix(unlist(fkx), ncol = k, byrow = TRUE)) * 
                log(matrix(unlist(fkx), ncol = k, byrow = TRUE)))
            i <- 1
            teta <- list()
            teta[[1]] <- list(iteration = 1, Mean = meani, Sigma = Sigki, 
                Proba = pki, log_like = li, Q = Qi, Bic = BIC, 
                Icl = ICL, proba_individu = matrix(unlist(fkx), 
                  ncol = k, byrow = TRUE), part_MAP = tkx)
            x11()
            repeat {
                pairs.panels(XX, gap = 0, bg = Y, pch = 21)
                Y <- apply(tkx, MARGIN = 1, which.max)
                Nk <- nk(tkx)
                Meank <- meank(Nk, tkx, XX, Y)
                Sigk <- sigmak(Nk, tkx, XX, Meanki)
                Pk <- pk(Nk, XX)
                fkx <- fk_quanti2(XX, Meank, Sigk, Y)
                f <- fx(Pk, fkx, Y)
                lf <- l(f)
                tkx <- tk(fkx, Pk)
                Qf <- Qt(tkx, fkx)
                BIC <- Bic(f, Y, Xquanti = XX)
                ICL <- BIC - sum(tkx * log1p(tkx))
                teta[[i + 1]] <- list(iteration = i + 1, Mean = Meank, 
                  Sigma = Sigk, Proba = Pk, log_like = lf, Q = Qf, 
                  Bic = BIC, Icl = ICL, proba_individu = matrix(unlist(fkx), 
                    ncol = k, byrow = TRUE), part_MAP = tkx)
                if (!is.null(nbr_iteration)) {
                  if (i >= nbr_iteration) 
                    break
                }
                if (abs(li - lf) < e) {
                  break
                }
                else {
                  Meanki <- as.data.frame(Meank)
                  Sigki <- Sigk
                  Pki <- Pk
                  Alphaki <- Alphak
                  li <- lf
                  i <- i + 1
                }
            }
        }
        else {
            Meanki <- meani
            Sigki <- sigmai
            Pki <- pki
            fkx <- fk_quanti(Xquanti, meani, sigmai, Y)
            f <- fx(Pki, fkx, Y)
            tkx <- tk(fkx, Pki)
            li <- l(f)
            Qi <- Qt(tkx, fkx, Pki)
            BIC <- Bic(f, Y, Xquanti = Xquanti)
            ICL <- BIC - sum((matrix(unlist(fkx), ncol = k, byrow = TRUE)) * 
                log(matrix(unlist(fkx), ncol = k, byrow = TRUE)))
            i <- 1
            teta <- list()
            repeat {
                pairs.panels(Xquanti, gap = 0, bg = Y, pch = 21)
                Y <- apply(tkx, MARGIN = 1, which.max)
                Nk <- nk(tkx)
                Meank <- meank(Nk, tkx, Xquanti, Y)
                Sigk <- sigmak(Nk, tkx, Xquanti, Meanki)
                Pk <- pk(Nk, Xquanti)
                fkx <- fk_quanti(Xquanti, Meank, Sigk, Y)
                f <- fx(Pk, fkx, Y)
                lf <- l(f)
                tkx <- tk(fkx, Pk)
                Qf <- Qt(tkx, fkx, Pk)
                BIC <- Bic(f, Y, Xquanti = Xquanti)
                ICL <- BIC - sum(tkx * log1p(tkx))
                teta[[i]] <- list(iteration = i, Mean = Meank, 
                  Sigma = Sigk, Proba = Pk, log_like = lf, Q = Qf, 
                  Bic = BIC, Icl = ICL, proba_individu = matrix(unlist(fkx), 
                    ncol = k, byrow = TRUE), part_MAP = tkx)
                if (!is.null(nbr_iteration)) {
                  if (i >= nbr_iteration) 
                    break
                }
                if (abs(li - lf) < e) {
                  break
                }
                else {
                  Meanki <- Meank
                  Sigki <- Sigk
                  Pki <- Pk
                  li <- lf
                  i <- i + 1
                }
            }
        }
    }
    return(teta)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
