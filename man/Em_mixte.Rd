\name{Em_mixte}
\alias{Em_mixte}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do  function to do de clustering a Dataset containing  guaussian futures.and categorial futures ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Em_mixte(Xquanti, Xquali, meani = NULL, sigmai = NULL, alpha_ki = NULL, Xijh_i = NULL, pki = NULL, Y = NULL, methode = "init", e = 0.1, nbr_iteration = NULL, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xquanti}{
%%     ~~Describe \code{Xquanti} here~~
}
  \item{Xquali}{
%%     ~~Describe \code{Xquali} here~~
}
  \item{meani}{
%%     ~~Describe \code{meani} here~~
}
  \item{sigmai}{
%%     ~~Describe \code{sigmai} here~~
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
function (Xquanti, Xquali, meani = NULL, sigmai = NULL, alpha_ki = NULL, 
    Xijh_i = NULL, pki = NULL, Y = NULL, methode = "init", e = 0.1, 
    nbr_iteration = NULL, k) 
{
    if (methode == "kmeans") {
        Y <- kmeans(Xquanti, k)$cluster
    }
    else {
        if (methode == "kmodes") {
            Y <- kmodes(Xquali, k)$cluster
        }
        else {
            Y <- init_class(iris, k)
        }
    }
    if (is.null(meani)) {
        meani <- aggregate(Xquanti, by = list(Y), FUN = mean)
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
            Alphaki <- alpha_ki
            fkx <- fk2(XX, meani, sigmai, Y, Xijh_i, alpha_ki, 
                p = length(Xquali))
            f <- fx(Pki, fkx, Y)
            tkx <- tk(fkx, Pki)
            li <- l(f)
            Qi <- Qt(tkx, fkx, Pki)
            BIC <- Bic(f, Y, Xquanti = XX)
            ICL <- BIC - sum((matrix(unlist(fkx), ncol = k, byrow = TRUE)) * 
                log(matrix(unlist(fkx), ncol = k, byrow = TRUE)))
            i <- 1
            teta <- list()
            x11()
            repeat {
                pairs.panels(XX, gap = 0, bg = Y, pch = 21)
                Y <- apply(tkx, MARGIN = 1, which.max)
                Nk <- nk(tkx)
                Meank <- meank(Nk, tkx, XX, Y)
                Sigk <- sigmak(Nk, tkx, XX, Meanki)
                Alphak <- alphak(tkx, Xijh_i, Y, Xquali)
                Pk <- pk(Nk, XX)
                fkx <- fk2(XX, Meank, Sigk, Y, Xijh_i, Alphak, 
                  p = ncol(Xquali))
                f <- fx(Pk, fkx, Y)
                lf <- l(f)
                tkx <- tk(fkx, Pk)
                Qf <- Qt(tkx, fkx, Pk)
                BIC <- Bic(f, Y, Xquanti = XX)
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
            Alphaki <- alpha_ki
            fkx <- fk(Xquanti, meani, sigmai, Y, Xijh_i, alpha_ki, 
                p = length(Xquali))
            f <- fx(Pki, fkx, Y)
            tkx <- tk(fkx, Pki)
            li <- l(f)
            Qi <- Qt(tkx, fkx, Pki)
            BIC <- Bic(f, Y, Xquali = Xquali, Xquanti = Xquanti)
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
                Alphak <- alphak(tkx, Xijh_i, Y, Xquali)
                fkx <- fk(Xquanti, Meank, Sigk, Y, Xijh_i, Alphak, 
                  p = ncol(Xquali))
                f <- fx(Pk, fkx, Y)
                lf <- l(f)
                tkx <- tk(fkx, Pk)
                Qf <- Qt(tkx, fkx, Pk)
                BIC <- Bic(f, Y, Xquali = Xquali, Xquanti = Xquanti)
                ICL <- BIC - sum(tkx * log1p(tkx))
                teta[[i]] <- list(iteration = i, Mean = Meank, 
                  Sigma = Sigk, Proba = Pk, Alphak = Alphak, 
                  log_like = lf, Q = Qf, Bic = BIC, Icl = ICL, 
                  proba_individu = matrix(unlist(fkx), ncol = k, 
                    byrow = TRUE), part_MAP = tkx)
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
                  Alphaki <- Alphak
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
            Alphaki <- alpha_ki
            fkx <- fk2(XX, meani, sigmai, Y, Xijh_i, alpha_ki, 
                p = length(Xquali))
            f <- fx(Pki, fkx, Y)
            tkx <- tk(fkx, Pki)
            li <- l(f)
            Qi <- Qt(tkx, fkx, Pki)
            BIC <- Bic(f, Y, Xquanti = XX)
            ICL <- BIC - sum((matrix(unlist(fkx), ncol = k, byrow = TRUE)) * 
                log(matrix(unlist(fkx), ncol = k, byrow = TRUE)))
            i <- 1
            teta <- list()
            x11()
            repeat {
                pairs.panels(XX, gap = 0, bg = Y, pch = 21)
                Y <- apply(tkx, MARGIN = 1, which.max)
                Nk <- nk(tkx)
                Meank <- meank(Nk, tkx, XX, Y)
                Sigk <- sigmak(Nk, tkx, XX, Meanki)
                Alphak <- alphak(tkx, Xijh_i, Y, Xquali)
                Pk <- pk(Nk, XX)
                fkx <- fk2(XX, Meank, Sigk, Y, Xijh_i, Alphak, 
                  p = ncol(Xquali))
                f <- fx(Pk, fkx, Y)
                lf <- l(f)
                tkx <- tk(fkx, Pk)
                Qf <- Qt(tkx, fkx, Pk)
                BIC <- Bic(f, Y, Xquanti = XX)
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
            Alphaki <- alpha_ki
            fkx <- fk(Xquanti, meani, sigmai, Y, Xijh_i, alpha_ki, 
                p = length(Xquali))
            f <- fx(Pki, fkx, Y)
            tkx <- tk(fkx, Pki)
            li <- l(f)
            Qi <- Qt(tkx, fkx, Pki)
            BIC <- Bic(f, Y, Xquali = Xquali, Xquanti = Xquanti)
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
                Alphak <- alphak(tkx, Xijh_i, Y, Xquali)
                fkx <- fk(Xquanti, Meank, Sigk, Y, Xijh_i, Alphak, 
                  p = ncol(Xquali))
                f <- fx(Pk, fkx, Y)
                lf <- l(f)
                tkx <- tk(fkx, Pk)
                Qf <- Qt(tkx, fkx, Pk)
                BIC <- Bic(f, Y, Xquali = Xquali, Xquanti = Xquanti)
                ICL <- BIC - sum(tkx * log1p(tkx))
                teta[[i]] <- list(iteration = i, Mean = Meank, 
                  Sigma = Sigk, Proba = Pk, Alphak = Alphak, 
                  log_like = lf, Q = Qf, Bic = BIC, Icl = ICL, 
                  proba_individu = matrix(unlist(fkx), ncol = k, 
                    byrow = TRUE), part_MAP = tkx)
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
                  Alphaki <- Alphak
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
