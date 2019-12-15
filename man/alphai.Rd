\name{alphai}
\alias{alphai}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
alphai(Y, Xquali)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
function (Y, Xquali) 
{
    Xquali_k <- as.list(unique(Y))
    if (is.null(dim(Xquali))) {
        j <- 1
        for (i in unique(Y)) {
            Xquali_k[[j]] <- data.frame(MAT = Xquali[Y == i])
            j <- j + 1
        }
        alpha_k <- as.list(unique(Y))
        for (k in 1:length(alpha_k)) {
            variables <- vector("list", 1)
            X <- as.numeric(table(Xquali_k[[k]])/nrow(Xquali_k[[k]]))
            variables[[1]] <- data.frame(proba = X, modalite = unlist(attributes(table(Xquali_k[[k]])/nrow(Xquali_k[[k]]))$dimnames))
            for (h in unique(Xquali)) {
                if (!(h \%in\% variables[[1]]$modalite)) {
                  mod_miss <- data.frame(proba = 0, modalite = h)
                  variables[[1]] <- rbind(variables[[1]], mod_miss)
                }
            }
            variables[[1]]$modalite <- as.character(variables[[1]]$modalite)
            variables[[1]] <- variables[[1]][order(variables[[1]]$modalite, 
                decreasing = F), ]
            alpha_k[[k]] <- variables
        }
    }
    else {
        j <- 1
        for (i in unique(Y)) {
            Xquali_k[[j]] <- data.frame(MAT = Xquali[Y == i, 
                ])
            j <- j + 1
        }
        alpha_k <- as.list(unique(Y))
        for (k in 1:length(alpha_k)) {
            variables <- as.list(colnames(Xquali))
            for (j in 1:length(variables)) {
                X <- as.numeric(table(Xquali_k[[k]][, j])/nrow(Xquali_k[[k]]))
                variables[[j]] <- data.frame(proba = X, Var = variables[[j]], 
                  modalite = unlist(attributes(table(Xquali_k[[k]][, 
                    j])/nrow(Xquali_k[[k]]))$dimnames))
                for (h in unique(Xquali[, j])) {
                  if (!(h \%in\% variables[[j]]$modalite)) {
                    mod_miss <- data.frame(proba = 0, Var = unique(variables[[j]]$Var), 
                      modalite = h)
                    variables[[j]] <- rbind(variables[[j]], mod_miss)
                  }
                }
                variables[[j]]$modalite <- as.character(variables[[j]]$modalite)
                variables[[j]] <- variables[[j]][order(variables[[j]]$modalite, 
                  decreasing = F), ]
            }
            alpha_k[[k]] <- variables
        }
    }
    return(alpha_k)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
