\name{Plot_quali}
\alias{Plot_quali}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Plot_quali(Dataquali, Em_results, all = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Dataquali}{
%%     ~~Describe \code{Dataquali} here~~
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
function (Dataquali, Em_results, all = TRUE) 
{
    tkx <- Em_results[[length(Em_results)]][["part_MAP"]]
    Y <- apply(tkx, MARGIN = 1, which.max)
    if (all) {
        layout(matrix(1:ncol(Dataquali), ncol = floor(sqrt(ncol(Dataquali)))))
        names_vars <- colnames(Dataquali)
        for (i in 1:ncol(Dataquali)) {
            my_variable <- table(Dataquali[, i], Y)
            my_variable2 <- prop.table(my_variable, 1)
            barplot(my_variable2, beside = TRUE, legend = paste("cluster", 
                rownames(my_variable2)), ylab = "proportion du cluster ", 
                xlab = names_vars[i], col = 1:length(unique(Y)))
        }
    }
    else {
        names_vars <- colnames(Dataquali)
        for (i in 1:ncol(Dataquali)) {
            my_variable <- table(Y, Dataquali[, i])
            my_variable2 <- prop.table(my_variable, 1)
            barplot(my_variable2, beside = TRUE, legend = paste("cluster", 
                rownames(my_variable2)), ylab = "proportion du cluster ", 
                xlab = names_vars[i], col = 1:length(unique(Y)))
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
