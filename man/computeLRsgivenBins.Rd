\name{computeLRsgivenBins}
\alias{computeLRsgivenBins}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Likelihood Ratios (LR) for testing set
}
\description{
This is a helper function for OncoSigNB. After training on a training set, this function computes LRs on a new testing set. Note that labels must be provided for the testing set as well.
}
\usage{
computeLRsgivenBins(df_1, the_bins, the_bins_info)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df_1}{
%%     ~~Describe \code{df_1} here~~
}
  \item{the_bins}{
%%     ~~Describe \code{the_bins} here~~
}
  \item{the_bins_info}{
%%     ~~Describe \code{the_bins_info} here~~
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
function (df_1, the_bins, the_bins_info) 
{
    the_min = -1e+45
    the_max = 1e+43
    df_1[is.na(df_1)] <- -1e+23
    df_1_copy = df_1
    for (i in 3:ncol(df_1)) {
        df_1_copy[, i] = .bincode(df_1_copy[, i], c(the_min, 
            the_bins[i - 2][[1]], the_max), right = FALSE)
    }
    for (i in 3:ncol(df_1)) {
        df_1_copy[, i] = replaceBinswithLR(df_1_copy[, i], the_bins_info[i - 
            2][[1]])
    }
    return(df_1_copy)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
