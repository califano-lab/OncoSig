\name{NaiveBayesBin}
\alias{NaiveBayesBin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bins features according to bin parameters
}
\description{
This is a helper function for OncoSigNB() that takes in a dataframe and bins for each feature and returns, for each bin, the corresponding likelihood Ratios (
Where LR=p(1|bin)/(p(1))
}
\usage{
NaiveBayesBin(df_1, the_bins)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df_1}{
%%     ~~Describe \code{df_1} here~~
}
  \item{the_bins}{
%%     ~~Describe \code{the_bins} here~~
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
function (df_1, the_bins) 
{
    the_min = -1e+45
    the_max = 1e+43
    df_1[is.na(df_1)] <- -1e+23
    df_1_copy = df_1
    for (i in 3:ncol(df_1)) {
        df_1_copy[, i] = .bincode(df_1_copy[, i], c(the_min, 
            the_bins[i - 2][[1]], the_max), right = FALSE)
    }
    new_bin_info = list()
    for (i in 3:ncol(df_1)) {
        new_bin_info = append(new_bin_info, list(getLRsgivenBin_info(df_1_copy[, 
            i], the_bins[i - 2][[1]], df_1_copy[, 2])))
    }
    return(new_bin_info)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
