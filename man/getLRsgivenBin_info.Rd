\name{getLRsgivenBin_info}
\alias{getLRsgivenBin_info}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return vector of Likelihood Ratios for each bin.
}
\description{
This is a helper function for OncoSigNB. Given a bin vector and the gold standard vector (i.e. the two vectors of the same length), return the Likelihood Ratio vector.
}
\usage{
getLRsgivenBin_info(bin_vector, the_bin, label_vector)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bin_vector}{
%%     ~~Describe \code{bin_vector} here~~
}
  \item{the_bin}{
%%     ~~Describe \code{the_bin} here~~
}
  \item{label_vector}{
%%     ~~Describe \code{label_vector} here~~
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
function (bin_vector, the_bin, label_vector) 
{
    the_bins_new = the_bin
    bin_vector_new = bin_vector
    prior = table(label_vector)[2]/table(label_vector)[1]
    bin_vector_2 = unique(sort(bin_vector))
    bin_vector_3 = rep(0, length(bin_vector_2))
    for (i in 1:length(bin_vector_2)) {
        the_num = bin_vector_2[i]
        ratio_1 = table(label_vector[bin_vector == i])[2]/table(label_vector[bin_vector == 
            i])[1]
        LR = ratio_1/prior
        bin_vector_3[i] = LR
    }
    names(bin_vector_3) = bin_vector_2
    return(bin_vector_3)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
