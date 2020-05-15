\name{replaceBinswithLR}
\alias{replaceBinswithLR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Replace bins with LRs
}
\description{
This is a helper function for OnocoSigNB. Given the bin info, and bined data, replace each bin with the corresponding Likelihood ratio.
}
\usage{
replaceBinswithLR(bin_vector, the_bin_info)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bin_vector}{
%%     ~~Describe \code{bin_vector} here~~
}
  \item{the_bin_info}{
%%     ~~Describe \code{the_bin_info} here~~
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
function (bin_vector, the_bin_info) 
{
    new_bin_vector = lapply(1:len(bin_vector), function(x) {
        the_bin = as.character(bin_vector[x])
        bin_value = the_bin_info[the_bin]
        bin_value
    })
    new_bin_vector = unlist(new_bin_vector)
    return(new_bin_vector)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
