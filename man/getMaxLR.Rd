\name{getMaxLR}
\alias{getMaxLR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get maximum LR for correlated features
}
\description{
 This is a helper function for OncosigNB;Given the specified columns as a dataframe, return the maximum LR for each case. Note that these columns should be highly correlated.
}
\usage{
getMaxLR(df_1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df_1}{
%%     ~~Describe \code{df_1} here~~
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
function (df_1) 
{
    the_max_results = lapply(1:nrow(df_1), function(x) {
        max(df_1[x, ])
    })
    the_max_results = unlist(the_max_results)
    return(the_max_results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
