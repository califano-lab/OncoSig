\name{getFinalLR}
\alias{getFinalLR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multipy LRs to get final LR
}
\description{
This is a helper function for OncoSigNB(). This function gets the final LR given the input  dataframe. The first two columns are the name and the response variable, columns 3 and up feature values.
}
\usage{
getFinalLR(df_1)
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
    to_return = lapply(1:nrow(df_1), function(x) {
        prod(df_1[x, 3:ncol(df_1)])
    })
    return(to_return)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
