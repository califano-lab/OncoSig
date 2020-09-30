\name{listToMatrix}
\alias{listToMatrix}
\title{
Convert Adjancency list to adjancency matrix
}
\description{
Convert Adjancency list (as a dataframe) to adjancency matrix. The first column will become the columns
the second column the rows, and the third column the weights of the edges. 
}

\usage{

listToMatrix(df)
}
\arguments{
  \item{df}{

}
}
\details{
}
\value{
returns a adjacency matrix as a matrix object
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
Network_location="Input_data_files/LUAD/original_network.txt"
Network=read.delim(Network_location,header=F)
Network$V1=as.character(Network$V1)
Network$V2=as.character(Network$V2)
Network$V3=as.numeric(Network$V3)
Network=as.matrix(Network)
Gold_Standard_location="Input_data_files/LUAD/10_oncogene_pathways/KRAS/total.txt"
Gold_Standard=read.delim(Gold_Standard_location,header=F)
Gold_Standard$V1=as.character(Gold_Standard$V1)

Network[,3]=as.numeric(Network[,3])
Network_matrix=listToMatrix(Network)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
