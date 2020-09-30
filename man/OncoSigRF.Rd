\name{OncoSigRF}
\alias{OncoSigRF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create and Generate Predictions Using the OncoSig Random Forest Classifier
}
\description{
This function creates the OncoSig Random Forest classifer, and returns predictions generated using Monte-Carlo cross validation. Optionally, the forests generated may be saved, which can be used to generate further predictions.
}
\usage{
OncoSigRF(Network_matrix_df, Gold_Standard_in_Network_names, Fraction_Gold_sample = NULL, ntrees = NULL, max_iterations = NULL, balance = NULL, to_save = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Network_matrix_df}{
  A network feature matrix, where the rownames are the points (e..g protein or gene names) to build the classifier on, and the columns are features.
  }
  \item{Gold_Standard_in_Network_names}{
  A list of proteins/genes in the gold standard. 
}
  \item{Fraction_Gold_sample}{
  Fraction of gold standard to sample in each Monte-Carlo run. Default:.5
}
  \item{ntrees}{
  Number of trees to create in each Random Forest (default:50)
}
  \item{max_iterations}{
  Number of iterations of Monte-Carlo samplings to run (i.e. number of forests ot create) (Default:20)
}
  \item{balance}{
  Ratio of  proteins not in the gold sample to sample in each run.Balance of 1 corresponds to an equal number of proteins in the gold sample and not. (Default:1) 
}
  \item{to_save}{
  Whether to save the forests created (in a file called "All_forests.R"). This argument must be set to 1 if unsupervised OncoSig is to be used.
}
}
\details{

}
\value{
Returns a dataframe corresponding to predictions from the Monte-Carlo cross-validation. Higher scores in the first column correspond to higher confidence predictions to be part of the oncogene centric map
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

}
\examples{
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
library (randomForest)
#Process the network
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


#Convert Matrix to Dataframe for future steps
Network_matrix_df=as.data.frame(Network_matrix)
#Remove Members of Gold Standard Not in the Network:
Gold_Standard_in_Network_names=intersect(rownames(Network_matrix_df),Gold_Standard$V1)
Negative_Set_names=setdiff(rownames(Network_matrix_df),Gold_Standard_in_Network_names)
#Create Negative_standard
#save(Network_matrix,file="Network_matrix.Rda")
remove(Network_matrix)

Query_output_results=OncoSigRF(Network_matrix_df,Gold_Standard_in_Network_names,max_iterations = 5)
Query_output_results_scores=Query_output_results[[1]]
}
\keyword{ Random Forest }% use one of  RShowDoc("KEYWORDS")
\keyword{ OncoSig }% __ONLY ONE__ keyword per line
\keyword{ Monte Carlo }% __ONLY ONE__ keyword per line
