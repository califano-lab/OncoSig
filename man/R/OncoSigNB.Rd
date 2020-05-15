un the OncoSig Naive Bayes Classifier
\title{OncoSig Naive Bayes}
\description{
  Run the OncoSig Naive Bayes Classifier
  \code{save}.
}
\usage{
load(file, envir = parent.frame())
}
\arguments{
  \item{training_set}{a dataframe containing the training set}
  \item{testing_set}{a dataframe of the testing set}
  \item{the_bins}{a list of list of the binning parameters}
  \item{correlated_features}{a list of correlated features that are statistically dependent. Pass empty list if none}
}
\seealso{
  \code{\link{save}}.
}
\examples{
## set bins;get correlated features
df_1=read.delim("Input_data_files/Naive_Bayes_evidences_set_1.txt",header=TRUE)
df_2=read.delim("Input_data_files/Naive_Bayes_evidences_set_2.txt",header=TRUE)
the_bins=list(c(0,40,200,1200),c(0,.1),c(-2,-0.15,-0.02,0.0925),c(1,2,6),c(0,0.25),c(1,3,20),c(1,4,20),c(1,4,20),c(0,0.0001,0.9999),c(0,0.01,0.05))
predictions=OncoSigNB(training_set = df_1,testing_set = df_2,the_bins=the_bins,correlated_features =list())

\keyword{file}
