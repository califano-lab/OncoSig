\name{OncoSigNB}
\alias{OncoSigNB}
\title{
  Run the OncoSig Naive Bayes Classifier
}
\description{
  This function runs the OncoSign Naive Bayes Classifier, utilizing the user provided binning parameters. Optionally, allows the user to provide a list of features that are statistically dependent (and thus violate the assumption of Naive Bayes). The output of this function is a dataframe of predictions in the testing set whith scores based on training the classifier on the training set with corresponding likelihood ratios. Higher scores correspond to higher confidence predictions to be part of the oncogene-centric map.
}
\usage{
OncoSigNB(training_set, testing_set, the_bins, correlated_features)
}
\arguments{
  \item{training_set}{
      \code{a dataframe containing the training set}
}
  \item{testing_set}{
      \code{a dataframe containing the testing set}
}
  \item{the_bins}{
      \code{a list of list of the binning parameters. This list of list must be in the same order of the features/columns in the training and testing dataframes} 
}
  \item{correlated_features}{
      \code{a list of correlated features that are statistically dependent. Pass empty list if none}
}
}
\details{
In both the training and testing set, the first column should be a unique string identifying the datapoint (e.g. a protein id), and the second column is the label (0 or 1).
}
\value{
   returns a dataframe that is the predictions of classifier on the testing set

}
\references{
}
\author{

}
\note{


}

\seealso{

}
\examples{
#set bins
df_1=read.delim("Input_data_files/Naive_Bayes_evidences_set_1.txt",header=TRUE)
df_2=read.delim("Input_data_files/Naive_Bayes_evidences_set_2.txt",header=TRUE)
the_bins=list(c(0,40,200,1200),c(0,.1),c(-2,-0.15,-0.02,0.0925),c(1,2,6),c(0,0.25),c(1,3,20),c(1,4,20),c(1,4,20),c(0,0.0001,0.9999),c(0,0.01,0.05))
#specify correlated features
predictions=OncoSigNB(training_set = df_1,testing_set = df_2,the_bins=the_bins,correlated_features =list(correlated_features))
}
\keyword{ ~Naive Bayes }
