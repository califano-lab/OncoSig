\name{runNaiveBayesClassifier}
\alias{runNaiveBayesClassifier}
\title{
Run the OncoSig Naive Bayes Classifier as used in the acoompanying paper.
}
\description{
}
\usage{
runNaiveBayesClassifier()
}
\details{

}
\value{
}
\references{
}
\author{
}
\note{
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
runNaiveBayesClassifier()
## The function is currently defined as
function () 
{
    df_1 = read.delim("Input_data_files/Naive_Bayes_evidences_set_1.txt", 
        header = TRUE)
    df_2 = read.delim("Input_data_files/Naive_Bayes_evidences_set_2.txt", 
        header = TRUE)
    the_bins = list(c(0, 40, 200, 1200), c(0, 0.1), c(-2, -0.15, 
        -0.02, 0.0925), c(1, 2, 6), c(0, 0.25), c(1, 3, 20), 
        c(1, 4, 20), c(1, 4, 20), c(0, 1e-04, 0.9999), c(0, 0.01, 
            0.05))
    correlated_features = grep("MS_", colnames(df_1), value = TRUE)
    message("Calculating LR_posterior for fold two holdout set\n")
    the_results_set_1 = OncoSigNB(df_1, df_2, the_bins, correlated_features)
    message("Calculating LR_posterior for fold one holdout set\n")
    the_results_set_2 = OncoSigNB(df_2, df_1, the_bins, correlated_features)
    the_results_set_2_rank = cbind(the_results_set_2, rank(-the_results_set_2))
    the_results_set_1_rank = cbind(the_results_set_1, rank(-the_results_set_1))
    temp = rbind(the_results_set_1_rank, the_results_set_2_rank)
    temp = as.data.frame(temp)
    colnames(temp) = c("LR_post", "Rank")
    cross_validated_predictions = temp[order(temp$Rank), ]
    return(cross_validated_predictions)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
