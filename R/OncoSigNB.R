OncoSigNB <- function (training_set,testing_set,the_bins,correlated_features){
	message("\tBinning features based on holdout set\n")
	the_bins_info=NaiveBayesBin(training_set,the_bins)
	testing_set=computeLRsgivenBins(testing_set,the_bins,the_bins_info)
	#Get the maximum of correlated features if they were passed
	if (len(correlated_features)>0){
		mass_spec_features=correlated_features
		the_mass_spec_features=testing_set[mass_spec_features]
		#print(the_mass_spec_features)
		message("\tCorrecting for correlated features \n")
		#message("test")
		ms_max=getMaxLR(the_mass_spec_features)

		testing_set[mass_spec_features] <- 1
		testing_set$MS_max=ms_max
	}
	message("\tExtracting predicted LR_posterior holdout set\n")
	the_results=unlist(getFinalLR(testing_set))
	names(the_results)=testing_set$V1

	return(the_results)
}
