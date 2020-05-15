source("R/rFunctions.R")
#These sets of functions allow the creation of a binary Naive Bayes classifer, given two sets of data. 
#1. The feature matrix
#2. The bin parameters. 

#1.Feature Matrix
#The first column in the feature matrix is the name of the instance (e.g. the name of a protein)
#The second column is the feature value (either a 1 or 0)
#All later columns are feature values, which can be either discrete or continous numeric values, or "NA"
#For example:
# 	V1	df_labels	feature_1 feature_2
#1	Q16539	1	6.76840e+02	1.000e+00	
#2	P78383	1	NA	NA	
#3	P30281	1	NA	1.000e+00

#2. Bin Parameters
#bin parameters are passed to the function as lists of lists, with one list per feature: For example
#the_bins=list(c(0,40,200,1200),c(0,.1)
#will bin the first feature into bins corresponding to: 0-40,40-200,200-1200,1200-Inf. "NA's" (i.e. no feature present) are given a seperate bin for each feature.
#NOTE: Due to the way that "NA"'s are imputed, features should not have any features less than -99999999999999999999999, or greater than 9999999999999999999999999999999999999999.
#NOTE: These functions assume that Naive Bayes is used as a binary classifier, where there are only two labels in the response vector:1 or 0.

#This function takes in a dataframe and a bins for each feature and returns, for each bin, the corresponding likelihood Ratios (LR)
#Where LR=p(1|bin)/(p(1))
NaiveBayesBin <-function (df_1,the_bins){
	#impute all NAs with a very small value temporarily.
	 the_min=-999999999999999999999999999999999999999999999
	 the_max=9999999999999999999999999999999999999999999
	 df_1[is.na(df_1)] <- -99999999999999999999999
	 df_1_copy=df_1

	 #go through the data frame with the assigned breaks and create the bins
	 #the first two columns must be the name of the entry and the labels; all future columns are features
	 #for each feature, find the proper bins, which is 2 less than the according column
	 for (i in 3:ncol(df_1)){
	 	df_1_copy[,i]=.bincode(df_1_copy[,i],c(the_min,the_bins[i-2][[1]],the_max),right=FALSE)
	 	#print(i)
	 }
	 new_bin_info=list()
	 for (i in 3:ncol(df_1)){
	 	new_bin_info=append(new_bin_info,list(getLRsgivenBin_info(df_1_copy[,i],the_bins[i-2][[1]],df_1_copy[,2])))

	 }
	return(new_bin_info)


}
#After training on a training set, this function computes LRs on a new testing set. Note that labels must be provided for the testing set as well.
computeLRsgivenBins <- function (df_1,the_bins,the_bins_info){
#impute all NAs with a very small value temporarily.
	 the_min=-999999999999999999999999999999999999999999999
	 the_max=9999999999999999999999999999999999999999999
	 df_1[is.na(df_1)] <- -99999999999999999999999
	 df_1_copy=df_1

	 #go through the data frame with the assigned breaks and create the bins
	 #the first two columns must be the name of the entry and the labels; all future columns are features
	 #for each feature, find the proper bins, which is 2 less than the according column
	 for (i in 3:ncol(df_1)){
	 	df_1_copy[,i]=.bincode(df_1_copy[,i],c(the_min,the_bins[i-2][[1]],the_max),right=FALSE)
	 	#print(i)
	 }
	 #replace bins with LR
	 for (i in 3:ncol(df_1)){
	 	df_1_copy[,i]=replaceBinswithLR(df_1_copy[,i],the_bins_info[i-2][[1]])
	 	#print(i)
	 }
	 return(df_1_copy)
}

#given a bin vector and the gold standard vector (i.e. the two vectors of the same length), return the Likelihood Ratio vector
getLRsgivenBin_info <- function (bin_vector,the_bin,label_vector){
	#get bins info
	the_bins_new=the_bin
	bin_vector_new=bin_vector
	prior=table(label_vector)[2]/table(label_vector)[1]
	bin_vector_2=unique(sort(bin_vector))
	bin_vector_3=rep(0,length(bin_vector_2))
	for (i in 1:length(bin_vector_2)) {
		the_num=bin_vector_2[i]
		ratio_1=table(label_vector[bin_vector==i])[2]/table(label_vector[bin_vector==i])[1]
		LR=ratio_1/prior
		bin_vector_3[i]=LR

	}
	names(bin_vector_3)=bin_vector_2
	#for (i in 1:length(bin_vector_new)) {
	#	the_bin_value=bin_vector_new[i]
	#	bin_vector_new[i]=bin_vector_3[the_bin_value]
	#}

	return(bin_vector_3)

}

#get the final LR given the dataframe. Columns starting at 3 are feature values
getFinalLR <- function(df_1){
	to_return=lapply(1:nrow(df_1),
		function(x){
			prod(df_1[x,3:ncol(df_1)])

		}

	
	)
	return (to_return)

}

#given the specified columns, return the maximum LR for each case
getMaxLR <- function(df_1){
	the_max_results=lapply(1:nrow(df_1),
		function(x){
			max(df_1[x,])
		}
		)
	the_max_results=unlist(the_max_results)
	return(the_max_results)
}
#Given the bin info, and bined data, replace each bin with the corresponding Likelihood ratio.
replaceBinswithLR <- function(bin_vector,the_bin_info){
	new_bin_vector=lapply(1:len(bin_vector),
		function(x){
			the_bin=as.character(bin_vector[x])
			bin_value=the_bin_info[the_bin]
			#print(x)
			bin_value


		}

		)
	new_bin_vector=unlist(new_bin_vector)
	return(new_bin_vector)

}
