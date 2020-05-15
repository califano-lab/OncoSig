source("./R/rFunctions.R")
OncoSigUnsup <- function(Network_location, forest_location){
	load(forest_location,verbose=T)

	Network=read.delim(Network_location,header=F)
	Network$V1=as.character(Network$V1)
	Network$V2=as.character(Network$V2)
	Network$V3=as.numeric(Network$V3)
	Network=as.matrix(Network)

	#Convert to Matrix. Inputes missing values as 0, so make sure your scores range from  greater than zero to higher!
	Network[,3]=as.numeric(Network[,3])
	Network_matrix=listToMatrix(Network)
        
	result_matrix=matrix(nrow=nrow(Network_matrix),ncol=length(all_forests))
	rownames(result_matrix)=rownames(Network_matrix)
	for (i in 1:length(all_forests)){
        	Query_results=predict(all_forests[[i]],newdata = Network_matrix,type="prob")
        	result_matrix[,i]=Query_results[,2]
}
	the_means=rowMeans(result_matrix)
	the_means_df=as.data.frame(the_means)
	return(the_means_df)
}
