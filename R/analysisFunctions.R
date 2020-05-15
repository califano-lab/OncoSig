#This function invokes runs the Naive Bayes OncoSig classifier to replicate the results presented in
#the accompanying paper
source("R/rFunctions.R")

runNaiveBayesClassifier <- function(){
	#read in training and testing set
	df_1=read.delim("Input_data_files/Naive_Bayes_evidences_set_1.txt",header=TRUE)
	df_2=read.delim("Input_data_files/Naive_Bayes_evidences_set_2.txt",header=TRUE)
	
	#set binning parameters
	the_bins=list(c(0,40,200,1200),c(0,.1),c(-2,-0.15,-0.02,0.0925),c(1,2,6),c(0,0.25),c(1,3,20),c(1,4,20),c(1,4,20),c(0,0.0001,0.9999),c(0,0.01,0.05))
	correlated_features=grep("MS_",colnames(df_1),value = TRUE)
	
	#perform two fold cross validation
	message("Calculating LR_posterior for fold two holdout set\n")
	the_results_set_1=OncoSigNB(df_1,df_2,the_bins,correlated_features)
	message("Calculating LR_posterior for fold one holdout set\n")
	the_results_set_2=OncoSigNB(df_2,df_1,the_bins,correlated_features)

	#rank the results
	the_results_set_2_rank=cbind(the_results_set_2,rank(-the_results_set_2))
	the_results_set_1_rank=cbind(the_results_set_1,rank(-the_results_set_1))
	temp=rbind(the_results_set_1_rank,the_results_set_2_rank)
	temp=as.data.frame(temp)
	colnames(temp)=c("LR_post","Rank")
	cross_validated_predictions=temp[order(temp$Rank),]
	return(cross_validated_predictions)
	#function for performing Naive Bayes Classification  

}


#This function calls a script that generates a ROC Curve
generateROCcurve<- function(object_to_create,column_to_use,predictions_file,gold_standard,pdf_outfile){
	setwd("Output_files")
	cmd=paste("../scripts/quickROC.pl -s",object_to_create,"-c",column_to_use,predictions_file,gold_standard,pdf_outfile,sep=" ")
	system(cmd)
	setwd("..")
	#print(cmd)
		
}

#For a dataframe containing geneids and Log fold change values, this function finds the p-value of each Log Fold Change
#This loop assigns a p-value to each individual shRNA (note that there are multiple shRNas targeting each gene)
getPvalueofLogFC <- function (df_1,density_null){
	for (i in 1:nrow(df_1)) {
		number=df_1[i,2]
		Avg.pos <- number;
		xt <- diff(density_null$x[density_null$x < Avg.pos]);
		#integrate over the density
		yt <- rollmean(density_null$y[density_null$x < Avg.pos  ],2);
	 	pvalue=sum(xt*yt)
	 	df_1[i,3]=pvalue
	#print(i)
	}
	#Due to errors in integration rounding, some p-values may be above 1, set those to 1.
	above_1=which((df_1[,3]) > 1) 
	df_1[above_1,3]=1
	return(df_1)
}

#This function Integrates pvalues of the same genes using fisher integration:
#Set maximum pvalues to 1, there are some above 1 do to rounding errors in the integration. The input is a dataframe of genes and raw p-values
integratePvaluesbyGene <- function(df_1){
	gene_ids=unique(sort(df_1$Gene))
	Integrated_pvalues=data.frame(row.names = gene_ids)
	#Integrate the values using fisher integration
	for (i in rownames(Integrated_pvalues)) {
		nums=df_1[which(df_1$Gene==i),2]
		Integrated_pvalue=fisherIntegration(nums)
		Integrated_pvalues[i,1]=Integrated_pvalue
		#print(Integrated_pvalue)
	}
	Integrated_pvalues_2=Integrated_pvalues[order(Integrated_pvalues$V1),,drop=F]
	return(Integrated_pvalues_2)

}

generateROCcurve10OncogenePathways <- function(){
	system("scripts/generateROC_curves_OncosigRF.sh")
}

#this function gets all pairwise pearson correlations between two dataframe columns, and returns it as a vector

getPairwiseCordataframes <- function(df_1,df_2){
	to_return=list()
	for (i in colnames(df_1)) {
		for (j in colnames(df_2)){
			z=cor.test(df_1[,i],df_2[,j])
			#print(c(i,j,z))
			to_return=append(to_return,z$estimate)
		}
	}
	to_return=unlist(to_return)
}



