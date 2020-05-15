#For a ROC Curve (blah), given a False positive threshold (num) and the number of positives (num_pos), tell me on the that ROC curve the correspondoing true positive rate
#example: Roc_FPR(.01,perf_ELRON_no_interactions_KRB.R,250)
#This function reports back the FPR threshold used (which will be close to the threshold you input), the Trupe postive rate, the number of true positive found at this threshold and the number of false pistives at this threshold
Roc_FPR <- function(num,blah,num_pos)
{
	closest=1
	guess=abs(blah@x.values[[1]][1] -num)
	for (i in 2:length(blah@x.values[[1]])){
	guess2=abs(blah@x.values[[1]][[i]]-num)
		if (guess2<guess){
			guess=guess2
			closest=i
		}

	}
	#print(closest)
	#print(guess)
	print(blah@x.values[[1]][closest])
	print(blah@y.values[[1]][closest])


	print(blah@y.values[[1]][closest]*num_pos)
	print(blah@x.values[[1]][closest]*(19548-num_pos))
}
#Given Roc Curve (from ROCR, get the AUC)

Roc_AUC <- function(roc_object)
{
	library(sfsmisc)
	out=integrate.xy(roc_object@x.values[[1]],roc_object@y.values[[1]])
	return(out)
}

#Remove NA's from a dataframe row by row, and replace with mean of Row. 
replaceNARowMean <- function (df)
{
	for (i in 1:nrow(df)){
		my_mean=mean(as.numeric(df[i,]),na.rm=TRUE)
     	df[i,is.na(df[i,])] <- my_mean
     	message(i)
	}
	return(df)
}

#Plot many ROC Curves given the ROC curves, the colors to plot with, and the xlim
plotROC_curves <- function (ROC_curves,colors,xlim)
{
	plot(ROC_curves[1],col=colors[1],xlim=xlim)
	for (i in 2:length(ROC_curves))
	{
		par(new=T)
		plot(ROC_curves[i],col=colors[i],xlim=c(0,1))
	}
}

#An alias for View for Rstudio, since less is in bash
less <- function (blah)
{
	View(blah)
}

#Write a description of the figure in the margins in a small font
description <-function (blah)
{
	mtext(blah,cex=.3,outer =TRUE)
}
#alias for bash cd
cd <- function (blah)
{
	setwd(blah)
	
}
#alias for bash pwd
pwd <- function ()
{
	getwd()
	
}
#given a vector, add some random noise to each number. first vector is the number, second number is maximum number of noise to add (either up or down)
jiggle <- function (vector, noise)
{
	noise_1=runif(length(vector,min=-noise,max=noise))
	return(vector+noise_1)


}
#Integrate pvalues using fisher's method, given a list of pvalue as a list. Make sure you exclude NAs yourelf first:
fisherIntegration  <- function (vector){
	my_length=length(vector)
	deg_free=my_length*2
	y=-2*sum(log(vector))
	p.val <- 1-pchisq(y, df = deg_free);
	p.val=as.numeric(p.val);
	return(p.val)
}
#Integrate pvalues using stouffer's method, given a list of pvalue as a list. Make sure you exclude NAs yourelf first:
#Stouffers method takes a z-score, note the raw pvalue, so remember to convert to z-score first
stoufferIntegration <- function (vector,lower.tail=NULL){
	if(is.null(lower.tail)) {lower.tail=TRUE}
	my_length=length(vector)
	#zscores=scale(vector)
	p.val=sum(vector)/sqrt(my_length)
	p.val=pnorm(p.val,lower.tail=lower.tail)
	return(p.val)
}
#list_1 is a list, num_1 is how many fold validation you want to do (e.g. 10 fold validation)
#Returns 2 lists of lists, where the first one is a lists of all the training sets and the second one is the list of all the testing sets. These lists have the same length obviously
#A particular training set list can be retrived as   as: training_and_testing_folds[1][1][[1]][[1]]
#A particular testing set can be retrieved as: training_and_testing_folds[2][1][[1]][[1]]
kfolds <-function (list_1,num_1){
	list_1=sample(list_1)
	training_list_of_lists=list()
	testing_lists_of_lists=list()
	folds <- cut(seq(1,length(list_1)),breaks=num_1,labels=FALSE)
	folds_indices=unique(sort(folds))
	training_matrix=combn(folds_indices,num_1-1)
	for (i in 1:ncol(training_matrix)){
		training_folds=training_matrix[,i]
		testing_fold=setdiff(folds_indices,training_matrix[,i])
		training_folds_2=list_1[folds %in% training_folds]
		testing_fold_2=list_1[folds %in% testing_fold]
		training_list_of_lists[[i]]=training_folds_2
		testing_lists_of_lists[[i]]=testing_fold_2

	}

	return(list(training_list_of_lists,testing_lists_of_lists))

}

#plot to a plot folder as pdf
myPlot <- function (to_plot,file_name){
	file.pdf=paste("/ifs/home/c2b2/ac_lab/jb3401/Plots",file_name,sep="")
	#file.pdf=paste("~/Work_Files/Post-Doc/Plots/",file_name,sep="")
	pdf(file.pdf,width=8,height=6,paper='special') 
	plot(to_plot)
}

#create an empyt dataframe with the row and column names passed (as lists) to the function
namedrowscolsDataframe <- function (rownames,colnames){
	df_1=data.frame(matrix(NA, nrow = length(rownames), ncol = length(colnames)))
	rownames(df_1)=rownames
	colnames(df_1)=colnames
	return (df_1)

}
#this tries to evaluate the expression and returns it, but if it does not work, returns the alternative value. It is simply a wrapper for trycatch. This is similair to python's try except
#e.g.
#the_value=tryExcept(tissue_1[[1]][i,3],0)
#This will return the value of  tissue_1[[1]][i,3], unless it does not exists, in which case it will teturn zeron
tryExcept <- function(expr,alternative){
	tryCatch(expr,error=function(e) alternative)
}
#reorder a data frame based on a list. The returned dataframe is the same, except the rows are reordered. Note that this can screw up if your dataframe has only one column
reorderDataframe<-function(df_1,list_1){
	new_df=df_1[list_1,]
	return(new_df)
}


#python-like aliases:
len <- length
flattenList <- unlist
shape <- dim
index <-rownames
pwd <- getwd
intersection <-intersect
sorted <- sort
setDiff <- setdiff

