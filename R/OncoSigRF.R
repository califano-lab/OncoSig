
OncoSigRF <- function(Network_matrix_df,Gold_Standard_in_Network_names, Fraction_Gold_sample=NULL, ntrees=NULL, max_iterations=
NULL,balance=NULL,to_save=NULL){
        #This function runs the Random Forest Learner using the Network Provided and the gold standard provided

        if(is.null(Fraction_Gold_sample)) {Fraction_Gold_sample=.5}
        if(is.null(ntrees)) {ntrees=50}
        if(is.null(max_iterations)) {max_iterations=20}
        if(is.null(balance)) {balance=1}
        if(is.null(to_save)) {to_save=0}
        message("Running OncoSig")
        message("Fraction of Gold Standard to train on for each Random Forest: ", Fraction_Gold_sample,sep="")
        message("Number of Trees per Iteration: ", ntrees,sep="")
        message("Number of Iterations: ", max_iterations,sep="")
        message("Balance: ", balance,sep="")
        #Sample gold_standard, user can Change this
        #ntrees=50 #Give each Random Fores 50 Trees, user can change this
        #max_iterations=100 #How many Iterations to do
        #Number to sample from the sets:
        Num_to_sample=floor(Fraction_Gold_sample*length(Gold_Standard_in_Network_names))
        #balance=3 #Change this depending on whether you want a balanced classifier or not, 1 means a balanced classifier, this will create more errors overall
        message("Number of positive results to sample: ", Num_to_sample,sep="")

        Num_to_sample_negative=Num_to_sample*balance;
        message("Number of negative results to sample: ", Num_to_sample_negative,sep="")
        QueryResults_scores=data.frame(row.names=rownames(Network_matrix_df))
        importance_df=data.frame()
        all_forests=list()
        for (i in 1:max_iterations){
                Gold_sample=sample(Gold_Standard_in_Network_names,Num_to_sample)
                Negative_sample=sample(Negative_Set_names,Num_to_sample_negative)
                label_vector=c(rep(1,Num_to_sample),rep(0,Num_to_sample_negative))
                label_vector=as.factor(label_vector)

                Not_in_Gold_or_Negative_Sample=setdiff(rownames(Network_matrix_df),c(Gold_sample,Negative_sample))
                #You can Cadd dotrace=TRUE if you want to see the trace of the random Forests
                message("Performing Random Forest",sep="")
                #Testing with fast by only passing it part of the matrix in the first place
                #my_col_sample=sample(colnames(Network_matrix_df),3000)


                #result=randomForest(Network_matrix_df[c(Gold_sample,Negative_sample),],label_vector,ntree = ntrees,importance=TRUE,do.trace=FALSE)
                #set mtry

                mtry=floor(ncol(Network_matrix_df)**.5)
                message("mtry equals ",mtry)
                result=randomForest(Network_matrix_df[c(Gold_sample,Negative_sample),],label_vector,ntree = ntrees,importance=TRUE,do.trace=TRUE,mtry=mtry)
                if (to_save==1){
                        all_forests[[i]]=result
                }
                Query_results=predict(result,type="prob",newdata=Network_matrix_df[Not_in_Gold_or_Negative_Sample,])
                QueryResults_scores[Not_in_Gold_or_Negative_Sample,i]=Query_results[Not_in_Gold_or_Negative_Sample,2]

                #Get the Importance, using mean decrease accuracy
                importance=as.data.frame(result$importance);importance=importance[order(importance$MeanDecreaseAccuracy,decreasing=T),]
                importance_vector=importance$MeanDecreaseAccuracy
                names(importance_vector)=rownames(importance)
                importance_df[names(importance_vector),i]=importance_vector
                #If the matrix is very large and you cannot query the results all at once, do it in chunks

                #How Converged are we if i>1
                if (i>2){
                        #QueryResults_scores_Complete_cases=QueryResults_scores[complete.cases(QueryResults_scores),]
                        #old_ones=rowMeans(QueryResults_scores_Complete_cases[,c(1:c(i-1))])
                        #old_plus_new=rowMeans(QueryResults_scores_Complete_cases[,c(1:c(i))])
                       #old_plus_new=rowMeans(QueryResults_scores_Complete_cases[,c(1:c(i))])
                        #Turn off warning for correlaiton, otherwise it spits back tie-related errors
                        #options(warn=-1)
                        #Correlation=cor.test(old_ones,old_plus_new,method="spearman")$estimate
                        #options(warn=0)
                        message("At iteration ", i,sep="")

                }
                Gold_sample_old=Gold_sample
                Negative_sample_old=Negative_sample
        }
        QueryResults_scores_average=as.data.frame(rowMeans(QueryResults_scores,na.rm=TRUE))
        colnames(QueryResults_scores_average)=c("Score")
        QueryResults_scores_average=as.data.frame(QueryResults_scores_average)
        if (to_save==1){
                        save(all_forests,file="All_forests.r")
                }
        return(list(QueryResults_scores_average,QueryResults_scores, importance_df))

}
