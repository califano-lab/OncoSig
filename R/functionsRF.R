listToMatrix <- function(df){
        message("Converting Network to Adjacency Matrix...")
        na_impute=0
        mat <- matrix(0, length(unique(unlist(df[,2]))), length(unique(unlist(df[,1]))))
        #mat <- Matrix(0, nrow = length(unique(unlist(df[,2]))), ncol = length(unique(unlist(df[,1]))),sparse=TRUE)
        rownames(mat)=sort(unique(unlist(df[,2])))
        colnames(mat)=sort(unique(unlist(df[,1])))
        #mat[]=na_impute
        #z=nrow(df)
        #lapply(1:nrow(df),function (x){mat[df[x,2],df[x,1]]=df[x,3];y=x/z;message (y)})
        #for (x in  1:nrow(df)){mat[df[x,2],df[x,1]]=df[x,3];y=x/z;message (y)}
        mat[df[,2:1]] <- as.numeric(df[,3])
        message("Done.")
        return(mat)
        #for (x in 1:nrow(df)){ mat[as.character(df[x,2]),as.character(df[x,1])]=df[x,3];message (x)}
}
