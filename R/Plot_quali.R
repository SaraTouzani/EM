Plot_quali <-
function(Dataquali,Em_results,all=TRUE){
  tkx <-Em_results[[length(Em_results)]][["part_MAP"]]
  Y<-apply(tkx, MARGIN = 1 ,which.max)
  if(all){
    layout(matrix(1:ncol(Dataquali),ncol = floor(sqrt(ncol(Dataquali)))))
    names_vars<-colnames(Dataquali)
    for (i in 1:ncol(Dataquali)) {
      my_variable<-table(Dataquali[,i],Y)
      my_variable2<-prop.table(my_variable,1)
      barplot(my_variable2,beside=TRUE,legend=paste("cluster",rownames(my_variable2)), ylab="proportion du cluster ",xlab=names_vars[i],col=1:length(unique(Y)))
    }

  }
  else{
    names_vars<-colnames(Dataquali)
    for (i in 1:ncol(Dataquali)) {
      my_variable<-table(Y,Dataquali[,i])
      my_variable2<-prop.table(my_variable,1)
      barplot(my_variable2,beside=TRUE,legend=paste("cluster",rownames(my_variable2)), ylab="proportion du cluster ",xlab=names_vars[i],col=1:length(unique(Y)))
    }
  }
  
}
