fk_quali <-
function(Y,Xijh,alpha_ki,k=length(alpha_ki),n=nrow(Xquali), p=length(colnames(Xquali)) ){
  if(p==0){
    p<-1
    n<-length(Xquali)
  }
  fk_quali_ik<-matrix(1,nrow = n, ncol = length(unique(Y)) )
  K=length(unique(Y))
  fkx_q<-as.list(seq_len(length(Y)))
  for (i in 1:n) {
    for (k in 1:K) {
      for (j in 1:p) {
        fk_quali_ik[i,k]<-fk_quali_ik[i,k] * (as.numeric(alpha_ki[[k]][[j]][["proba"]]%*%as.numeric (Xijh[[j]][i,-length(Xijh[[j]][i,])]))) 
      }
    }
    fkx_q[[i]]<-fk_quali_ik[i,]
  }

   return(fkx_q)
 }
