fk_quanti <-
function(Xquanti,mean,sigma,Y,k=length(sigma),n=nrow(Xquanti), p=length(colnames(Xquanti)) ){
   K=length(unique(Y))
   fk_quanti_ik<-matrix(1,nrow = n, ncol = length(unique(Y)) )
   fkx_q<-as.list(seq_len(length(Y)))
   for (i in 1:n) {
     for (k in 1:K) {
       
       fk_quanti_ik[i,k]<-(1/ ( ((2*pi)**(p/2))* (det(as.matrix(sigma[[k]]))**(0.5))  ) ) *exp((-0.5)*as.matrix(Xquanti[i,]- as.numeric( mean[k,-1])) %*% ginv(as.matrix(sigma[[k]])) %*% t( as.matrix(Xquanti[i,]- as.numeric( mean[k,-1]))))
     }
     fkx_q[[i]]<-fk_quanti_ik[i,]
   }
   return(fkx_q)
 }
