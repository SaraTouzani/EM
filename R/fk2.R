fk2 <-
function(Xquanti,mean,sigma,Y,Xijh,alpha_k,K=length(sigma),n=nrow(Xquanti), p=length(colnames(Xquali)) ){
   if(p==0){
     p<-1
   }
   fk_quali_ik<-matrix(1,nrow = n, ncol = length(unique(Y)) )
   K=length(unique(Y))
   for (i in 1:n) {
     for (k in 1:K) {
       for (j in 1:p) {
         fk_quali_ik[i,k]<-fk_quali_ik[i,k] * (as.numeric(alpha_k[[k]][[j]][["proba"]]%*%as.numeric (Xijh[[j]][i,-length(Xijh[[j]][i,])]))) 
       }
     }
   }
   # alpha[[1]][["proba"]] 
   # alpha pour cluster k, variable j :on calcule la proba de chaque modalite 
   #alpha_k[[1]][[2]][["proba"]] les probas de chaque modalit<U+00E9> de la variable 2 par rapport au cluster 1 
   # as.numeric (Xijh[[1]][1,-length(Xijh[[1]][1,])])) #Xijh[[1]] represente la variable 1 
   #  Xijh[[2]][1,] presente (les appartenance par rapport au modalit<U+00E9>s) 
   #de l'observation i=1 pour la deuxieme variable
   fk_quati_ik<-matrix(1,nrow = n, ncol = length(unique(Y)) )
   p=length(colnames(Xquanti))
   if(p==0){
     p<-1
   }
   for (i in 1:n) {
     for (k in 1:K) {
       
       fk_quati_ik[i,k]<- (1/ ( ((2*pi)**(p/2))* (det(as.matrix(sigma[[k]]))**(0.5))  ) )* exp((Xquanti[i,]- as.numeric( mean[k,-1]))%*%ginv(as.matrix(sigma[[k]])) %*%as.matrix(Xquanti[i,]- as.numeric( mean[k,-1]) ) *(-0.5)) 
     }
   }
   
   fkxi<-as.list(seq_len(nrow(Xquanti)))
   for (i in 1:n) {
     fkxi[[i]]<-data.frame(densite=fk_quali_ik[i,]*fk_quati_ik[i,])
   } 
   return(fkxi)
 }
