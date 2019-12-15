#install.packages("qlcMatrix")  
library(psych)
library(qlcMatrix)
library(FactoMineR)
library(MASS)
library(klaR)
library(ggplot2)
##############################

##########
data("iris")

iris$Sepal.Length_cat[iris$Sepal.Length<=5]<-" tres petit"
iris$Sepal.Length_cat[iris$Sepal.Length<=6 & iris$Sepal.Length>5]<-"petit"
iris$Sepal.Length_cat[iris$Sepal.Length<=7 & iris$Sepal.Length>6]<-"moyen"
iris$Sepal.Length_cat[iris$Sepal.Length>7]<-"grand"


Xquali<-iris[,5:6]
Xquanti<-iris[,1:4]
#######################


Em_mixte <- function(Xquanti, Xquali, meani=NULL, sigmai=NULL, alpha_ki=NULL,Xijh_i=NULL, pki=NULL,Y=NULL, methode="init",e=0.1,nbr_iteration=NULL,k) {
  #  Yobs<-Y[!is.na(Y)]
  #  Ymis<-Y[is.na(Y)]
  #  n<-length(c(Yobs,Ymis))
  #  o<-length(Yobs)
  if(methode=="kmeans"){
    Y<-kmeans(Xquanti,k)$cluster
  }
  else{
    if(methode=="kmodes"){
      Y<-kmodes(Xquali,k)$cluster
    }
    else{
      Y<-init_class(iris,k)
    } 
  }
    if(is.null(meani)){
    meani<-aggregate(Xquanti, by=list(Y), FUN=mean)
  }
  if(is.null(alpha_ki)){
    alpha_ki<-alphai(Y,Xquali)
  }
  if(is.null(Xijh_i)){
    Xijh_i<-xijh(Xquali)
  }
  if(is.null(pki)){
    pki<-as.numeric(table(Y)/length(Y))
  }
  
  if(is.null(sigmai)){
    sigmai<-sigma_init(Xquanti,Y)
    Test<-TRUE
    for(i in 1:length(sigmai)){
      if((det(as.matrix(sigmai[[i]]))**(0.5)) < 1e-3 || is.na(det(as.matrix(sigmai[[2]]))**0.5) ) Test<-FALSE
    }
    if(!Test){
      res<-PCA(Xquanti, scale.unit = TRUE, ncp = Inf, graph = FALSE)
      X<-(res$ind[["coord"]])
      sigmai<-sigma_init(as.data.frame(X),Y) 
      o<-0
      vars<-colnames(Xquanti)
      VAR<-NULL
      XQQ<-NULL
      for(j in 1:ncol(Xquanti)){
        t<-TRUE
        for (i in 1:k) {
          if(sum(abs(sigmai[[i]][[j]]))<0.1){
            # Xquanti[,j]<-NULL
            # j<-j+1
            o<-o+1
            t<-FALSE
          }
        }
        if(t){
          XQQ<-cbind(XQQ,X[,j])
          #VAR<-c(VAR,vars[j])
        }
      }
      for (j in 2:ncol(XQQ)) {
        sigmaii<-sigma_init(as.data.frame(XQQ[,1:j]),Y)
        t<-TRUE
        for (i in 1:k) {
          if(!is.na(det(as.matrix(sigmaii[[k]]))**(0.5))){
            if((det(as.matrix(sigmaii[[k]]))**(0.5))<0.1 ){
              
              t<-FALSE
            }       
          }
          else{
            t<-FALSE
          }
          
          if(t){
            XX<-XQQ[,1:j]
          }
          
        }
        if(!t){
          break
        }
      }
      sigmaiii<-sigma_init(as.data.frame(XX),Y)
      sigmai<-sigma_init(as.data.frame(XQQ),Y)
      sigmai<-sigmaiii
      meani<-aggregate(XX, by=list(Y), FUN=mean)
      #(1/ ( ((2*pi)**(p/2))* (det(as.matrix(sigma[[k]]))**(0.5))  ) ) *exp(t(as.matrix ((Xquanti[i,]- as.numeric( mean[k,-1])))) %*% ginv(as.matrix(sigma[[k]])) %*% (as.matrix ((Xquanti[i,]- as.numeric( mean[k,-1])))))
      if(is.null(pki)){
        pki<-as.numeric(table(Y)/length(Y))
      }
      Meanki <- meani
      Sigki <- sigmai
      Pki <- pki
      Alphaki <- alpha_ki
      
      
      fkx<-fk2(XX, meani, sigmai, Y, Xijh_i , alpha_ki, p=length(Xquali) )
        #fk_quanti2(XX,meani,sigmai,Y)
      f<-fx(Pki,fkx,Y)
      tkx<-tk(fkx,Pki)  
      li <- l(f)
      Qi<-Qt(tkx,fkx,Pki)
      BIC<-Bic(f,Y,Xquanti = XX)
      ICL<-BIC-sum((matrix(unlist(fkx),ncol=k, byrow = TRUE))*log(matrix(unlist(fkx),ncol=k, byrow = TRUE)))
      i <- 1
      #  teta_init<-list(iteration=i,
      #                   Mean=Meanki,
      #                   Sigma =Sigki,
      #                   Proba=Pki,
      #                   Alphak = Alphaki,
      #                   log_like=li,
      #                   Q=Qi,
      #                   Bic=BIC,
      #                   ICL=ICL,
      #                   proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
      #                   part_MAP=tkx  
      # )
      #plot(Xquanti$Sepal.Length,Xquanti$Sepal.Width, col = Y, pch = 19)
      #Y[1]<-which.max(tkx[1,])
      teta<-list()
      # teta[[1]]<-list(iteration=1,
      #                 Mean=meani,
      #                 Sigma =Sigki,
      #                 Proba=pki,
      #                 log_like=li,
      #                 Q=Qi,
      #                 Bic=BIC,
      #                 Icl=ICL,
      #                 proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
      #                 part_MAP=tkx          
      # )
      x11()
      #  teta[[1]]<-teta_init
      repeat {
        pairs.panels(XX,gap=0, bg=Y,pch = 21)
        Y<-apply(tkx, MARGIN = 1 ,which.max)
        Nk<-nk(tkx)
        Meank <- meank(Nk,tkx,XX,Y)
        Sigk <- sigmak(Nk,tkx,XX,Meanki)
        Alphak <- alphak(tkx,Xijh_i,Y,Xquali)
        Pk<-pk(Nk,XX)
        #   tk <- tk(pk, alphak)
        fkx<-fk2(XX, Meank, Sigk, Y, Xijh_i , Alphak, p=ncol(Xquali) )
          #fk_quanti2(XX,Meank,Sigk,Y)
        f<-fx(Pk,fkx,Y)
        lf<-l(f)
        tkx<-tk(fkx,Pk)  
        Qf<-Qt(tkx,fkx,Pk)
        BIC<-Bic(f,Y,Xquanti = XX)
        ICL<-BIC-sum(tkx*log1p(tkx))
        teta[[i]]<-list(iteration=i,
                          Mean=Meank,
                          Sigma =Sigk,
                          Proba=Pk,
                          log_like=lf,
                          Q=Qf,
                          Bic=BIC,
                          Icl=ICL,
                          proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
                          part_MAP=tkx          
        )
        if(!is.null(nbr_iteration)){
          if(i>=nbr_iteration) break
        }
        if (abs(li - lf) < e)  {
          break
        }
        else{
          Meanki <- as.data.frame(Meank)
          Sigki <- Sigk
          Pki <- Pk
          Alphaki <- Alphak
          
          li <- lf
          i<-i+1
        }
        
        
      }
      
    }
    
    else{
      Meanki <- meani
      Sigki <- sigmai
      Pki <- pki
      Alphaki <- alpha_ki
      
      fkx<-fk(Xquanti, meani, sigmai, Y, Xijh_i , alpha_ki, p=length(Xquali) )
      f<-fx(Pki,fkx,Y)
      tkx<-tk(fkx,Pki)  
      li <- l(f)
      Qi<-Qt(tkx,fkx,Pki)
      BIC<-Bic(f,Y,Xquali = Xquali,Xquanti = Xquanti)
      ICL<-BIC-sum((matrix(unlist(fkx),ncol=k, byrow = TRUE))*log(matrix(unlist(fkx),ncol=k, byrow = TRUE)))
      i <- 1
      #  teta_init<-list(iteration=i,
      #                   Mean=Meanki,
      #                   Sigma =Sigki,
      #                   Proba=Pki,
      #                   Alphak = Alphaki,
      #                   log_like=li,
      #                   Q=Qi,
      #                   Bic=BIC,
      #                   ICL=ICL,
      #                   proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
      #                   part_MAP=tkx  
      # )
      #plot(Xquanti$Sepal.Length,Xquanti$Sepal.Width, col = Y, pch = 19)
      #Y[1]<-which.max(tkx[1,])
      teta<-list()
      #  teta[[1]]<-teta_init
      repeat {
        pairs.panels(Xquanti,gap=0, bg=Y,pch = 21)
        Y<-apply(tkx, MARGIN = 1 ,which.max)
        Nk<-nk(tkx)
        Meank <- meank(Nk,tkx,Xquanti,Y)
        Sigk <- sigmak(Nk,tkx,Xquanti,Meanki)
        Pk<-pk(Nk,Xquanti)
        Alphak <- alphak(tkx,Xijh_i,Y,Xquali)
        #   tk <- tk(pk, alphak)
        fkx<-fk(Xquanti, Meank, Sigk, Y, Xijh_i , Alphak, p=ncol(Xquali) )
        f<-fx(Pk,fkx,Y)
        lf<-l(f)
        tkx<-tk(fkx,Pk)  
        Qf<-Qt(tkx,fkx,Pk)
        BIC<-Bic(f,Y,Xquali = Xquali,Xquanti = Xquanti)
        ICL<-BIC-sum(tkx*log1p(tkx))
        teta[[i]]<-list(iteration=i,
                        Mean=Meank,
                        Sigma =Sigk,
                        Proba=Pk,
                        Alphak = Alphak,
                        log_like=lf,
                        Q=Qf,
                        Bic=BIC,
                        Icl=ICL,
                        proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
                        part_MAP=tkx          
        )
        if(!is.null(nbr_iteration)){
          if(i>=nbr_iteration) break
        }
        if (abs(li - lf) < e)  {
          break
        }
        else{
          Meanki <- Meank
          Sigki <- Sigk
          Pki <- Pk
          Alphaki <- Alphak
          
          li <- lf
          i<-i+1
        }
        
        
      }
    }
    
    
  }
  else{
    Test<-TRUE
    for(i in 1:length(sigmai)){
      if((det(as.matrix(sigmai[[i]]))**(0.5)) < 1e-3 || is.na(det(as.matrix(sigmai[[2]]))**0.5) ) Test<-FALSE
    }
    if(!Test){
      res<-PCA(Xquanti, scale.unit = TRUE, ncp = Inf, graph = FALSE)
      X<-(res$ind[["coord"]])
      sigmai<-sigma_init(as.data.frame(X),Y) 
      o<-0
      vars<-colnames(Xquanti)
      VAR<-NULL
      XQQ<-NULL
      for(j in 1:ncol(Xquanti)){
        t<-TRUE
        for (i in 1:k) {
          if(sum(abs(sigmai[[i]][[j]]))<0.1){
            # Xquanti[,j]<-NULL
            # j<-j+1
            o<-o+1
            t<-FALSE
          }
        }
        if(t){
          XQQ<-cbind(XQQ,X[,j])
          #VAR<-c(VAR,vars[j])
        }
      }
      for (j in 2:ncol(XQQ)) {
        sigmaii<-sigma_init(as.data.frame(XQQ[,1:j]),Y)
        t<-TRUE
        for (i in 1:k) {
          if(!is.na(det(as.matrix(sigmaii[[k]]))**(0.5))){
            if((det(as.matrix(sigmaii[[k]]))**(0.5))<0.1 ){
              
              t<-FALSE
            }       
          }
          else{
            t<-FALSE
          }
          
          if(t){
            XX<-XQQ[,1:j]
          }
          
        }
        if(!t){
          break
        }
      }
      sigmaiii<-sigma_init(as.data.frame(XX),Y)
      sigmai<-sigma_init(as.data.frame(XQQ),Y)
      sigmai<-sigmaiii
      meani<-aggregate(XX, by=list(Y), FUN=mean)
      #(1/ ( ((2*pi)**(p/2))* (det(as.matrix(sigma[[k]]))**(0.5))  ) ) *exp(t(as.matrix ((Xquanti[i,]- as.numeric( mean[k,-1])))) %*% ginv(as.matrix(sigma[[k]])) %*% (as.matrix ((Xquanti[i,]- as.numeric( mean[k,-1])))))
      if(is.null(pki)){
        pki<-as.numeric(table(Y)/length(Y))
      }
      Meanki <- meani
      Sigki <- sigmai
      Pki <- pki
      Alphaki <- alpha_ki
      
      
      fkx<-fk2(XX, meani, sigmai, Y, Xijh_i , alpha_ki, p=length(Xquali) )
      #fk_quanti2(XX,meani,sigmai,Y)
      f<-fx(Pki,fkx,Y)
      tkx<-tk(fkx,Pki)  
      li <- l(f)
      Qi<-Qt(tkx,fkx,Pki)
      BIC<-Bic(f,Y,Xquanti = XX)
      ICL<-BIC-sum((matrix(unlist(fkx),ncol=k, byrow = TRUE))*log(matrix(unlist(fkx),ncol=k, byrow = TRUE)))
      i <- 1
      #  teta_init<-list(iteration=i,
      #                   Mean=Meanki,
      #                   Sigma =Sigki,
      #                   Proba=Pki,
      #                   Alphak = Alphaki,
      #                   log_like=li,
      #                   Q=Qi,
      #                   Bic=BIC,
      #                   ICL=ICL,
      #                   proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
      #                   part_MAP=tkx  
      # )
      #plot(Xquanti$Sepal.Length,Xquanti$Sepal.Width, col = Y, pch = 19)
      #Y[1]<-which.max(tkx[1,])
      teta<-list()
      # teta[[1]]<-list(iteration=1,
      #                 Mean=meani,
      #                 Sigma =Sigki,
      #                 Proba=pki,
      #                 log_like=li,
      #                 Q=Qi,
      #                 Bic=BIC,
      #                 Icl=ICL,
      #                 proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
      #                 part_MAP=tkx          
      # )
      x11()
      #  teta[[1]]<-teta_init
      repeat {
        pairs.panels(XX,gap=0, bg=Y,pch = 21)
        Y<-apply(tkx, MARGIN = 1 ,which.max)
        Nk<-nk(tkx)
        Meank <- meank(Nk,tkx,XX,Y)
        Sigk <- sigmak(Nk,tkx,XX,Meanki)
        Alphak <- alphak(tkx,Xijh_i,Y,Xquali)
        Pk<-pk(Nk,XX)
        #   tk <- tk(pk, alphak)
        fkx<-fk2(XX, Meank, Sigk, Y, Xijh_i , Alphak, p=ncol(Xquali) )
        #fk_quanti2(XX,Meank,Sigk,Y)
        f<-fx(Pk,fkx,Y)
        lf<-l(f)
        tkx<-tk(fkx,Pk)  
        Qf<-Qt(tkx,fkx,Pk)
        BIC<-Bic(f,Y,Xquanti = XX)
        ICL<-BIC-sum(tkx*log1p(tkx))
        teta[[i]]<-list(iteration=i,
                        Mean=Meank,
                        Sigma =Sigk,
                        Proba=Pk,
                        log_like=lf,
                        Q=Qf,
                        Bic=BIC,
                        Icl=ICL,
                        proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
                        part_MAP=tkx          
        )
        if(!is.null(nbr_iteration)){
          if(i>=nbr_iteration) break
        }
        if (abs(li - lf) < e)  {
          break
        }
        else{
          Meanki <- as.data.frame(Meank)
          Sigki <- Sigk
          Pki <- Pk
          Alphaki <- Alphak
          
          li <- lf
          i<-i+1
        }
        
        
      }
      
    }
    
    else{
      Meanki <- meani
      Sigki <- sigmai
      Pki <- pki
      Alphaki <- alpha_ki
      
      fkx<-fk(Xquanti, meani, sigmai, Y, Xijh_i , alpha_ki, p=length(Xquali) )
      f<-fx(Pki,fkx,Y)
      tkx<-tk(fkx,Pki)  
      li <- l(f)
      Qi<-Qt(tkx,fkx,Pki)
      BIC<-Bic(f,Y,Xquali = Xquali,Xquanti = Xquanti)
      ICL<-BIC-sum((matrix(unlist(fkx),ncol=k, byrow = TRUE))*log(matrix(unlist(fkx),ncol=k, byrow = TRUE)))
      i <- 1
      #  teta_init<-list(iteration=i,
      #                   Mean=Meanki,
      #                   Sigma =Sigki,
      #                   Proba=Pki,
      #                   Alphak = Alphaki,
      #                   log_like=li,
      #                   Q=Qi,
      #                   Bic=BIC,
      #                   ICL=ICL,
      #                   proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
      #                   part_MAP=tkx  
      # )
      #plot(Xquanti$Sepal.Length,Xquanti$Sepal.Width, col = Y, pch = 19)
      #Y[1]<-which.max(tkx[1,])
      teta<-list()
      #  teta[[1]]<-teta_init
      repeat {
        pairs.panels(Xquanti,gap=0, bg=Y,pch = 21)
        Y<-apply(tkx, MARGIN = 1 ,which.max)
        Nk<-nk(tkx)
        Meank <- meank(Nk,tkx,Xquanti,Y)
        Sigk <- sigmak(Nk,tkx,Xquanti,Meanki)
        Pk<-pk(Nk,Xquanti)
        Alphak <- alphak(tkx,Xijh_i,Y,Xquali)
        #   tk <- tk(pk, alphak)
        fkx<-fk(Xquanti, Meank, Sigk, Y, Xijh_i , Alphak, p=ncol(Xquali) )
        f<-fx(Pk,fkx,Y)
        lf<-l(f)
        tkx<-tk(fkx,Pk)  
        Qf<-Qt(tkx,fkx,Pk)
        BIC<-Bic(f,Y,Xquali = Xquali,Xquanti = Xquanti)
        ICL<-BIC-sum(tkx*log1p(tkx))
        teta[[i]]<-list(iteration=i,
                        Mean=Meank,
                        Sigma =Sigk,
                        Proba=Pk,
                        Alphak = Alphak,
                        log_like=lf,
                        Q=Qf,
                        Bic=BIC,
                        Icl=ICL,
                        proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
                        part_MAP=tkx          
        )
        if(!is.null(nbr_iteration)){
          if(i>=nbr_iteration) break
        }
        if (abs(li - lf) < e)  {
          break
        }
        else{
          Meanki <- Meank
          Sigki <- Sigk
          Pki <- Pk
          Alphaki <- Alphak
          
          li <- lf
          i<-i+1
        }
        
        
      }
    }
  }
  



  # Alphak = Alphak   ,               #   tk <- tk(pk, alphak)
  #fkx<-fk(Xquanti, Meank, Sigk, Y, Xijh_i , Alphak, p=ncol(Xquali) )
  #f<-fx(Pk,fkx,Y)
  return(teta)
}

Em_quali <- function( Xquali,  alpha_ki=NULL,Xijh_i=NULL, pki=NULL,Y=NULL, methode="init",e=0.1,nbr_iteration=NULL,k) {
  #  Yobs<-Y[!is.na(Y)]
  #  Ymis<-Y[is.na(Y)]
  #  n<-length(c(Yobs,Ymis))
  #  o<-length(Yobs)
  if(methode=="kmodes"){
    Y<-kmodes(Xquali,k)$cluster
  }
  else{
    Y<-init_class(iris,k)
  }
  if(is.null(alpha_ki)){
    alpha_ki<-alphai(Y,Xquali)
  }
  if(is.null(Xijh_i)){
    Xijh_i<-xijh(Xquali)
  }
  if(is.null(pki)){
    pki<-as.numeric(table(Y)/length(Y))
  }
  Pki <- pki
  Alphaki <- alpha_ki
  p=length(colnames(Xquali))
  n<-nrow(Xquali)
  fkx<-fk_quali(Y,Xijh_i,alpha_ki, p=p, n = n)  
  #fkx<-fk_quali(Y,Xijh_i,alpha_ki)
  f<-fx(Pki,fkx,Y)
  tkx<-tk(fkx,Pki)  
  li <- l(f)
  Qi<-Qt(tkx,fkx,Pki)
  BIC<-Bic(f,Y,Xquali = Xquali)
  ICL<-BIC-sum(tkx*log(tkx))
  i <- 1
#  teta_init<-list(iteration=i,
#                  Proba=Pki,
#                  Alphak = Alphaki,
#                  log_like=li,
#                  Q=Qi,
#                  Bic=BIC,
#                  ICL=ICL,
#                  proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
#                  part_MAP=tkx  
#  )
  #plot(Xquanti$Sepal.Length,Xquanti$Sepal.Width, col = Y, pch = 19)
  #Y[1]<-which.max(tkx[1,])
  teta<-list()
  K<-k
 # teta[[1]]<-teta_init
#  x11()
  repeat {
    #layout(matrix(1:6, ncol =K ,nrow=2))
 #   vars<-colnames(Xquali)
#    #par(mfrow=c(length(Xquali),1))
 #   for (j in 1:length(Xquali)) {
  #    XQ_K<-NULL
   #   for(k in unique(Y)){
    #      XQ_K<-rbind(XQ_K,table(Xquali[Y==k,j]))
          #pie(table(Xquali[Y==k,j]),main =paste(c("Variable:", vars[j]), collapse=" ") )
         # text(-1, 0, paste(c("  Freq:", table(Xquali[Y==k,j])), collapse=" "), col = "black")
     # }
      #barplot(XQ_K,beside = FALSE,legend.text = TRUE,xlab = paste(c("Variable:", vars[j]), collapse=" ") )
    #}

    Y<-apply(tkx, MARGIN = 1 ,which.max)
    Nk<-nk(tkx)
    Pk<-pk(Nk,Xquali)
    Alphak <- alphak(tkx,Xijh_i,Y,Xquali)
    #   tk <- tk(pk, alphak)
    fkx<-fk_quali(Y,Xijh_i,Alphak, p=p, n = n)
    f<-fx(Pk,fkx,Y)
    lf<-l(f)
    tkx<-tk(fkx,Pk)  
    Qf<-Qt(tkx,fkx,Pk)
    BIC<-Bic(f,Y,Xquali = Xquali)
    ICL<-BIC-sum(tkx*log(tkx))
    teta[[i]]<-list(iteration=i,
                      Proba=Pk,
                      Alphak = Alphak,
                      log_like=lf,
                      Q=Qf,
                      Bic=BIC,
                      ICL=ICL,
                      proba_individu=fkx,
                      part_MAP=tkx          
    )
    if(!is.null(nbr_iteration)){
      if(i>=nbr_iteration) break
    }
    if (abs(li - lf) < e)  {
      break
    }
    else{

      Pki <- Pk
      Alphaki <- Alphak
      
      li <- lf
      i<-i+1
    }
    
    
  }
  
  # Alphak = Alphak   ,               #   tk <- tk(pk, alphak)
  #fkx<-fk(Xquanti, Meank, Sigk, Y, Xijh_i , Alphak, p=ncol(Xquali) )
  #f<-fx(Pk,fkx,Y)
  return(teta)
}

Em_quanti <- function(Xquanti, meani=NULL, sigmai=NULL, pki=NULL,Y=NULL, methode="init",e=0.1,nbr_iteration=NULL,k) {
  #  Yobs<-Y[!is.na(Y)]
  #  Ymis<-Y[is.na(Y)]
  #  n<-length(c(Yobs,Ymis))
  #  o<-length(Yobs)
  if(methode=="kmeans"){
    Y<-kmeans(Xquanti,k)$cluster
  }
  else{
    Y<-init_class(iris,k)
  }
  if(is.null(pki)){
    pki<-as.numeric(table(Y)/length(Y))
  }
  if(is.null(meani)){
    meani<-aggregate(Xquanti, by=list(Y), FUN=mean)
  }
  if(is.null(sigmai)){
    sigmai<-sigma_init(Xquanti,Y)
    Test<-TRUE
    for(i in 1:length(sigmai)){
      if((det(as.matrix(sigmai[[i]]))**(0.5)) < 1e-3 || is.na(det(as.matrix(sigmai[[2]]))**0.5) ) Test<-FALSE
    }
    #### vars problems 
    if(!Test){
      res<-PCA(Xquanti, scale.unit = TRUE, ncp = Inf, graph = FALSE)
      X<-(res$ind[["coord"]])
      sigmai<-sigma_init(as.data.frame(X),Y) 
      o<-0
      vars<-colnames(Xquanti)
      VAR<-NULL
      XQQ<-NULL
      for(j in 1:ncol(Xquanti)){
        t<-TRUE
        for (i in 1:k) {
          if(sum(abs(sigmai[[i]][[j]]))<0.1){
            # Xquanti[,j]<-NULL
            # j<-j+1
            o<-o+1
            t<-FALSE
          }
        }
        if(t){
          XQQ<-cbind(XQQ,X[,j])
          #VAR<-c(VAR,vars[j])
        }
      }
      for (j in 2:ncol(XQQ)) {
        sigmaii<-sigma_init(as.data.frame(XQQ[,1:j]),Y)
        t<-TRUE
        for (i in 1:k) {
          if(!is.na(det(as.matrix(sigmaii[[k]]))**(0.5))){
            if((det(as.matrix(sigmaii[[k]]))**(0.5))<0.1 ){
              
              t<-FALSE
            }       
          }
          else{
            t<-FALSE
          }
          
          if(t){
            XX<-XQQ[,1:j]
          }
          
        }
        if(!t){
          break
        }
      }
      sigmaiii<-sigma_init(as.data.frame(XX),Y)
      sigmai<-sigma_init(as.data.frame(XQQ),Y)
      sigmai<-sigmaiii
      meani<-aggregate(XX, by=list(Y), FUN=mean)
      #(1/ ( ((2*pi)**(p/2))* (det(as.matrix(sigma[[k]]))**(0.5))  ) ) *exp(t(as.matrix ((Xquanti[i,]- as.numeric( mean[k,-1])))) %*% ginv(as.matrix(sigma[[k]])) %*% (as.matrix ((Xquanti[i,]- as.numeric( mean[k,-1])))))
      if(is.null(pki)){
        pki<-as.numeric(table(Y)/length(Y))
      }
      Meanki <- meani
      Sigki <- sigmai
      Pki <- pki
      
      fkx<-fk_quanti2(XX,meani,sigmai,Y)
      f<-fx(Pki,fkx,Y)
      tkx<-tk(fkx,Pki)  
      li <- l(f)
      Qi<-Qt(tkx,fkx,Pki)
      BIC<-Bic(f,Y,Xquanti = XX)
      ICL<-BIC-sum((matrix(unlist(fkx),ncol=k, byrow = TRUE))*log(matrix(unlist(fkx),ncol=k, byrow = TRUE)))
      i <- 1
      #  teta_init<-list(iteration=i,
      #                   Mean=Meanki,
      #                   Sigma =Sigki,
      #                   Proba=Pki,
      #                   Alphak = Alphaki,
      #                   log_like=li,
      #                   Q=Qi,
      #                   Bic=BIC,
      #                   ICL=ICL,
      #                   proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
      #                   part_MAP=tkx  
      # )
      #plot(Xquanti$Sepal.Length,Xquanti$Sepal.Width, col = Y, pch = 19)
      #Y[1]<-which.max(tkx[1,])
      teta<-list()
      teta[[1]]<-list(iteration=1,
                      Mean=meani,
                      Sigma =Sigki,
                      Proba=pki,
                      log_like=li,
                      Q=Qi,
                      Bic=BIC,
                      Icl=ICL,
                      proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
                      part_MAP=tkx          
      )
      x11()
      #  teta[[1]]<-teta_init
      repeat {
        pairs.panels(XX,gap=0, bg=Y,pch = 21)
        Y<-apply(tkx, MARGIN = 1 ,which.max)
        Nk<-nk(tkx)
        Meank <- meank(Nk,tkx,XX,Y)
        Sigk <- sigmak(Nk,tkx,XX,Meanki)
        Pk<-pk(Nk,XX)
        #   tk <- tk(pk, alphak)
        fkx<-fk_quanti2(XX,Meank,Sigk,Y)
        f<-fx(Pk,fkx,Y)
        lf<-l(f)
        tkx<-tk(fkx,Pk)  
        Qf<-Qt(tkx,fkx,Pk)
        BIC<-Bic(f,Y,Xquanti = XX)
        ICL<-BIC-sum(tkx*log1p(tkx))
        teta[[i+1]]<-list(iteration=i+1,
                        Mean=Meank,
                        Sigma =Sigk,
                        Proba=Pk,
                        log_like=lf,
                        Q=Qf,
                        Bic=BIC,
                        Icl=ICL,
                        proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
                        part_MAP=tkx          
        )
        if(!is.null(nbr_iteration)){
          if(i>=nbr_iteration) break
        }
        if (abs(li - lf) < e)  {
          break
        }
        else{
          Meanki <- as.data.frame(Meank)
          Sigki <- Sigk
          Pki <- Pk
          Alphaki <- Alphak
          
          li <- lf
          i<-i+1
        }
        
        
      }
      
    }
    ### 
    else{
      Meanki <- meani
      Sigki <- sigmai
      Pki <- pki
      
      fkx<-fk_quanti(Xquanti,meani,sigmai,Y)
      f<-fx(Pki,fkx,Y)
      tkx<-tk(fkx,Pki)  
      li <- l(f)
      Qi<-Qt(tkx,fkx,Pki)
      BIC<-Bic(f,Y,Xquanti = Xquanti)
      ICL<-BIC-sum((matrix(unlist(fkx),ncol=k, byrow = TRUE))*log(matrix(unlist(fkx),ncol=k, byrow = TRUE)))
      i <- 1
      #  teta_init<-list(iteration=i,
      #                   Mean=Meanki,
      #                   Sigma =Sigki,
      #                   Proba=Pki,
      #                   Alphak = Alphaki,
      #                   log_like=li,
      #                   Q=Qi,
      #                   Bic=BIC,
      #                   ICL=ICL,
      #                   proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
      #                   part_MAP=tkx  
      # )
      #plot(Xquanti$Sepal.Length,Xquanti$Sepal.Width, col = Y, pch = 19)
      #Y[1]<-which.max(tkx[1,])
      teta<-list()
      #  teta[[1]]<-teta_init
      repeat {
        pairs.panels(Xquanti,gap=0, bg=Y,pch = 21)
        Y<-apply(tkx, MARGIN = 1 ,which.max)
        Nk<-nk(tkx)
        Meank <- meank(Nk,tkx,Xquanti,Y)
        Sigk <- sigmak(Nk,tkx,Xquanti,Meanki)
        Pk<-pk(Nk,Xquanti)
       #   tk <- tk(pk, alphak)
        #fk_quanti(Xquanti,Meank,Sigk,Y)
        fkx<-fk_quanti(Xquanti,Meank,Sigk,Y)
        f<-fx(Pk,fkx,Y)
        lf<-l(f)
        tkx<-tk(fkx,Pk)  
        Qf<-Qt(tkx,fkx,Pk)
        BIC<-Bic(f,Y,Xquanti = Xquanti)
        ICL<-BIC-sum(tkx*log1p(tkx))
        teta[[i]]<-list(iteration=i,
                        Mean=Meank,
                        Sigma =Sigk,
                        Proba=Pk,
                        log_like=lf,
                        Q=Qf,
                        Bic=BIC,
                        Icl=ICL,
                        proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
                        part_MAP=tkx          
        )
        if(!is.null(nbr_iteration)){
          if(i>=nbr_iteration) break
        }
        if (abs(li - lf) < e)  {
          break
        }
        else{
          Meanki <- Meank
          Sigki <- Sigk
          Pki <- Pk

          li <- lf
          i<-i+1
        }
        
        
      }
    }

  }
  else{
    Test<-TRUE
    for(i in 1:length(sigmai)){
      if((det(as.matrix(sigmai[[i]]))**(0.5)) < 1e-3 || is.na(det(as.matrix(sigmai[[2]]))**0.5) ) Test<-FALSE
    }
    if(!Test){
      res<-PCA(Xquanti, scale.unit = TRUE, ncp = Inf, graph = FALSE)
      X<-(res$ind[["coord"]])
      sigmai<-sigma_init(as.data.frame(X),Y) 
      o<-0
      vars<-colnames(Xquanti)
      VAR<-NULL
      XQQ<-NULL
      for(j in 1:ncol(Xquanti)){
        t<-TRUE
        for (i in 1:k) {
          if(sum(abs(sigmai[[i]][[j]]))<0.1){
            # Xquanti[,j]<-NULL
            # j<-j+1
            o<-o+1
            t<-FALSE
          }
        }
        if(t){
          XQQ<-cbind(XQQ,X[,j])
          #VAR<-c(VAR,vars[j])
        }
      }
      for (j in 2:ncol(XQQ)) {
        sigmaii<-sigma_init(as.data.frame(XQQ[,1:j]),Y)
        t<-TRUE
        for (i in 1:k) {
          if(!is.na(det(as.matrix(sigmaii[[k]]))**(0.5))){
            if((det(as.matrix(sigmaii[[k]]))**(0.5))<0.1 ){
              
              t<-FALSE
            }       
          }
          else{
            t<-FALSE
          }
          
          if(t){
            XX<-XQQ[,1:j]
          }
          
        }
        if(!t){
          break
        }
      }
      sigmaiii<-sigma_init(as.data.frame(XX),Y)
      sigmai<-sigma_init(as.data.frame(XQQ),Y)
      sigmai<-sigmaiii
      meani<-aggregate(XX, by=list(Y), FUN=mean)
      #(1/ ( ((2*pi)**(p/2))* (det(as.matrix(sigma[[k]]))**(0.5))  ) ) *exp(t(as.matrix ((Xquanti[i,]- as.numeric( mean[k,-1])))) %*% ginv(as.matrix(sigma[[k]])) %*% (as.matrix ((Xquanti[i,]- as.numeric( mean[k,-1])))))
      if(is.null(pki)){
        pki<-as.numeric(table(Y)/length(Y))
      }
      Meanki <- meani
      Sigki <- sigmai
      Pki <- pki
      
      fkx<-fk_quanti2(XX,meani,sigmai,Y)
      f<-fx(Pki,fkx,Y)
      tkx<-tk(fkx,Pki)  
      li <- l(f)
      Qi<-Qt(tkx,fkx,Pki)
      BIC<-Bic(f,Y,Xquanti = XX)
      ICL<-BIC-sum((matrix(unlist(fkx),ncol=k, byrow = TRUE))*log(matrix(unlist(fkx),ncol=k, byrow = TRUE)))
      i <- 1
      #  teta_init<-list(iteration=i,
      #                   Mean=Meanki,
      #                   Sigma =Sigki,
      #                   Proba=Pki,
      #                   Alphak = Alphaki,
      #                   log_like=li,
      #                   Q=Qi,
      #                   Bic=BIC,
      #                   ICL=ICL,
      #                   proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
      #                   part_MAP=tkx  
      # )
      #plot(Xquanti$Sepal.Length,Xquanti$Sepal.Width, col = Y, pch = 19)
      #Y[1]<-which.max(tkx[1,])
      teta<-list()
      teta[[1]]<-list(iteration=1,
                      Mean=meani,
                      Sigma =Sigki,
                      Proba=pki,
                      log_like=li,
                      Q=Qi,
                      Bic=BIC,
                      Icl=ICL,
                      proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
                      part_MAP=tkx          
      )
      x11()
      #  teta[[1]]<-teta_init
      repeat {
        pairs.panels(XX,gap=0, bg=Y,pch = 21)
        Y<-apply(tkx, MARGIN = 1 ,which.max)
        Nk<-nk(tkx)
        Meank <- meank(Nk,tkx,XX,Y)
        Sigk <- sigmak(Nk,tkx,XX,Meanki)
        Pk<-pk(Nk,XX)
        #   tk <- tk(pk, alphak)
        fkx<-fk_quanti2(XX,Meank,Sigk,Y)
        f<-fx(Pk,fkx,Y)
        lf<-l(f)
        tkx<-tk(fkx,Pk)  
        Qf<-Qt(tkx,fkx)
        BIC<-Bic(f,Y,Xquanti = XX)
        ICL<-BIC-sum(tkx*log1p(tkx))
        teta[[i+1]]<-list(iteration=i+1,
                          Mean=Meank,
                          Sigma =Sigk,
                          Proba=Pk,
                          log_like=lf,
                          Q=Qf,
                          Bic=BIC,
                          Icl=ICL,
                          proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
                          part_MAP=tkx          
        )
        if(!is.null(nbr_iteration)){
          if(i>=nbr_iteration) break
        }
        if (abs(li - lf) < e)  {
          break
        }
        else{
          Meanki <- as.data.frame(Meank)
          Sigki <- Sigk
          Pki <- Pk
          Alphaki <- Alphak
          
          li <- lf
          i<-i+1
        }
        
        
      }
      
    }
    ### 
    else{
      Meanki <- meani
      Sigki <- sigmai
      Pki <- pki
      
      fkx<-fk_quanti(Xquanti,meani,sigmai,Y)
      f<-fx(Pki,fkx,Y)
      tkx<-tk(fkx,Pki)  
      li <- l(f)
      Qi<-Qt(tkx,fkx,Pki)
      BIC<-Bic(f,Y,Xquanti = Xquanti)
      ICL<-BIC-sum((matrix(unlist(fkx),ncol=k, byrow = TRUE))*log(matrix(unlist(fkx),ncol=k, byrow = TRUE)))
      i <- 1
      #  teta_init<-list(iteration=i,
      #                   Mean=Meanki,
      #                   Sigma =Sigki,
      #                   Proba=Pki,
      #                   Alphak = Alphaki,
      #                   log_like=li,
      #                   Q=Qi,
      #                   Bic=BIC,
      #                   ICL=ICL,
      #                   proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
      #                   part_MAP=tkx  
      # )
      #plot(Xquanti$Sepal.Length,Xquanti$Sepal.Width, col = Y, pch = 19)
      #Y[1]<-which.max(tkx[1,])
      teta<-list()
      #  teta[[1]]<-teta_init
      repeat {
        pairs.panels(Xquanti,gap=0, bg=Y,pch = 21)
        Y<-apply(tkx, MARGIN = 1 ,which.max)
        Nk<-nk(tkx)
        Meank <- meank(Nk,tkx,Xquanti,Y)
        Sigk <- sigmak(Nk,tkx,Xquanti,Meanki)
        Pk<-pk(Nk,Xquanti)
        #   tk <- tk(pk, alphak)
        #fk_quanti(Xquanti,Meank,Sigk,Y)
        fkx<-fk_quanti(Xquanti,Meank,Sigk,Y)
        f<-fx(Pk,fkx,Y)
        lf<-l(f)
        tkx<-tk(fkx,Pk)  
        Qf<-Qt(tkx,fkx,Pk)
        BIC<-Bic(f,Y,Xquanti = Xquanti)
        ICL<-BIC-sum(tkx*log1p(tkx))
        teta[[i]]<-list(iteration=i,
                        Mean=Meank,
                        Sigma =Sigk,
                        Proba=Pk,
                        log_like=lf,
                        Q=Qf,
                        Bic=BIC,
                        Icl=ICL,
                        proba_individu=matrix(unlist(fkx),ncol=k, byrow = TRUE),
                        part_MAP=tkx          
        )
        if(!is.null(nbr_iteration)){
          if(i>=nbr_iteration) break
        }
        if (abs(li - lf) < e)  {
          break
        }
        else{
          Meanki <- Meank
          Sigki <- Sigk
          Pki <- Pk
          
          li <- lf
          i<-i+1
        }
        
        
      }
    }
    }
  # Alphak = Alphak   ,               #   tk <- tk(pk, alphak)
  #fkx<-fk(Xquanti, Meank, Sigk, Y, Xijh_i , Alphak, p=ncol(Xquali) )
  #f<-fx(Pk,fkx,Y)
  return(teta)
}




teta<-Em_mixte(Xquanti,Xquali,mean,sigma,alpha_ki,Xijh, Pk,Y,k=3)
teta_sans_init<-Em_mixte(Xquanti,Xquali,e=0.1,methode = "kmeans",k=3)
teta_sans_init_visaP<-Em_mixte(data_continuous,data_categ,e=10,methode = "kmeans",k=2)



teta_sans_init_quanti<-Em_quanti(Xquanti,methode = "kmeans",e=0.1,k=3)
teta_sans_init_quanti_visaP<-Em_quanti(data_continuous,methode = "kmeans",e=1,k=2)

teta_sans_init_quali<-Em_quali(Xquali,e=0.1,methode = "kmodes",k=3)
teta_sans_init_quali_visaP<-Em_quali(data_categ,e=10,nbr_iteration = 3,methode = "kmodes",k=2)


## engageml,nbcbptar,nbcb,mteparte,nbeparte, mteparlt,nbeparlt ,engageml,mtrejet : Y==1
## mtrejet,engageml, nbeparit,mteparit,nbeparte,mteparte,nbcb,nbcbptar Y==1
## mtrejet,nbeparte,mteparte, Y==2

# Xquanti$mtrejet<-NULL
# Xquanti$nbeparte <-NULL
# Xquanti$mteparte <-NULL

Q<-NULL
for (i in 1:length(teta_sans_init)) {
  Q<-c(Q,teta_sans_init[[i]][["Q"]])
}
plot(1:length(teta_sans_init),Q,type = 'l', main = "Q(teta,teta(q)) pour k=3 Données mixe d'iris",xlab = "Itérations : q")

#################
par(mfrow=c(1,1))
Q_quali<-NULL
for (i in 1:length(teta_sans_init_quali)) {
  Q_quali<-c(Q_quali,teta_sans_init_quali[[i]][["Q"]])
}
plot(1:length(teta_sans_init_quali),Q_quali,type = 'l', main = "Q(teta,teta(q)) pour les Variables catégorielles et  k=3",xlab = "Itérations : q")
#################
Q_quanti<-NULL
for (i in 1:length(teta_sans_init_quanti)) {
  Q_quanti<-c(Q_quanti,teta_sans_init_quanti[[i]][["Q"]])
}
plot(1:length(teta_sans_init_quanti),Q_quanti,type = 'l', main = "Q(teta,teta(q)) pour les Variables quantitatives et  k=3",xlab = "Itérations : q")
#############
Q_quanti_visa<-NULL
for (i in 1:length(teta_sans_init_quanti_visaP)) {
  Q_quanti_visa<-c(Q_quanti_visa,teta_sans_init_quanti_visaP[[i]][["Q"]])
}
plot(1:length(teta_sans_init_quanti_visaP),Q_quanti_visa,type = 'l', main = "Q(teta,teta(q)) pour les Variables quantitatives de VisaPremier et  k=2",xlab = "Itérations : q")
############
Q_quali_visa<-NULL
for (i in 1:length(teta_sans_init_quali_visaP)) {
  Q_quali_visa<-c(Q_quali_visa,teta_sans_init_quali_visaP[[i]][["Q"]])
}
plot(1:length(teta_sans_init_quali_visaP),Q_quali_visa,type = 'l', main = "Q(teta,teta(q)) pour les Variables qualitative de VisaPremier et  k=2",xlab = "Itérations : q")

##########
Q_mixte_visa<-NULL
for (i in 1:length(teta_sans_init_visaP)) {
  Q_mixte_visa<-c(Q_mixte_visa,teta_sans_init_visaP[[i]][["Q"]])
}
plot(1:length(teta_sans_init_visaP),Q_mixte_visa,type = 'l', main = "Q(teta,teta(q)) pour les Variables mixte de VisaPremier et  k=2",xlab = "Itérations : q")

ICL<-NULL
BIC<-NULL
K<-5
for (k in 2:K) {
  teta_sans_init<-Em_mixte(Xquanti,Xquali,e=0.1,methode = "kmeans",k=k)
  BIC<-c(BIC,teta_sans_init[[length(teta_sans_init)]][["Bic"]])
  ICL<-c(ICL,teta_sans_init[[length(teta_sans_init)]][["Icl"]])
 # kmeans(Xquanti,k=)
}
  
#> BIC
#[1] 304.60314 172.22355  78.10263 -63.36916
#> ICL
#[1]  200.63106   68.25147  -24.90362 -167.27449
#[1] 304.60314 172.22355  78.10263 -63.36916
plot(2:K,BIC,type = 'l', main = "BIC en fonction du nombre de clusters pour Em",xlab = "nombre de clusters : k")
plot(2:K,ICL,type = 'l', main = "ICL en fonction du nombre de clusters pour Em",xlab = "nombre de clusters : k")


mod <- Mclust(iris[,1:4])
icl(mod)

# loop through the grouping species to create initiale classes
#determining each elm class randomly
init_class<-function(data,k){
  p<-vector("list",1)
  p[[1]][1:k]<- (1/k)
  return(sample(k,nrow(data), replace = TRUE, prob =unlist(p) ))
}
Y<-kmeans(Xquanti,3)$cluster
Y<-init_class(iris,k=3)
#j<-1
#for(i in unique(iris$Species)) {
 # iris$class[iris$Species==i]<-j
  #j<-j+1
    
#}



#Xquanti_k<-as.list(unique(iris$class))
#j<-1
#for(i in unique(iris$class)) {
#  Xquanti_k[[j]]<-data.frame(MAT = Xquanti[iris$class==i,])
#  j<-j+1
#}

#Xquali_k1<-as.list(unique(Y))
#j<-1
#for(i in unique(Y)) {
#  Xquali_k1[[j]]<-data.frame(MAT = Xquali[Y==i,])
#  j<-j+1
#}
#################

#for (j in 1:length(alpha)) {
#  X<-as.numeric(table(Xquali[,j])/nrow(iris))
#  alpha[[j]]<-data_frame(proba=X,Var=alpha[[j]])
#}
#Y<-iris[,7]
mean<-aggregate(Xquanti, by=list(Y), FUN=mean) # mean par group meaninitiale
##### cov  par group de class sigmak initiale
sigma_init<-function(data,Y){
  sigma<-as.list(unique(Y))
  j<-1
  data$class<-Y
  for(i in unique(Y)) {
    Mat_cov<-(cov(data[data$class==i, -match("class", names(data))][,]))
    sigma[[j]]<-data.frame(MAT = Mat_cov)
    j<-j+1
  }
  return(sigma)
}
sigma<-sigma_init(Xquanti,Y)

####### alpha est une data frame qui contient pour chaque variable les probas de chaque modalités de cette dernière 
alphai<-function(Y, Xquali){
  #repartition des observations de variables qualitatives suivant
  Xquali_k<-as.list(unique(Y))
  if(is.null(dim(Xquali))) { # une seule variable 
    j<-1
    for(i in unique(Y)) {
      Xquali_k[[j]]<-data.frame(MAT = Xquali[Y==i])
      j<-j+1
    }
    
    alpha_k<-as.list(unique(Y))
    for (k in 1:length(alpha_k)) {
      variables<-vector("list",1)
      X<-as.numeric(table( Xquali_k[[k]] )/ nrow(Xquali_k[[k]]))
      variables[[1]]<-data.frame(proba=X,
                                 modalite=unlist(attributes(table( Xquali_k[[k]] )/ nrow(Xquali_k[[k]]))$dimnames))
      for (h in levels(Xquali)) {
        if(!(h %in% variables[[1]]$modalite)){
          mod_miss<-data.frame(proba=0,
                               modalite=h)
          variables[[1]]<-rbind( variables[[1]],mod_miss)
        }
      }
      variables[[1]]$modalite<-as.character( variables[[1]]$modalite)
      variables[[1]]<- variables[[1]][order( variables[[1]]$modalite , decreasing = F),]
      
      alpha_k[[k]]<- variables
      
    } 
  }
  else{
    j<-1
    for(i in unique(Y)) {
      Xquali_k[[j]]<-data.frame(MAT = Xquali[Y==i,])
      j<-j+1
    }
    alpha_k<-as.list(unique(Y))
    for (k in 1:length(alpha_k)) {
      variables<-as.list(colnames(Xquali))
      for (j in 1:length(variables)) {
        X<-as.numeric(table( Xquali_k[[k]][,j] )/ nrow(Xquali_k[[k]]))
        variables[[j]]<-data.frame(proba=X,
                                   Var=variables[[j]],
                                   modalite=unlist(attributes(table( Xquali_k[[k]][,j] )/ nrow(Xquali_k[[k]]))$dimnames))
        for (h in levels(Xquali[,j])) {
          if(!(h %in% variables[[j]]$modalite)){
            mod_miss<-data.frame(proba=0,
                                 Var=unique(variables[[j]]$Var),
                                 modalite=h)
            variables[[j]]<-rbind(variables[[j]],mod_miss)
          }
          
        }
        variables[[j]]$modalite<-as.character(variables[[j]]$modalite)
        variables[[j]]<-variables[[j]][order(variables[[j]]$modalite , decreasing = F),]
        
        
      }
      alpha_k[[k]]<-variables
      
    }  
    
  }
  return(alpha_k)
}
alpha_ki<-alphai(Y,Xquali)

###### matrice Xijh qui verifie que Xij=h tel que j est une variable qualitative et i est un element 
# donc on marquera 1 sinn 0
#h sont les modalités allant de 1 a mj de la variable j
xijh<-function(Xquali){
  if(is.null(dim(Xquali))){
    p<-1
    Xijh<-vector("list",1)
    Xquali<-as.matrix(Xquali)
  }
  else{
    p<-length(colnames(Xquali))
    Xijh<-as.list(colnames(Xquali))
  }
  for (j in 1:p) {
    
    Xih<-matrix(0, ncol = length (levels(Xquali[,j])), nrow = nrow(Xquali))
    modalite_ij<-matrix("",  nrow = nrow(Xquali))
    h_index<-1
    for (h in sort(levels(Xquali[,j]))) {
      for (i in 1:nrow(Xquali)) {
        if((Xquali[i,j]==h)){
          Xih[i,h_index]<-1
          modalite_ij[i]<-h
        }
        else{
          Xih[i,h_index]<-0
        }
        
      }
      h_index<-h_index+1
      
    }
    Xijh[[j]]<-data.frame(proba=Xih,
                          modalite=modalite_ij)
    
  }
  
  return(Xijh)
}


Xijh<-xijh(Xquali)

########### pk proba par classes
Pk<-as.numeric(table(Y)/length(Y))

############### Zik appartenance de i au cluster k
#Zik<-function(){
#  zik<-matrix(0,nrow = nrow(iris),ncol = length(unique(iris$class)))
#  for (i in 1:nrow(iris)) {
    #K<-1
   # for (k in unique(Y)) {
  #    if(Y[i]==k){
     #   zik[i,K]<-1
    #  }
   #   K<-K+1
  #  }
 # }
#  return(zik)
#}

#Z<-Zik()


############## fk(x) densité par cluster
 fk<-function(Xquanti,mean,sigma,Y,Xijh,alpha_k,K=length(sigma),n=nrow(Xquanti), p=length(colnames(Xquali)) ){
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
   #alpha_k[[1]][[2]][["proba"]] les probas de chaque modalité de la variable 2 par rapport au cluster 1 
   # as.numeric (Xijh[[1]][1,-length(Xijh[[1]][1,])])) #Xijh[[1]] represente la variable 1 
   #  Xijh[[2]][1,] presente (les appartenance par rapport au modalités) 
   #de l'observation i=1 pour la deuxieme variable
   fk_quati_ik<-matrix(1,nrow = n, ncol = length(unique(Y)) )
   p=length(colnames(Xquanti))
   if(p==0){
     p<-1
   }
   for (i in 1:n) {
     for (k in 1:K) {
        
       fk_quati_ik[i,k]<-(1/ ( ((2*pi)**(p/2))* (det(as.matrix(sigma[[k]]))**(0.5))  ) ) *exp((-0.5)*as.matrix(Xquanti[i,]- as.numeric( mean[k,-1])) %*% ginv(as.matrix(sigma[[k]])) %*% t( as.matrix(Xquanti[i,]- as.numeric( mean[k,-1]))))
     }
   }
   
   fkxi<-as.list(seq_len(nrow(Xquanti)))
   for (i in 1:n) {
     fkxi[[i]]<-data.frame(densite=fk_quali_ik[i,]*fk_quati_ik[i,])
   } 
  return(fkxi)
  }

 fkx<-fk(Xquanti,mean,sigma,Y,Xijh,alpha_ki )

 fk2<-function(Xquanti,mean,sigma,Y,Xijh,alpha_k,K=length(sigma),n=nrow(Xquanti), p=length(colnames(Xquali)) ){
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
   #alpha_k[[1]][[2]][["proba"]] les probas de chaque modalité de la variable 2 par rapport au cluster 1 
   # as.numeric (Xijh[[1]][1,-length(Xijh[[1]][1,])])) #Xijh[[1]] represente la variable 1 
   #  Xijh[[2]][1,] presente (les appartenance par rapport au modalités) 
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
 
fk_quali<-function(Y,Xijh,alpha_ki,k=length(alpha_ki),n=nrow(Xquali), p=length(colnames(Xquali)) ){
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
 
 fkx_qual_ik<-fk_quali(Y,Xijh,alpha_ki)
 
 fk_quanti<-function(Xquanti,mean,sigma,Y,k=length(sigma),n=nrow(Xquanti), p=length(colnames(Xquanti)) ){
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
 
 fkx_quanti_ik<-fk_quanti(Xquanti,mean,sigma,Y)

 fk_quanti2<-function(Xquanti,mean,sigma,Y,k=length(sigma),n=nrow(Xquanti), p=length(colnames(Xquanti)) ){
   K=length(unique(Y))
   fk_quanti_ik<-matrix(1,nrow = n, ncol = length(unique(Y)) )
   fkx_q<-as.list(seq_len(length(Y)))
   for (i in 1:n) {
     for (k in 1:K) {
       
       #fk_quanti_ik[i,k]<-  (1/ ( ((2*pi)**(p/2))* (det(as.matrix(sigma[[k]]))**(0.5))  ) ) *exp(t(as.matrix ((Xquanti[i,]- as.numeric( mean[k,-1])))) %*% ginv(as.matrix(sigma[[k]])) %*% (as.matrix ((Xquanti[i,]- as.numeric( mean[k,-1])))))
       fk_quanti_ik[i,k]<- (1/ ( ((2*pi)**(p/2))* (det(as.matrix(sigma[[k]]))**(0.5))  ) )* exp((Xquanti[i,]- as.numeric( mean[k,-1]))%*%ginv(as.matrix(sigma[[k]])) %*%as.matrix(Xquanti[i,]- as.numeric( mean[k,-1]) ) *(-0.5))
       }
     fkx_q[[i]]<-fk_quanti_ik[i,]
   }
   return(fkx_q)
 }
 
  
 fx<-function(Pk,fkx,Y){
   f<-matrix(0, nrow = length(Y) )
   for (i in 1:length(Y)) {
     f[i]<- as.numeric(unlist(fkx[[i]]))%*%Pk
   }
   
   return(f)
 }

 f<-fx(Pk,fkx,Y)
 ########## tk
 tk<-function(fkx,pk){
   K<-length(pk)
   n<-length(fkx)
   tk<-matrix(0,ncol = K,nrow = n)
   for (i in 1:n) {
     tk[i,]<-(as.numeric(unlist(fkx[[i]]))*pk)/(as.numeric(unlist(fkx[[i]]))%*%pk)[1,1] 
     
   }
   return(tk)
 }

tkx<-tk(fkx,Pk)  

######### nk nbr delements par cluster
nk<-function(tkx){
  k<-ncol(tkx)
  nk<-list(0, k)
  for (i in 1:k) {
    nk[i]<-sum(tkx[,i])
  }
  return(unlist(nk))
}
Nk<-nk(tkx)

########### pk proabilité de chaque cluster
pk<-function(Nk,Xquanti){
  return(Nk/nrow(Xquanti))
}

Pk<-pk(Nk,Xquanti)
############# mean k 

meank<-function(Nk,tkx,Xquanti,Y){
  K<-length(Nk)
  p<-ncol(Xquanti)
  moyenne<-matrix(0, ncol = p, nrow =K )
  moyenne[,1]<-unique(Y)
  for (k in 1:K) {
    for (j in 1:p) {
      moyenne[k,j]<-sum((tkx[,k]*Xquanti)[,j])/Nk[k]
    }
  }
  moyenne<-cbind(sort(unique(Y)),moyenne)
  return(moyenne)
}
Meank<-meank(Nk,tkx,Xquanti)

#### sigmak
sigmak<-function(Nk,tkx,Xquanti,Meank){
  K<-length(Nk)
  n<-nrow(Xquanti)
  p<-ncol(Xquanti)
  sigma<-as.list(seq_len(K))
  for (k in 1:K) {
    sig<-matrix(0,ncol = p , nrow = p)
    for (i in 1:n) {
      sig<-sig + ((tkx[i,k]) * t(Xquanti[i,]- Meank[k,-1])%*%as.matrix(Xquanti[i,]-Meank[k,-1])) 
    }
    sigma[[k]]<-sig/Nk[k]
  }
  return(sigma)
  
}

Sigmak<-sigmak(Nk,tkx,Xquanti,mean) # on utilise l'ancienne moyenne

############
alphak<-function(tkx,Xijh,Y,Xquali){
  
  alpha_k<-as.list(unique(Y))
  
  for (k in 1:length(alpha_k)) {
    variables<-as.list(seq_len(length(Xijh)))
    for (j in 1:length(variables)) {
      
      X<-vector("list",length(levels(Xquali[,j])))
      h_index<-1
      for (h in sort(levels(Xquali[,j]))) {
        s<-0
        for (i in 1:nrow(Xquali)) {
          s<-s+ tkx[i,k]*Xijh[[j]][i  ,-ncol(Xijh[[j]])][h_index]
        }
        X[[h_index]]<-s
        h_index<-h_index+1
      }
      variables[[j]]<-data.frame(proba=unlist(X)/sum(unlist(X))
                                 #,Var=variables[[j]],
                                # modalite=unlist(attributes(table( Xquali_k[[k]][,j] )/ nrow(Xquali_k[[k]]))$dimnames)
                                )
      
    }
    alpha_k[[k]]<-variables
    
  }
  return(alpha_k)
  
  
}

Alphak<-alphak(tkx,Xijh,Y,Xquali)
############### log_lik de teta
l<-function(f){
  return(sum(log(f)))
}
####### Q log de vraisemblance completé


Qt<-function(tkx,fkx,Pk){
  n<-nrow(tkx)
  K<-ncol(tkx)
  somme<-0
    for (i in 1:n) {
#      somme<-somme + as.numeric(unlist(fkx[[i]])%*%tkx[i,]) 
      somme<-somme+((log1p(Pk)+log1p(unlist(fkx[[i]])))%*%tkx[i,])
    }
   return(somme)
}

Q<-Qt(tkx,fkx,Pk)



Bic<-function(f,Y=NULL,nbr_par=NULL,Xquali=NULL,Xquanti=NULL,K=length(unique(Y))) {
  nbr_par<-0
  if(!is.null(Xquali)){
    p<-ncol(Xquali)
    
    for (j in 1:p) {
      nbr_par<-nbr_par+length(levels(Xquali[,j])) -1
    }
    nbr_par<-nbr_par*K
  }
  if(!is.null(Xquanti)){
    p<-ncol(Xquanti)
    nbr_par<- K*(p*(p+1)/2+p+1)-1 + nbr_par    
  }

  return(-2*sum(log(f)) + nbr_par*log(length(Y)))
}

BIC<-Bic(f,Y,Xquali,Xquanti)
ICL<-BIC-sum(tkx*log(tkx))

Plot_quanti<-function(Dataquanti,Em_results,all=TRUE){
  
  tkx <-Em_results[[length(Em_results)]][["part_MAP"]]
  Y<-apply(tkx, MARGIN = 1 ,which.max)
  if(all){
                                                                                                                                                                                                                                                                                                                                            pairs.panels(Dataquanti,gap=0, bg=Y,pch = 21)
  }
  else{
    i<-1
    while(i< ((ncol(Dataquanti))-6) ){
      pairs.panels(Dataquanti[,i:(i+5)],gap=0, bg=Y,pch = 21)
      i<-i+6
    }
    pairs.panels(Dataquanti[,i:ncol(Dataquanti)],gap=0, bg=Y,pch = 21)
    
  }

  
  }
Plot_quanti(data_continuous,teta_sans_init_visaP,all=FALSE)

Plot_quali<-function(Dataquali,Em_results,all=TRUE){
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
Plot_quali(data_categ,teta_sans_init_visaP,all = FALSE)
