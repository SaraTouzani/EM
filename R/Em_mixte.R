Em_mixte <-
function(Xquanti, Xquali, meani=NULL, sigmai=NULL, alpha_ki=NULL,Xijh_i=NULL, pki=NULL,Y=NULL, methode="init",e=0.1,nbr_iteration=NULL,k) {
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
