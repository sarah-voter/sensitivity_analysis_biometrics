rm(list=ls())
library(MASS)
library(tidyverse)
##### Functions ####
{
  invlogit <- function(x)exp(x)/(1+exp(x))
  logit <- function(x)log(x/(1-x))
  lossfxn <- function(y,yhat)(y-yhat)^2
  
  makedata <- function(npts, nparams,coefs){
    covar <- matrix(0, ncol = nparams, nrow = nparams)
    covar <- coefs$rho^abs(row(covar)-col(covar))
    covar <- covar*coefs$scale_factor
    X <- mvrnorm(n=npts, mu = coefs$xmeans, Sigma=covar)
    #X <- mvrnorm(n=npts, mu = coefs$xmeans, Sigma=diag(coefs$rho*coefs$scale_factor,ncol=nparams ,nrow=nparams )) #TODO TODOHP!! change back!!
    
    X <- cbind(1,X)
    S <- rbinom(npts,1,prob = invlogit(X %*% coefs$sx - sum(c(1,coefs$xmeans)*coefs$sx)))
    #Y <- X %*% coefs$yx_mu + rnorm(npts,sd=(X%*%coefs$yx_sigma)^2)
    xfxn_yprob <- function(x)x#^2
    Yprob <- invlogit(
      xfxn_yprob(X%*%coefs$yx_sigma 
                 - sum(c(1,coefs$xmeans)*coefs$yx_sigma) 
      ))
    tilted_Yprob <- rep(NA,nrow(X))
    tilted_Yprob[S==0]  <- Yprob[S==0]*exp(coefs$eta)/((Yprob[S==0] *exp(coefs$eta))+ 1 - Yprob[S==0] )
    tilted_Yprob[S==1]  <- Yprob[S==1]       
    
    
    Y <-rbinom(npts,1,
               prob = tilted_Yprob)
    table(Y,S)
    df <- data.frame(Y,Yprob,tilted_Yprob,X,S)
    hx <- glm(
      formula(paste0("Y ~ ",paste0("X",1:nparams,collapse="+")," - 1"))#no intercept because X1 is already the intercept
      , data=df[df$S==1,],family=binomial)
    
    startplots <- F
    if(startplots){
      
      #X distribution S=0 vs S=1
      df %>%pivot_longer(cols=c(X2,X3))%>%#,X4,X5,X6,X7))%>%
        ggplot()+
        geom_histogram(aes(value,fill=as.factor(S)),alpha=0.5)+
        facet_wrap(~name,scales="free")
      
      
      
      #boxplots etc
      df%>% pivot_longer(cols=c(Yprob,tilted_Yprob))%>%
        
        ggplot(aes(y=(value),x=as.factor(S)))+
        geom_jitter(alpha=0.2)+
        geom_violin(alpha=0.5)+
        facet_wrap(~name)+
        ggtitle("Yprobs")
      
      
      
      #X2 vs Yprob
      df %>% pivot_longer(cols=c(Yprob,tilted_Yprob))%>%
        ggplot(aes(y=(value),x=X2,color=as.factor(S)))+
        geom_point(alpha=0.2)+
        facet_wrap(~name)+
        
        ggtitle("Yprobs vs X")
      
      
      #S vs loss
      df %>%group_by(S)%>%
        summarize(proploss=mean(loss))%>%
        ungroup()%>%
        ggplot(aes(x=as.factor(S),y=proploss))+
        geom_bar(stat="identity")
    }
    
    
    return(list(df=df, hx=hx))
  }
  
  est_risk <- function(eta, df){
    
    b_x_eta_preds <- unlist(sapply(which(df$S==0), function(i) 
      (lossfxn(1,df$hxpred_raw[i])*exp(eta)*df$gxpred[i] + 
         lossfxn(0,df$hxpred_raw[i])*(1-df$gxpred[i]))/(1+df$gxpred[i]*(exp(eta)-1))
    ))
    head(b_x_eta_preds)
    return(mean(b_x_eta_preds))
  }
  bootstrap_CI_risk <- function(eta, df, nboot=500,CIbounds=c(0.025,0.975),plotting=F){
    n0 <- sum(df$S==0)
    sampleMatrix <- sapply(1:nboot, function(i) sample(1:n0,n0,replace = T))  
    bs_ests <- sapply(1:nboot, function(b) est_risk(eta,df[df$S==0,][sampleMatrix[,b],]))
    if(plotting)hist(bs_ests,main=paste("nboot=",nboot))
    quantile(bs_ests,probs=CIbounds)
  }
}
  

##### set params ####
#npts <- 1000#2000#500#1e3
nptsBig <- 1e5 #for big dataset to numerically estimate true risk (only performed once)
nparams <- 6

nsims <- 1000

coefs <- list(
  #X correlation matrix structure (see makedata function for data generating mechanism)
  rho = 0.25,
  scale_factor=1,
  
  #X means
  xmeans = rep(1*c(-1,1),length.out=nparams),
  # yx_mu = c(0,rep(2*c(1,-1),length.out=nparams)),
  
  #P(S==1| invlogit(Xj) )
  sx = c(0, rep(0.5*c(1,-1),length.out=nparams)) ,#TODO add correlation structure

  # P(Y==1| invlogit(Xj) ) (prior to tilting)
  yx_sigma =c(0, rep(0.25*c(0.5,1.5,-2),length.out=nparams)),
  eta= c(-1, -.5,0,.5,1)# transportability assumption satisfied when equal to zero
)




savelist <- list()
plotting <- F

qq=1
nptsvec <- c(500,1000,2000)#,3000)#c(250,500,750,1000,1500)
set.seed(1)
for(i_npts in 1:length(nptsvec)){
  if(length(coefs$eta)>0){
    for(e_eta in 1:length(coefs$eta)){
      
      #set.seed(2)
      tmpcoefs <- coefs
      tmpcoefs$eta <- coefs$eta[e_eta]
      for(i_sim in 1:nsims ){
        {
          
          #make data
          # if(i_sim==30){stop()}
          #for psi_true
          if(i_sim==1){ #only estimate psi_true using the big dataset one time
            
            tmp<- makedata(nptsBig, nparams, tmpcoefs)
            dfBig <- tmp$df
            #define hx only once here for the big dataset
            hx <- tmp$hx
            
            dfBig <- dfBig %>%
            
              mutate(hxpred_raw = predict(hx,newdata=.,type="response")
              )%>% 
              mutate(loss = lossfxn(hxpred_raw,Y))#TODO ask jon whether it's okay to use prob
            
            riskS0Big <- mean(dfBig$loss[dfBig$S==0])
            
          }
          
          tmp <- makedata(nptsvec[i_npts], nparams, tmpcoefs)
          df <- tmp$df
          #hx <- tmp$hx
          
          # estimate risk in target pop ####
          
          df <- df %>%
            mutate(hxpred_raw = predict(hx,newdata=.,type="response")
            )%>% 
            mutate(loss = lossfxn(hxpred_raw,Y))
          #table(df$S)
          
          head(df)
          startplots <- F
          if(startplots){
            
            #X distribution S=0 vs S=1
            df %>%pivot_longer(cols=c(X2,X3))%>%#,X4,X5,X6,X7))%>%
              ggplot()+
              geom_histogram(aes(value,fill=as.factor(S)),alpha=0.5)+
              facet_wrap(~name,scales="free")
            
            #boxplots etc
            df%>% pivot_longer(cols=c(Yprob,tilted_Yprob))%>%
              
              ggplot(aes(y=(value),x=as.factor(S)))+
              geom_jitter(alpha=0.2)+
              geom_violin(alpha=0.5)+
              facet_wrap(~name)+
              ggtitle("Yprobs")
            
            #X2 vs Yprob
            df %>% pivot_longer(cols=c(Yprob,tilted_Yprob))%>%
              ggplot(aes(y=(value),x=X2,color=as.factor(S)))+
              geom_point(alpha=0.2)+
              facet_wrap(~name)+
              
              ggtitle("Yprobs vs X")
            
            
            #S vs loss
            df %>%group_by(S)%>%
              summarize(proploss=mean(loss))%>%
              ungroup()%>%
              ggplot(aes(x=as.factor(S),y=proploss))+
              geom_bar(stat="identity")
          }
          
        }
        
        {#estimate model performance
          
          riskS0 <- mean(df$loss[df$S==0])
          
          #g(X)= Phat(Y=1|X,S=1)
          gx <- glm(formula(paste0("Y ~ ",paste0("X",1:nparams,collapse="+")," - 1")),
                    family="binomial",
                    data=filter(df,S==1))
          gxpredsS0 <- predict(gx, newdata=filter(df,S==0),type="response") 
          df$gxpred <- NA
          df$gxpred[df$S==0] <- gxpredsS0#
          
          gxpredsS1 <- predict(gx, newdata=filter(df,S==1),type="response")
          df$gxpred[df$S==1] <- gxpredsS1
          
          estimate_eta <- F
          if(estimate_eta){ #code below is to estimate eta. otherwise assume it's known
            
            muS0 <- mean(df$Y[df$S==0])
            muS1 <- mean(df$Y[df$S==1])
            #estimate eta
            est_mu <- function(eta,gxpredsS0,muS0){
              muhats <- exp(eta)*gxpredsS0/(exp(eta)*gxpredsS0 + 1 - gxpredsS0)
              sum(muhats-muS0)
            }
            eta_vec <- seq(-2,0,length.out=50)
            eta_ests <- data.frame(t(sapply(eta_vec, function(eta_i) c("eta"=eta_i,"stat"=est_mu(eta_i, gxpredsS0, muS0 )))))
            selected_eta <- eta_ests$eta[which(abs(eta_ests$stat)== min(abs(eta_ests$stat)))]
            
            if(plotting) plot(eta_ests$eta, abs(eta_ests$stat));abline(v=coefs$eta,col="blue"); abline(v=selected_eta,col="red")
            selected_eta
            bias_eta <- coefs$eta-selected_eta
          }
          
          sensitivity_analysis <- T
          if(sensitivity_analysis){
            risk_est_true_eta <- est_risk(eta=tmpcoefs$eta, df)
            risk_est_true_eta
            est_risk_augmented <- function(eta,df){#expression 8 on p 12 of sensitivity analysis paper
              truncate<- TRUE;tprob=0.1
              if(truncate){
                df$gxpred <- pmax(tprob,df$gxpred)
              }
              df$bhatx_preds <- unlist(sapply(1:nrow(df), function(i) 
                (lossfxn(1,df$hxpred_raw[i])*exp(eta)*df$gxpred[i] + #checked
                   lossfxn(0,df$hxpred_raw[i])*(1-df$gxpred[i]))/(1+df$gxpred[i]*(exp(eta)-1))#checked
              ))
              #phat(X) is estimator for P(S=1|X)
              phatx <- glm(formula(paste0("S ~ ",paste0("X",1:nparams,collapse="+")," - 1")), #checked
                           data=df)
              df$phatx_preds<- NA 
             
              df$phatx_preds[df$S==1] <- predict(phatx,newdata = filter(df,S==1),type="response" )#checked
              if(truncate){
                df$phatx_preds[df$S==1] <- pmax(tprob, df$phatx_preds[df$S==1] )
              }
              df$aug_est_preds <- NA
              df$aug_est_preds[df$S==0] <- df$bhatx_preds[df$S==0]#checked
              df$aug_est_preds[df$S==1] <-sapply(which(df$S==1), function(i)
                ( ( (1-df$phatx_preds[i])*exp(eta*df$Y[i]) ) / #checked
                    ( df$phatx_preds[i]*(exp(eta)*df$gxpred[i] + 1 - df$gxpred[i]))  ) *#checked
                  ( lossfxn(df$Y[i], df$hxpred_raw[i]) - df$bhatx_preds[i] )#checked
              )
              aug_est <- sum(df$aug_est_preds)/sum(df$S==0)#checked
              return(aug_est)
            }
            risk_est_AUG_true_eta <- est_risk_augmented(eta=tmpcoefs$eta, df)
            risk_est_AUG_true_eta
            #print(data.frame(big=riskS0Big, thisdata=riskS0, est_true_eta=risk_est_true_eta))
            
            #sensitivity analysis (iterate over range of etas)
            full_sensitivity_analysis <- F #perform bootstrapping to get CIs for a set of etas in a range (slow)
            
            if(full_sensitivity_analysis){
              etas <- c(tmpcoefs$eta,seq(-5,0, length.out=10))
              
              risk_est_eta_range <- sapply(etas, function(eta_i) 
                est_risk(eta=eta_i,df))
              risk_CIs_eta_range <- sapply(etas, function(eta_i) 
                bootstrap_CI_risk(eta=eta_i,df=df, nboot=50 ))
              #plot(etas, risk_est_eta_range);abline(h=riskS0Big,col="green");abline(h=riskS0,col="cyan");abline(h=risk_est_true_eta,col="red");abline(v=coefs$eta,col='black')
              ggeta <- data.frame(eta=etas, m = risk_est_eta_range, 
                                  l=risk_CIs_eta_range[1,],h=risk_CIs_eta_range[2,])
              ggeta %>% ggplot()+
                geom_point(aes(x=eta, y=m))+
                geom_errorbar(aes(x=eta, ymin=l,ymax=h))+
                geom_hline(yintercept = riskS0Big,color="forestgreen",alpha=0.5)+
                #geom_hline(yintercept = risk_est_true_eta,color="blue",alpha=0.5)+
                geom_text(aes(y=riskS0Big,x=0),vjust=0,label=expression(psi[true]),color="forestgreen")
              #geom_text(aes(y=risk_est_true_eta,x=0),vjust=0,label=expression(widehat(psi[eta])),color="blue")
            }
            
            
          }
          
          
          
          #savelist[[i_sim]] <- data.frame(selected_eta=selected_eta, bias_eta=bias_eta, true_eta=tmpcoefs$eta)
          savelist[[qq]] <- data.frame(est_risk=risk_est_true_eta,
                                       est_risk_aug = risk_est_AUG_true_eta,
                                       true_risk_big = riskS0Big,
                                       true_risk_i = riskS0,
                                       i=i_sim,
                                       eta=tmpcoefs$eta,
                                       npts=nptsvec[i_npts]
                                      
          )
          
          if((i_sim/round(nsims/20) ==i_sim%/%round(nsims/20) ))print(paste("sim", i_sim, "of",nsims))
          qq=qq+1
        }
      }
    }
  }
  
}

resultsdf <- as.data.frame(do.call(rbind,savelist))
#resultsdf <- read_rds("~/Desktop/counterfactual_prediction_models2023/resultsdfsave1.RDS")
#saveRDS(resultsdf,"~/Desktop/counterfactual_prediction_models2023/resultssave3.rds")
head(resultsdf);tail(resultsdf)
risk_true <- resultsdf$true_risk_big[1]
summtable<- resultsdf %>% group_by(eta,npts)%>%
  mutate(mean_est_risk = mean(est_risk),
         median_est_risk = median(est_risk),
         sd_est_risk = sd(est_risk),
         mean_est_aug_risk = mean(est_risk_aug),
         median_est_aug_risk = median(est_risk_aug),
         sd_est_aug_risk = sd(est_risk_aug),
         )%>%
  ungroup()%>%
  select(starts_with("mean") | starts_with("median") | starts_with("sd") | contains("risk_big")| contains("eta") | contains("npts"))%>%
  distinct()%>%
 
  mutate(bias_est_risk = true_risk_big-mean_est_risk,
         bias_est_aug_risk= true_risk_big-mean_est_aug_risk,
         bias_MED_est_risk = true_risk_big-median_est_risk,
         bias_MED_est_aug_risk= true_risk_big-median_est_aug_risk)

bff <- summtable#%>% select("eta","bias_est_aug_risk","sd_est_aug_risk")

# ff <- round(ff,4)
# #ff <- as.character(ff)
# ff
# ff$sample_size <- npts
# write.csv(ff,"~/Downloads/tmp3.csv")
# summtable$bias_est_risk/summtable$sd_est_risk
# summtable$bias_est_aug_risk/summtable$sd_est_aug_risk
# 
# ff1 <- read.csv("~/Downloads/tmp1.csv")
# ff2 <- read.csv("~/Downloads/tmp2.csv")
# ff3 <- read.csv("~/Downloads/tmp3.csv")
# bff0 <- rbind(ff1,ff2)
# bff <- rbind(bff0,ff3)

bff
bff %>% ggplot(aes(y=log(sd_est_aug_risk),x=npts,color=as.factor(eta)))+
  geom_point()+
  geom_line()+ggtitle("SD augmented estimator vs. npts (log scale)")
bff %>% ggplot(aes(y=log(sd_est_risk),x=npts,color=as.factor(eta)))+
  geom_point()+
  geom_line()+ggtitle("SD conditional loss estimator vs. npts (log scale)")
bff %>% ggplot(aes(y=bias_est_aug_risk,x=npts,color=as.factor(eta)))+
  geom_point()+
  geom_line()
bff %>% ggplot(aes(y=bias_MED_est_aug_risk,x=npts,color=as.factor(eta)))+
  geom_point()+
  geom_line()
bff %>% ggplot(aes(y=bias_MED_est_risk,x=npts,color=as.factor(eta)))+
  geom_point()+
  geom_line()
resultsdf %>% ggplot()+
  geom_histogram(aes(x=est_risk_aug))+
  facet_wrap(~ eta + npts,scales="free",ncol=5)

bff <- bff %>% arrange(eta,npts)
bff <- bff %>% select(eta, npts, bias_est_aug_risk, sd_est_aug_risk)
bff <- bff %>% filter(npts<3000)

df <- bff#summtable

colnames(df) <- c("eta", "Sample Size", "Bias","SD")
df <- as.data.frame(df)
#df$test <- "a"
makelatextable <- function(df,roundto=4){
  # df2 <- (sapply(1:ncol(df),function(j){
  #   jx=df[,j]
  #   ifelse(is.numeric(unlist(jx)) ,round(jx,roundto),jx)
  # } ))
  df2 <- df
  colnames(df2)<- colnames(df)
  for(j in 1:ncol(df2)){
    df2[,j] <- format(round(df2[,j],roundto),scientific=F)
  }
  df2 <- sapply(1:nrow(df2), function(i) paste0(paste0(df2[i,],collapse = " & ")," \\"))
  for(i in 1:length(df2)){
    if(i==1){
      cat("\\begin{tabular}{|")
      for(ii in 1:ncol(df)){
        cat("c|")
      }
      cat("}")
      cat("\n")
      for(ii in 1:ncol(df)){
        cat("\\textbf{",gsub("_"," ",colnames(df)[ii]),"}")
        cat(ifelse(ii==ncol(df), " \\\\"," & "))
      }
      cat(" \n \\hline")
      # cat("\n")
      # cat("\\newline")
    }
    cat(df2[i])
    cat("\\")
    cat("\n")
    ff="\\hline"
    cat(ff)
    cat("\n")
    if(i==length(df2)){
      cat("\\end{tabular}")
    }
  }
}
#makelatextable(select(df, contains("eta")|contains("bias")|contains("sd")))
makelatextable(df)
longdf <- resultsdf%>% pivot_longer(cols=contains("est") )
longdf %>% 
  #filter(name=="est_risk_aug")%>%
  ggplot()+
  geom_histogram(aes(x=value,fill=name))+
  geom_vline(aes(xintercept = true_risk_big))+
  facet_wrap(~eta,scales="free")+
  xlim(0,.2)


tmpresultsdf <- filter(resultsdf,eta==1)
hist(tmpresultsdf$est_risk); abline(v=risk_true,col="blue");abline(v=c(median(tmpresultsdf$est_risk),mean(tmpresultsdf$est_risk)),col="red");
abline(v=c(median(tmpresultsdf$est_risk_aug),mean(tmpresultsdf$est_risk_aug)),col="orange");
hist(tmpresultsdf$est_risk_aug,200); abline(v=risk_true,col="blue");abline(v=c(median(tmpresultsdf$est_risk),mean(tmpresultsdf$est_risk)),col="red");
which(tmpresultsdf$est_risk_aug<0)
mean(risk_true - tmpresultsdf$est_risk) /sd(tmpresultsdf$est_risk)

#do 1, .5,-0.5, -1
# head(resultsdf)
# resultsdf$selected_eta
# hist(resultsdf$selected_eta); abline(v=resultsdf$true_eta,col="blue"); abline(v=mean(resultsdf$selected_eta),col="red")










#misc plots
df%>% pivot_longer(cols=c(Yprob,tilted_Yprob))%>%
  
  ggplot(aes(y=(value),x=as.factor(S)))+
  geom_jitter(alpha=0.2)+
  geom_violin(alpha=0.5)+
  facet_wrap(~name)+
  ggtitle("Yprobs")



#X2 vs Yprob
df %>% pivot_longer(cols=c(Yprob,tilted_Yprob))%>%
  ggplot(aes(y=(value),x=X2,color=as.factor(S)))+
  geom_point(alpha=0.2)+
  facet_wrap(~name)+
  
  ggtitle("Yprobs vs X")


head(df)
table(df$hxpred)
par(mfrow=c(2,2))
summary(hx)





#old code for continuous Y ####
{
  # {
  #   npts <- 1e4
  #   nparams <- 2
  #   coefs <- list(
  #     sx = c(0, rep(1*c(-1,1),length.out=nparams)) ,
  #     yx_mu = c(0,rep(2*c(1,-1),length.out=nparams)),
  #     yx_sigma =c(0.25, rep(2*c(1,-1),length.out=nparams)) 
  #   )
  #   makedata <- function(npts, nparams,coefs){
  #     X <- mvrnorm(n=npts, mu = rep(0,nparams), Sigma=diag(1,nparams))
  #     X <- cbind(1,X)
  #     S <- rbinom(npts,1,prob = invlogit(X %*% coefs$sx))
  #     Y <- X %*% coefs$yx_mu + rnorm(npts,sd=(X%*%coefs$yx_sigma)^2)
  #     df <- data.frame(Y,X,S)
  #     hx <- lm(
  #       formula(paste0("Y ~ ",paste0("X",1:nparams,collapse="+")))
  #       , data=df[df$S==1,])
  #     return(list(df=df, hx=hx))
  #   }
  #   
  #   #make data
  #   #for psi_true
  #   nptsBig <- 1e5
  #   dfBig <- makedata(nptsBig, nparams, coefs)$df
  #   
  #   tmp <- makedata(npts, nparams, coefs)
  #   df <- tmp$df
  #   hx <- tmp$hx
  #   
  #   #df %>% group_by(S) %>% summarise(meanY=mean(Y))
  #  
  #   # df %>% ggplot()+
  #   #   geom_point(aes(x=X2,y=Y))+
  #   #   facet_wrap(~S)
  #   
  
  #   # estimate risk in target pop ####
  #   df <- df %>%
  #     mutate(hxpred = predict(hx,newdata=.))%>%
  #     mutate(loss = lossfxn(hxpred,Y))
  #   dfBig <- dfBig %>%
  #     mutate(hxpred = predict(hx,newdata=.))%>%
  #     mutate(loss = lossfxn(hxpred,Y))
  #   riskS0Big <- mean(dfBig$loss[df$S==0])
  #   riskS0 <- mean(df$loss[df$S==0])
  #   
  #   riskS0
  #   riskS0Big
  #   
  #   bx <- lm(formula(paste0("loss ~ ",paste0("X",1:nparams,collapse="+"))),
  #            data=df)
  #   df <- df %>%
  #     mutate(losspred = predict(bx,newdata=df))
  #   hist(df$losspred);abline(v=riskS0Big,col='red',lwd=2);abline(v=riskS0,col='orange',lwd=2);abline(v=median(df$losspred),col='brown',lty='dotted')
  #   
  #   }
  # par(mfrow=c(2,2))
  # summary(hx)
  
}

