graph.B <- function(donnee, param, objReport=NULL, msyVal=NULL, tacProj=NA, fProj=NA, valProj=NULL, ylimLog=NULL, pl=list(), plsd=list(), langue=c('fr','en','bil'), residus=TRUE,
                    ajouterProc=FALSE){
  switch(langue,
         'fr' = {labAn <- 'Ann\U{00E9}e'; labBiom <- "Biomasse ('000t)"; labRes <- 'R\U{00E9}sidus standardis\U{00E9}s'},
         'en' = {labAn <- 'Year'; labBiom <- "Biomass ('000t)"; labRes <- 'Standardized residuals'},
         'bil' = {labAn <- 'Ann\U{00E9}e/Year'; labBiom <- "Biomasse/Biomass ('000t)"; labRes <- 'R\U{00E9}sidus/Residuals'})
  if(!any(is.na(tacProj)) | !any(is.na(fProj))){
    nbAnProj <- max(length(tacProj[[1]]), length(fProj), na.rm=TRUE)
  }else{
    nbAnProj <- 0
  }
  ##
  plot(c(min(donnee$anneesFittees)-1,donnee$anneesFittees),
       objReport$Bpred/1000/1000,
       type='l', xlim=range(donnee$anneesFittees)+c(0,nbAnProj),
       ylim=c(0,max(objReport$Bpred,donnee$Bobs$valeur/objReport$qRel[donnee$Bobs$source]))/1000/1000,
       xlab=labAn, ylab=labBiom)
  if(ajouterProc){
    for(i in 1:length(objReport$Bpred.proc)){
      ## graphics::lines(c(min(donnee$anneesFittees)-1,donnee$anneesFittees)[i+c(0,1)], c(objReport$Bpred[i], objReport$Bpred.proc[i])/1000/1000)
      graphics::polygon(c(min(donnee$anneesFittees)-1,donnee$anneesFittees)[i+c(0,1,1)], c(objReport$Bpred[i], objReport$Bpred.proc[i], objReport$Bpred[i+1])/1000/1000, col=5, border=NA)
      graphics::lines(c(min(donnee$anneesFittees)-1,donnee$anneesFittees)[i+c(1,1)], c(objReport$Bpred[i+1], objReport$Bpred.proc[i])/1000/1000, col=1, lwd=1)
      ## graphics::points(donnee$anneesFittees[i], objReport$Bpred.proc[i]/1000/1000, pch=3)
    }
    graphics::lines(c(min(donnee$anneesFittees)-1,donnee$anneesFittees), objReport$Bpred/1000/1000)
  }
  graphics::abline(h=0, col='grey')
  if(!is.null(msyVal)){ # des points de r\U{00E9}f\U{00E9}rence sont ajout\U{00E9}s \U{00E0} 0.4 et 0.8 msy
    graphics::abline(h=msyVal$B/1000/1000, lty=2) #Bmsy
    graphics::axis(4, at=msyVal$B/1000/1000, label='RMS')
    graphics::axis(4, at=as.numeric(msyVal$B)/1000/1000*c(0.4,0.8), label=rep('',2))
    graphics::abline(h=as.numeric(msyVal$B)/1000/1000*c(0.4,0.8), lty=3, col='grey70') #Blim et Bsup
    graphics::axis(4, at=utils::tail(objReport$Bpred/1000/1000,1), label=round(utils::tail(objReport$Bpred,1)/msyVal$B,2), las=1)
  }
  for(i in unique(donnee$Bobs$source)){
    temp <- subset(donnee$Bobs, source==i)
    graphics::lines(donnee$anneesFittees[temp$annee], temp$valeur/objReport$qRel[i]/1000/1000, col=i+1, lwd=0.5)
    graphics::points(donnee$anneesFittees[temp$annee], temp$valeur/objReport$qRel[i]/1000/1000, pch=21, bg=i+1)
  }
  graphics::legend('topleft', inset=0.03, legend=paste0('qRel', 1:length(objReport$qRel), '= ',
                                              c(signif(objReport$qRel, 3))),
         lty=rep(1,length(objReport$qRel)), pch=rep(21,length(objReport$qRel)) , pt.bg=1:length(objReport$qRel) + 1)
  if(length(pl)>0){
    graphics::lines(c(donnee$anneesFittees[1]-1,donnee$anneesFittees), exp(pl$logBpred+2*plsd$logBpred)/1000/1000, lty=2, col=2)
    graphics::lines(c(donnee$anneesFittees[1]-1,donnee$anneesFittees), exp(pl$logBpred-2*plsd$logBpred)/1000/1000, lty=2, col=2)
  }
  ##
  ## projection
  if(nbAnProj > 0){
    ## print('projection a valider')
    ## exit()
    for(i in seq_along(tacProj)){
      graphics::points(max(donnee$anneesFittees)+seq(1,by=1,length.out=length(valProj[[i]]$Bproj)), valProj[[i]]$Bproj/1000/1000, pch=i, col=1)
      ## val <- projeter(objReport, c(6500,6500), fProj)
      ## graphics::points(max(donnee$anneesFittees)+seq(1,by=1,length.out=length(valProj[[i]]$Bproj)), valProj[[i]]$Bproj/1000/1000, pch=16, col=1)
      graphics::lines(max(donnee$anneesFittees)+seq(0,by=1,length.out=length(valProj[[i]]$Bproj)+1), c(utils::tail(objReport$Bpred,1), valProj[[i]]$Bproj)/1000/1000,
            lty=2)
    }
  }
  ## graphics::abline(v=2010-0.5)
  graphics::text(x=mean(graphics::graphics::par('usr')[c(1,2)]), y=diff(graphics::graphics::par('usr')[c(3,4)])*0.9 + graphics::graphics::par('usr')[3], label='A', cex=1.5)
  ##
  ##
  if(residus){
    ## B residus log
    plot(donnee$anneesFittees[subset(donnee$Bobs, source==1)$annee]-0.4,
         log(subset(donnee$Bobs, source==1)$valeur/objReport$qRel[1])-
           (log(objReport$Bpred)[subset(donnee$Bobs, source==1)$annee+1]),
         type='n', xlim=range(donnee$anneesFittees)+c(0,nbAnProj), ylim=ylimLog, xlab=labAn, ylab=labRes)
    Bres.ylim <- graphics::graphics::par('usr')[3:4]
    graphics::abline(h=0, col='grey')
    for(i in unique(donnee$Bobs$source)){
      temp <- subset(donnee$Bobs, source==i)
      for(j in 1:nrow(temp)){
        annee <- temp[j,'annee']
        graphics::lines(rep(donnee$anneesFittees[annee]+c(-0.4,0,0.4)[i], 2),
              ## c(0,log(temp[j,'valeur']/objReport$qRel[i])-log(objReport$Bpred[annee])),
              c(0,log(temp[j,'valeur']/objReport$qRel[i])-log(objReport$Bpred[annee+1]))/exp(fit$par[which(names(fit$par)=='logSigma_Bobs')][i]),
              lwd=2, col=i+1)
        graphics::points(donnee$anneesFittees[annee]+c(-0.4,0,0.4)[i],
               ## (log(temp[j,'valeur']/objReport$qRel[i])-log(objReport$Bpred[annee])),
               (log(temp[j,'valeur']/objReport$qRel[i])-log(objReport$Bpred[annee+1]))/exp(fit$par[which(names(fit$par)=='logSigma_Bobs')][i]),
               pch=21, bg=i+1)
      }
    }
    ## graphics::abline(v=2010-0.5-0.4)
    graphics::text(x=mean(graphics::graphics::par('usr')[c(1,2)]), y=diff(graphics::graphics::par('usr')[c(3,4)])*0.9 + graphics::graphics::par('usr')[3], label='B', cex=1.5)
    ##
    ## Erreur de processus sur B
    ## pour am\U{00E9}liorer, consulter la formation Halifax2024 par Anders Nielsen
    print('Erreur de processus sur B en developpement')
    ## plot(donnee$anneesFittees, (log(objReport$Bpred.proc)-log(utils::tail(objReport$Bpred,-1))), xlim=range(donnee$anneesFittees)+c(0,nbAnProj), ylim=Bres.ylim, yaxs='i',
    ##      pch=21, bg=7, xlab=labAn, ylab=labRes); graphics::abline(h=0)
    plot(donnee$anneesFittees, (log(utils::tail(objReport$Bpred,-1))-log(objReport$Bpred.proc))/exp(fit$par[which(names(fit$par)=='logSigma_Bproc')]), xlim=range(donnee$anneesFittees)+c(0,nbAnProj), ylim=Bres.ylim, yaxs='i',
         pch=21, bg=7, xlab=labAn, ylab=labRes); graphics::abline(h=0)
    for(i in 1:length(donnee$anneesFittees)){
      ## graphics::lines(rep(donnee$anneesFittees[i],2), c(0,log(objReport$Bpred.proc[i])-log(objReport$Bpred[i+1])), lwd=2, col=7)
      ## graphics::points(donnee$anneesFittees[i], (log(objReport$Bpred.proc[i])-log(objReport$Bpred[i+1])), pch=21, bg=7)
      graphics::lines(rep(donnee$anneesFittees[i],2), c(0,log(objReport$Bpred[i+1])-log(objReport$Bpred.proc[i]))/exp(fit$par[which(names(fit$par)=='logSigma_Bproc')]), lwd=2, col=7)
      graphics::points(donnee$anneesFittees[i], (log(objReport$Bpred[i+1])-log(objReport$Bpred.proc[i]))/exp(fit$par[which(names(fit$par)=='logSigma_Bproc')]), pch=21, bg=7)
    }
    graphics::text(x=mean(graphics::graphics::par('usr')[c(1,2)]), y=diff(graphics::graphics::par('usr')[c(3,4)])*0.9 + graphics::graphics::par('usr')[3], label='C', cex=1.5)
    ## graphics::axis(3, at=donnee$anneesFittees, labels=round((objReport$Bpred.proc-utils::tail(objReport$Bpred,-1))/utils::tail(objReport$Bpred,-1)*100), cex.axis=0.6)
  }
}
##
graph.R <- function(donnee, param, objReport=NULL, msyVal=NULL, tacProj=NA, fProj=NA, valProj=NULL, ylimLog=NULL, pl=list(), plsd=list(), langue=c('fr','en','bil'), residus=TRUE){
  switch(langue,
         'fr' = {labAn <- 'Ann\U{00E9}e'; labBiom <- "Biomasse ('000t)"; labRes <- 'R\U{00E9}sidus standardis\U{00E9}s'},
         'en' = {labAn <- 'Year'; labBiom <- "Biomass ('000t)"; labRes <- 'Standardized residuals'},
         'bil' = {labAn <- 'Ann\U{00E9}e/Year'; labBiom <- "Biomasse/Biomass ('000t)"; labRes <- 'R\U{00E9}sidus/Residuals'})
  if(!any(is.na(tacProj)) | !any(is.na(fProj))){
    nbAnProj <- max(length(tacProj), length(fProj), na.rm=TRUE)
  }else{
    nbAnProj <- 0
  }
  plot(donnee$anneesFittees, objReport$Rpred/1000, type='l',
       xlim=range(donnee$anneesFittees)+c(0,nbAnProj)+c(0,max(donneeInit$Robs$annee)-max(donnee$anneesFittees)),
       ylim=c(0,max(objReport$Rpred,donnee$Robs[,'valeur']/objReport$qRecru[donnee$Robs[,'source']], na.rm=TRUE))/1000,
       ylab="Nombre ('000)", xlab='Ann\U{00E9}e')
  ## graphics::abline(v=1992.5, col='grey50')
  ## graphics::abline(v=1992.5+9, col='grey50', lty=3)
  graphics::abline(h=0, col='grey')
  for(i in unique(donnee$Robs[,'source'])){
    temp <- subset(donnee$Robs, source==i)
    graphics::lines(donnee$anneesFittees[temp[,'annee']], temp[,'valeur']/objReport$qRecru[i]/1000, type='o', col=i+1, lty=3)
    temp2 <- subset(donneeInit$Robs, source==i&annee>=donnee$anneesFittees[max(temp$annee)])
    if(nrow(temp2)>0){
      graphics::lines(temp2[,'annee'], temp2[,'valeur']/objReport$qRecru[i]/1000, type='o', col=i+1, lty=3)
    }
  }
  graphics::lines(donnee$anneesFittees, objReport$Rpred/1000)
  graphics::legend('topleft', inset=0.03, legend=paste('q=', round(objReport$qRecru, 3)), pch=1, col=unique(donnee$Robs[,'source'])+1)
  nbAnMoy <- 15
  graphics::lines(utils::tail(donnee$anneesFittees,nbAnMoy), rep(mean(utils::tail(objReport$Rpred/1000,nbAnMoy)),nbAnMoy), col=2, lwd=2)
  graphics::axis(4, at=mean(utils::tail(objReport$Rpred/1000,nbAnMoy)), label=round(mean(utils::tail(objReport$Rpred/1000,nbAnMoy))))
  if(nbAnProj > 0){
    for(i in seq_along(tacProj)){
      graphics::points(max(donnee$anneesFittees)+seq(1,by=1,length.out=length(valProj[[i]]$Rproj)), valProj[[i]]$Rproj/1000, pch=i, col=1)
      graphics::lines(max(donnee$anneesFittees)+seq(0,by=1,length.out=length(valProj[[i]]$Rproj)+1), c(utils::tail(objReport$Rpred,1), valProj[[i]]$Rproj)/1000,
            lty=2)
    }
  }
  graphics::text(x=mean(graphics::graphics::par('usr')[c(1,2)]), y=diff(graphics::graphics::par('usr')[c(3,4)])*0.9 + graphics::graphics::par('usr')[3], label='A', cex=1.5)
  ##
  if(length(pl)>0){
    graphics::lines(donnee$anneesFittees, exp(pl$logRpred+2*plsd$logRpred)/1000, lty=2, col=2)
    graphics::lines(donnee$anneesFittees, exp(pl$logRpred-2*plsd$logRpred)/1000, lty=2, col=2)
  }
  ##
  ## residus
  if(residus){
    plot(donnee$anneesFittees[subset(donnee$Robs, source==1)$annee]-0.4,
         (log(subset(donnee$Robs, source==1)$valeur / objReport$qRecru[1]) - (log(objReport$Rpred)[subset(donnee$Robs, source==1)$annee])) / exp(fit$par[which(names(fit$par)=='logSigma_Robs')][1]),
         xlim=range(donnee$anneesFittees)+c(0,nbAnProj), ylim=ylimLog, type='n', xlab='Ann\U{00E9}e', ylab='R\U{00E9}sidus standardis\U{00E9}s')
    graphics::abline(h=0, col='grey')
    for(i in unique(donnee$Robs$source)){
      temp <- subset(donnee$Robs, source==i)
      for(j in 1:nrow(temp)){
        annee <- temp[j,'annee']
        sigma <- temp[j,'sigma']
        graphics::lines(rep(donnee$anneesFittees[annee]+c(-0.4,0,0.4)[i], 2),
              c(0,log(temp[j,'valeur']/objReport$qRecru[i])-log(objReport$Rpred[annee]))/exp(fit$par[which(names(fit$par)=='logSigma_Robs')][sigma]), lwd=2, col=i+1)
        graphics::points(donnee$anneesFittees[annee]+c(-0.4,0,0.4)[i], (log(temp[j,'valeur']/objReport$qRecru[i])-log(objReport$Rpred[annee]))/exp(fit$par[which(names(fit$par)=='logSigma_Robs')][sigma]),
               pch=21, bg=i+1)
      }
    }
    graphics::text(x=mean(graphics::graphics::par('usr')[c(1,2)]), y=diff(graphics::graphics::par('usr')[c(3,4)])*0.9 + graphics::graphics::par('usr')[3], label='B', cex=1.5)
    ##
    ## proportion de recures R/N
    plot(donnee$anneesFittees, objReport$Rpred/objReport$Npred, type='l', xlim=range(donnee$anneesFittees)+c(0,nbAnProj), ylim=c(0,0.5),
         xlab=labAn, ylab='Proportion')
    graphics::lines(donnee$anneesFittees, objReport$Rpred*donnee$omegaK$valeur/utils::tail(objReport$Bpred,-1), col=2)
    graphics::abline(h=0, col='grey')
    graphics::legend('topright', inset=0.03, legend=c('R/N','R/B'), col=c(1,2), lty=1)
    graphics::text(x=mean(graphics::graphics::par('usr')[c(1,2)]), y=diff(graphics::graphics::par('usr')[c(3,4)])*0.9 + graphics::graphics::par('usr')[3], label='C', cex=1.5)
    if(nbAnProj > 0){
      for(i in seq_along(tacProj)){
        graphics::points(max(donnee$anneesFittees)+seq(1,by=1,length.out=length(valProj[[i]]$Rproj)), valProj[[i]]$Rproj/valProj[[i]]$Nproj, pch=i, col=1)
        graphics::lines(max(donnee$anneesFittees)+seq(0,by=1,length.out=length(valProj[[i]]$Rproj)+1), c(utils::tail(objReport$Rpred,1), valProj[[i]]$Rproj)/c(utils::tail(objReport$Npred,1), valProj[[i]]$Nproj),
              lty=2)
        graphics::points(max(donnee$anneesFittees)+seq(1,by=1,length.out=length(valProj[[i]]$Rproj)), valProj[[i]]$Rproj*utils::tail(donnee$omegaK$valeur,1)/valProj[[i]]$Bproj, pch=i, col=2)
        graphics::lines(max(donnee$anneesFittees)+seq(0,by=1,length.out=length(valProj[[i]]$Rproj)+1), c(utils::tail(objReport$Rpred,1), valProj[[i]]$Rproj)*utils::tail(donnee$omegaK$valeur,1)/c(utils::tail(objReport$Bpred,1), valProj[[i]]$Bproj),
              lty=2, col=2)
      }
    }
  }
}
graph.kobe <- function(donnee, param, objReport=NULL, msyVal=NULL, tacProj=NA, fProj=NA, valProj=NULL, Bmax=NA, langue=c('fr','en','bil')){
  switch(langue,
         'fr' = {labAn <- 'Ann\U{00E9}e'; labBiom <- "Biomasse ('000t)"; labRes <- 'R\U{00E9}sidus standardis\U{00E9}s'},
         'en' = {labAn <- 'Year'; labBiom <- "Biomass ('000t)"; labRes <- 'Standardized residuals'},
         'bil' = {labAn <- 'Ann\U{00E9}e/Year'; labBiom <- "Biomasse/Biomass ('000t)"; labRes <- 'R\U{00E9}sidus/Residuals'})
  if(!any(is.na(tacProj)) | !any(is.na(fProj))){
    nbAnProj <- max(length(tacProj[[1]]), length(fProj), na.rm=TRUE)
  }else{
    nbAnProj <- 0
  }
  if(is.null(msyVal)){x.max <- max(objReport$Bpred,Bmax)/1000/1000}else{x.max <- max(objReport$Bpred, msyVal$B, valProj[[1]]$Bproj, Bmax, na.rm=TRUE)/1000/1000}
  if(is.null(msyVal)){y.max <- max(objReport$F)}else{y.max <- min(max(c(objReport$F,msyVal$F),na.rm=TRUE),2*msyVal$F)}
  plot(utils::tail(objReport$Bpred,-1)/1000/1000, objReport$F, type='o',
       xlim=c(0,x.max*1.1), ylim=c(0,y.max*1.1), xaxs='i', yaxs='i',
       xlab="B85+ ('000 t)", ylab="F"); graphics::abline(h=0, col='grey70'); graphics::abline(v=0, col='grey70')
  if(nbAnProj > 0){
    for(i in seq_along(tacProj)){
      graphics::points(valProj[[i]]$Bproj/1000/1000, valProj[[i]]$fProj, pch=i, col=1)
      graphics::lines(c(utils::tail(objReport$Bpred,1), valProj[[i]]$Bproj)/1000/1000, c(-log(1-utils::tail(objReport$tauxExp,1)), valProj[[i]]$fProj),
            lty=2)
    }
  }
  lesquels <- which(donnee$anneesFittees %in% pretty(donnee$anneesFittees))
  graphics::text((utils::tail(objReport$Bpred,-1)/1000/1000)[lesquels], (objReport$F)[lesquels], labels=donnee$anneesFittees[lesquels],
       pos=3)
  if(!is.null(msyVal)){
    graphics::abline(v=msyVal$B/1000/1000, lty=2)
    graphics::abline(v=as.numeric(msyVal$B)/1000/1000*c(0.4,0.8), lty=3, col='grey70')
    graphics::abline(h=msyVal$F, col=2)
    graphics::axis(4, at=msyVal$F, labels='Fmsy', col.axis=2, las=1)
    graphics::axis(3, at=msyVal$B/1000/1000*c(0.4,0.8,1), labels=c('0.4','0.8','Bmsy'), las=1)
  }
  graphics::legend('topright', inset=0.03, legend=paste0(c('B','F'), utils::tail(donnee$anneesFittees,1), c('/Bmsy = ','/Fmsy = '),
                                               c(round(utils::tail(objReport$Bpred,1)/msyVal$B,2),round(utils::tail(objReport$F,1)/msyVal$F,2))))
}
##
graph.omega <- function(donnee, param, objReport=NULL, msyVal=NULL, tacProj=NA, fProj=NA, valProj=NULL, ylimLog=NULL, langue=c('fr','en','bil'), residus=TRUE){
  switch(langue,
         'fr' = {labAn <- 'Ann\U{00E9}e'; labPoids <- "Poids individuel moyen (kg)"; labRes <- 'R\U{00E9}sidus standardis\U{00E9}'},
         'en' = {labAn <- 'Year'; labPoids <- "Mean individual weight (kg)"; labRes <- 'Strandardized residuals'},
         'bil' = {labAn <- 'Ann\U{00E9}e/Year'; labPoids <- "Poids individuel moyen/Mean indivudual weight (kg)"; labRes <- 'R\U{00E9}sidus/Residuals'})
  if(!any(is.na(tacProj)) | !any(is.na(fProj))){
    nbAnProj <- max(length(tacProj[[1]]), length(fProj), na.rm=TRUE)
  }else{
    nbAnProj <- 0
  }
  ##
  plot(donnee$anneesFittees, objReport$omegaPred, type='l', xlim=range(donnee$anneesFittees)+c(0,nbAnProj), ylim=c(7,max(donnee$omega$valeur)), xlab=labAn, ylab=labPoids)
  graphics::abline(h=0, col='grey70')
  for(i in unique(donnee$omega$source)){
    temp <- subset(donnee$omega, source==i)
    graphics::points(donnee$anneesFittees[temp$annee], temp$valeur, pch=16, col=i+1)
  }
  graphics::axis(4, at=as.numeric(donnee$lpAlpha)*pretty((graphics::graphics::par('usr')[3:4]/as.numeric(donnee$lpAlpha)) ^ (1/as.numeric(donnee$lpBeta)))^as.numeric(donnee$lpBeta),
       labels=pretty((graphics::graphics::par('usr')[3:4]/as.numeric(donnee$lpAlpha)) ^ (1/as.numeric(donnee$lpBeta))))
  if(nbAnProj > 0){
    for(i in seq_along(tacProj)){
      graphics::points(max(donnee$anneesFittees)+seq(1,by=1,length.out=length(valProj[[i]]$Rproj)), valProj[[i]]$Bproj/valProj[[i]]$Nproj, pch=i, col=1)
      graphics::lines(max(donnee$anneesFittees)+seq(0,by=1,length.out=length(valProj[[i]]$Rproj)+1), c(utils::tail(objReport$omegaPred,1), valProj[[i]]$Bproj/valProj[[i]]$Nproj),
            lty=2)
    }
  }
  graphics::legend('topleft', inset=0.03, legend=c('\U{00E9}chantillonneur \U{00E0} quai','Observateur en mer','Mod\U{00E9}lis\U{00E9}'), lty=c(NA,NA,1), pch=c(16,16,NA), col=c(2,3,1))
  graphics::text(x=mean(graphics::graphics::par('usr')[c(1,2)]), y=diff(graphics::graphics::par('usr')[c(3,4)])*0.9 + graphics::graphics::par('usr')[3], label='A', cex=1.5)
  ##
  ## residuels
  if(residus){
    if(is.null(ylimLog)){
      ylimLog <- range(c(log(subset(donnee$omega, source==1)$valeur)-log(objReport$omegaPred[subset(donnee$omega, source==1)$annee]),
                         log(subset(donnee$omega, source==2)$valeur)-log(objReport$omegaPred[subset(donnee$omega, source==2)$annee])))
    }
    plot(donnee$anneesFittees[subset(donnee$omega, source==1)$annee]-0.2,
         (subset(donnee$omega, source==1)$valeur-objReport$omegaPred[subset(donnee$omega, source==1)$annee])/exp(fit$par[which(names(fit$par)=='logSigma_oBar')][1]),
         ## (log(subset(donnee$omega, source==1)$valeur)-log(objReport$omegaPred[subset(donnee$omega, source==1)$annee])),
         xlim=range(donnee$anneesFittees)+c(0,nbAnProj),
         ylim=ylimLog,
         type='n', xlab=labAn, ylab=labRes)
    graphics::abline(h=0, col='grey')
    for(i in unique(donnee$omega$source)){
      temp <- subset(donnee$omega, source==i)
      for(j in 1:nrow(temp)){
        annee <- temp[j,'annee']
        sigma <- temp[j,'sigma']
        graphics::lines(rep(donnee$anneesFittees[annee]+c(-0.2,0,0.2)[i], 2), (c(0,temp[j,'valeur']-objReport$omegaPred[annee]))/exp(fit$par[which(names(fit$par)=='logSigma_oBar')][sigma]),
              lwd=2, col=i+1)
        ## graphics::lines(rep(donnee$anneesFittees[annee]+c(-0.2,0,0.2)[i], 2), (c(0,log(temp[j,'valeur'])-log(objReport$omegaPred[annee]))),
        ##       lwd=2, col=i+1)
        graphics::points(donnee$anneesFittees[annee]+c(-0.2,0,0.2)[i], (temp[j,'valeur']-objReport$omegaPred[annee])/exp(fit$par[which(names(fit$par)=='logSigma_oBar')][sigma]), pch=16, col=i+1)
        ## graphics::points(donnee$anneesFittees[annee]+c(-0.2,0,0.2)[i], (log(temp[j,'valeur'])-log(objReport$omegaPred[annee])), pch=16, col=i+1)
      }
    }
    graphics::text(x=mean(graphics::graphics::par('usr')[c(1,2)]), y=diff(graphics::graphics::par('usr')[c(3,4)])*0.9 + graphics::graphics::par('usr')[3], label='B', cex=1.5)
  }
}
##
graph.N <- function(donnee, param, objReport=NULL, msyVal=NULL, tacProj=NA, valProj=NULL, fProj=NA, langue=c('fr','en','bil')){
  switch(langue,
         'fr' = {labAn <- 'Ann\U{00E9}e'; labNb <- "Nombre"},
         'en' = {labAn <- 'Year'; labNb <- "Number"},
         'bil' = {labAn <- 'Ann\U{00E9}e/Year'; labNb <- "Nombre/Number"})
  if(!any(is.na(tacProj)) | !any(is.na(fProj))){
    nbAnProj <- max(length(tacProj[[1]]), length(fProj), na.rm=TRUE)
  }else{
    nbAnProj <- 0
  }
  ##
  plot(donnee$anneesFittees,objReport$Npred,type='l', range(donnee$anneesFittees)+c(0,nbAnProj), ylim=c(0,max(objReport$Npred)),
       xlab=labAn, ylab=labNb)
  graphics::abline(h=0, col='grey70')
  graphics::lines(donnee$anneesFittees,objReport$Rpred, col=2)
  graphics::legend('topleft', inset=0.03, legend=c('N total', 'Recrutement'), lty=1, col=c(1,2))
  if(nbAnProj > 0){
    for(i in seq_along(tacProj)){
      graphics::points(max(donnee$anneesFittees)+seq(1,by=1,length.out=length(valProj[[i]]$Nproj)), valProj[[i]]$Nproj, pch=i, col=1)
      graphics::lines(max(donnee$anneesFittees)+seq(0,by=1,length.out=length(valProj[[i]]$Nproj)+1), c(utils::tail(objReport$Npred,1), valProj[[i]]$Nproj),
            lty=2)
      graphics::points(max(donnee$anneesFittees)+seq(1,by=1,length.out=length(valProj[[i]]$Rproj)), valProj[[i]]$Rproj, pch=i, col=2)
      graphics::lines(max(donnee$anneesFittees)+seq(0,by=1,length.out=length(valProj[[i]]$Rproj)+1), c(utils::tail(objReport$Rpred,1), valProj[[i]]$Rproj),
            lty=2, col=2)
    }
  }
}
##
graph.C <- function(donnee, param, objReport=NULL, msyVal=NULL, tacProj=NA, fProj=NA, valProj=NULL, Cmax=NA, langue=c('fr','en','bil')){
  switch(langue,
         'fr' = {labAn <- 'Ann\U{00E9}e'; labDeb <- "D\U{00E9}barquements ('000t)"; labF <- 'Taux instantan\U{00E9}'},
         'en' = {labAn <- 'Year'; labDeb <- "Landings ('000t)"; labF <- 'Instantaneous rate'},
         'bil' = {labAn <- 'Ann\U{00E9}e/Year'; labDeb <- "D\U{00E9}barquements/Landings ('000t)"; labF <- 'Taux instantan\U{00E9}/Instantaneous rate'})
  if(!any(is.na(tacProj)) | !any(is.na(fProj))){
    nbAnProj <- max(length(tacProj[[1]]), length(fProj), na.rm=TRUE)
  }else{
    nbAnProj <- 0
  }
  ##
  ## if(nbAnProj > 0){
  ##     val <- projeter(objReport, tacProj, fProj)
  ## }else{
  ##     val <- list()
  ##     val$Cproj <- NA
  ## }
  plot(donnee$anneesFittees, objReport$Cpred/1000, type='l', xlim=range(donnee$anneesFittees)+c(0,nbAnProj), ylim=c(0,max(2500,Cmax,na.rm=TRUE)),
       xlab=labAn, ylab=labDeb)
  graphics::abline(h=0, col='grey70')
  graphics::lines(donnee$anneesFittees[donnee$Cobs[,'annee']], donnee$Cobs[,'valeur']/1000, col='pink')
  for(i in seq_along(donnee$Cobs[,'valeur'])){
    graphics::lines(rep(donnee$anneesFittees[donnee$Cobs[i,'annee']], 2), c(donnee$Cobs[i,'valeur']/1000, objReport$Cpred[i]/1000), lty=1, col=2, lwd=2)
  }
  graphics::lines(donnee$anneesFittees, objReport$Cpred/1000)
  if(!is.null(msyVal)){
    graphics::abline(h=msyVal$C/1000, col=2)
    graphics::axis(4, at=msyVal$C/1000, labels='msy', col.axis=2, las=1)
  }
  if(nbAnProj > 0){
    for(i in seq_along(tacProj)){
      graphics::points(max(donnee$anneesFittees)+seq(1,by=1,length.out=length(valProj[[i]]$Bproj)), tacProj[[i]], pch=i, col=1)
      graphics::lines(max(donnee$anneesFittees)+seq(0,by=1,length.out=length(valProj[[i]]$Bproj)+1), c(utils::tail(objReport$Cpred,1)/1000, tacProj[[i]]),
            lty=2)
    }
  }
  if(!is.null(msyVal)){
    graphics::legend('topleft', legend=paste(c('msy='), c(round(msyVal$C/1000))))
  }
  graphics::text(x=mean(graphics::graphics::par('usr')[c(1,2)]), y=diff(graphics::graphics::par('usr')[c(3,4)])*0.9 + graphics::graphics::par('usr')[3], label='A', cex=1.5)
  ##
}
##
graph.F <- function(donnee, param, objReport=NULL, msyVal=NULL, tacProj=NA, fProj=NA, valProj=NULL, Fmax=NA, pl=list(), plsd=list(), langue=c('fr','en','bil')){
  switch(langue,
         'fr' = {labAn <- 'Ann\U{00E9}e'; labDeb <- "D\U{00E9}barquements ('000t)"; labF <- 'Taux instantan\U{00E9}'},
         'en' = {labAn <- 'Year'; labDeb <- "Landings ('000t)"; labF <- 'Instantaneous rate'},
         'bil' = {labAn <- 'Ann\U{00E9}e/Year'; labDeb <- "D\U{00E9}barquements/Landings ('000t)"; labF <- 'Taux instantan\U{00E9}/Instantaneous rate'})
  if(!any(is.na(tacProj)) | !any(is.na(fProj))){
    nbAnProj <- max(length(tacProj[[1]]), length(fProj), na.rm=TRUE)
  }else{
    nbAnProj <- 0
  }
  ##
  plot(donnee$anneesFittees, objReport$F, type='l', xlim=range(donnee$anneesFittees)+c(0,nbAnProj), ylim=c(0,max(0.3,Fmax,na.rm=TRUE)),
       xlab=labAn, ylab=labF)
  graphics::abline(h=0, col='grey70')
  graphics::abline(h=donnee$M, col=4, lwd=2)
  graphics::axis(4, at=donnee$M, labels='M', col.axis=4, las=1)
  graphics::lines(donnee$anneesFittees[length(objReport$F)+c(-9,0)], rep(mean(utils::tail(objReport$F,10)),2), col=3, lwd=3)
  for(i in 1){
    temp <- subset(donnee$Bobs, source==i)
    graphics::points(donnee$anneesFittees[temp$annee],
           ## 1-exp(-donnee$Cobs[temp$annee,'valeur']/(temp$valeur/objReport$qRel[i])), pch=21, bg=i+1)
           1-exp(-donnee$Cobs[temp$annee,'valeur']/(temp$valeur/objReport$qRel[i]*exp(-donnee$M))), pch=21, bg=i+1)
  }
  ## if(!is.null(objReport$nTagCapt)) graphics::points(donnee$anneesFittees, 1-exp(-(objReport$NtagCapt/objReport$tauxRetour)/c(0,utils::tail(objReport$Ntag,-1))), pch=21, bg=5)
  nbAn <- 10
  if(!is.null(msyVal)){
    graphics::legend('topright', inset=0.03,
           legend=paste(c('M=','F 10ans','Fmsy='),
                        c(round(donnee$M,3),' ',round(msyVal$F,3))),
           lty=c(1,1,1,NA,NA), col=c(4,3,2,NA,NA), lwd=c(2,2,3,NA,NA))
    graphics::abline(h=msyVal$F, col=2, lwd=2)
    graphics::axis(4, at=msyVal$F, labels='Fmsy', col.axis=2, las=1)
  }else{
    graphics::legend('topright', inset=0.03,
           legend=paste(c('M=','F 10ans'), c(round(donnee$M,3),' ')),
           lty=c(1,1,NA,NA), col=c(4,3,NA,NA), lwd=c(2,3,NA,NA))
  }
  ## graphics::legend('topleft', inset=0.03, legend=c('Deb/Releve','tags'), pch=rep(16,2), col=c(2,5))
  graphics::axis(4, at=mean(utils::tail(objReport$F,10)), label=round(mean(utils::tail(objReport$F,10)),4), col.lab=3)
  if(length(pl)>0){
    graphics::lines(donnee$anneesFittees, 0.001 + 0.9*stats::plogis(pl$transTauxExp+2*plsd$transTauxExp), lty=2, col=2)
    graphics::lines(donnee$anneesFittees, 0.001 + 0.9*stats::plogis(pl$transTauxExp-2*plsd$transTauxExp), lty=2, col=2)
  }
  if(nbAnProj > 0){
    for(i in seq_along(tacProj)){
      graphics::points(max(donnee$anneesFittees)+seq(1,by=1,length.out=length(valProj[[i]]$fProj)), valProj[[i]]$fProj, pch=i, col=1)
      graphics::lines(max(donnee$anneesFittees)+seq(0,by=1,length.out=length(valProj[[i]]$fProj)+1), c(utils::tail(objReport$F,1), valProj[[i]]$fProj),
            lty=2)
    }
  }
  graphics::lines(donnee$anneesFittees, objReport$F)
  graphics::text(x=mean(graphics::graphics::par('usr')[c(1,2)]), y=diff(graphics::graphics::par('usr')[c(3,4)])*0.9 + graphics::graphics::par('usr')[3], label='B', cex=1.5)
}
##
graph.retourTag <- function(donnee, param, objReport=NULL, msyVal=NULL, tacProj=NA, fProj=NA, valProj=NULL, langue=c('fr','en','bil'), residus=TRUE){
  switch(langue,
         'fr' = {labAn <- 'Ann\U{00E9}e'; labNbRet <- 'Nombre de retours'; labTxRet <- 'Taux de retour'; labSurvie <- 'Survie post-marquage';
         labAnCapt='Ann\U{00E9}e de capture'; labAnMarq='Ann\U{00E9}es de marquage'},
         'en' = {labAn <- 'Year'; labNbRet <- 'Number of returns'; labTxRet <- 'Return rate'; labSurvie <- 'Post-tagging survival';
         labAnCapt='Year of capture'; labAnMarq='Year of tagging'},
         'bil' = {labAn <- 'Ann\U{00E9}e/Year'; labNbRet <- 'Nombre de retours/Number of returns'; labTxRet <- 'Taux de retour/Return rate';
         labSurvie <- 'Survie post-marquage/Post-tagging survival';
         labAnCapt='Capture'; labAnMarq='Marquage/Tagging'})
  if(!any(is.na(tacProj)) | !any(is.na(fProj))){
    nbAnProj <- max(length(tacProj[[1]]), length(fProj), na.rm=TRUE)
  }else{
    nbAnProj <- 0
  }
  ##
  ## temp <- donnee$nTagsRetour; temp[which(temp==0)] <- NA
  ## plot(donnee$anneesFittees, objReport$nTagRetourPred,
  ##      xlim=range(donnee$anneesFittees)+c(0,nbAnProj), ylim=c(0,max(c(objReport$nTagRetour,temp),na.rm=TRUE)),type='l',
  ##      xlab=labAn, ylab=labNbRet)#; graphics::abline(h=0, col='grey70')
  ## graphics::points(donnee$anneesFittees, temp, pch=21, bg=3)
  ## for(i in seq_along(temp)){
  ##     if(!is.na(temp[i])){
  ##         graphics::lines(rep(donnee$anneesFittees[i],2), c(temp[i],objReport$NtagCapt[i]), lty=1, col=2)
  ##     }
  ## }
  ## graphics::legend('topleft', inset=0.03,
  ##        legend=paste(c('tx ret.=','s. p-m='),
  ##                     c(round(objReport$tauxRetour,2), round(donnee$sPostMarquage,2))))
  ##
  plot(0,0,type='n',
       xlim=c(utils::head(which(donnee$nTagsPoses$valeur>0),1)+1,length(donnee$anneesFittees))+c(-0.5,nbAnProj+0.5), ylim=c(0,max(donnee$nTagsRetourObs$valeur,objReport$nTagRetourPred,na.rm=TRUE)),
       xlab=labAn, ylab=labNbRet, axes=FALSE)
  graphics::box(); graphics::axis(1, at=min(donnee$anneesFitteesID):(max(donnee$anneesFitteesID)+nbAnProj), labels=min(donnee$anneesFittees):(max(donnee$anneesFittees)+nbAnProj))
  graphics::axis(2); graphics::abline(h=0, col='grey70')
  quelle.annee <- subset(donnee$nTagsPoses, valeur>0)$annee
  for(i in seq_along(quelle.annee)){
    graphics::lines(1:ncol(objReport$nTagRetourPred)+rep(seq(-0.1,0.1,length.out=length(unique(donnee$nTagsRetourObs$anneeRecap))),100)[i],
          objReport$nTagRetourPred[quelle.annee[i],], type='o', col=i+1, lwd=2, pch=16)
    obs.temp <- subset(donnee$nTagsRetourObs, anneePose==quelle.annee[i])
    graphics::points(obs.temp$anneeRecap+rep(seq(-0.1,0.1,length.out=length(unique(donnee$nTagsRetourObs$anneeRecap))),100)[i], obs.temp[,'valeur'], pch=1, col=i+1)#pch=22, bg=i+1)
    for(j in 1:nrow(obs.temp)){
      graphics::lines(rep(obs.temp$anneeRecap[j]+rep(seq(-0.1,0.1,length.out=length(unique(donnee$nTagsRetourObs$anneeRecap))),100)[i],2),
            c(objReport$nTagRetourPred[quelle.annee[i],obs.temp[j,'anneeRecap']], obs.temp[j,'valeur']), col=i+1)
    }
  }
  ## graphics::axis(3, at=1:length(objReport$F), labels=round(objReport$F,3), tick=FALSE, line=-1, cex.axis=0.7)
  graphics::axis(3, at=1:ncol(objReport$nTagRetourPred), labels=round(apply(objReport$nTagRetourPred,2,sum,na.rm=TRUE),1), tick=FALSE, line=-1, cex.axis=0.7)
  temp <- stats::aggregate(donnee$nTagsRetourObs$valeur, by=donnee$nTagsRetourObs['anneeRecap'], FUN=sum)
  graphics::axis(3, at=temp$anneeRecap, labels=temp$x, tick=FALSE, line=-2, cex.axis=0.7)
  ## graphics::legend('topleft', inset=0.03,
  ##        legend=paste(c(labTxRet,labSurvie),
  ##                     c(round(objReport$tauxRetour,2), round(donnee$sPostMarquage,2))))
  if(nbAnProj > 0){
    for(i in seq_along(tacProj)){
      for(i in seq_along(quelle.annee)){
        graphics::points(max(donnee$anneesFitteesID)+seq(1,by=1,length.out=ncol(valProj[[i]]$NtagCaptProj)), valProj[[i]]$NtagCaptProj[quelle.annee[i],], col=i+1, lwd=2, pch=i)
        graphics::lines(max(donnee$anneesFitteesID)+seq(0,by=1,length.out=ncol(valProj[[i]]$NtagCaptProj)+1), c(objReport$nTagRetourPred[quelle.annee[i],ncol(objReport$nTagRetourPred)], valProj[[i]]$NtagCaptProj[quelle.annee[i],]),
              lty=2, col=i+1)
      }
      graphics::axis(3, at=max(donnee$anneesFitteesID)+seq(1,by=1,length.out=ncol(valProj[[i]]$NtagCaptProj)), labels=round(apply(valProj[[i]]$NtagCaptProj,2,sum),1), tick=FALSE, line=-1, cex.axis=0.7)
    }
  }
  ##
  ## r\U{00E9}siduels
  if(residus){
    plot(0,0,type='n',
         xlim=donnee$anneesFittees[range(donnee$nTagsRetourObs$anneeRecap)], ylim=donnee$anneesFittees[range(donnee$nTagsRetourObs$anneePose)]+c(-0.5,0.5),
         xlab=labAnCapt, ylab=labAnMarq)
    for(i in 1:nrow(donnee$nTagsRetourObs)){
      temp <- donnee$nTagsRetourObs[i,]
      res <- temp$valeur-objReport$nTagRetourPred[temp$anneePose,temp$anneeRecap]
      resStd <- (temp$valeur-objReport$nTagRetourPred[temp$anneePose,temp$anneeRecap])/exp(fit$par['logSigma_retourTag'])
      graphics::points(donnee$anneesFittees[temp$anneeRecap], donnee$anneesFittees[temp$anneePose], cex=2*sqrt(abs(resStd)), bg=ifelse(res>0,2,4), pch=21)
    }
  }
}
##
graph.SSR <- function(donnee, param, objReport=NULL, msyVal=NULL, bh=NULL, ylimLog=NULL, steep=TRUE, lesquels=1:4, langue=c('fr','en','bil')){
  switch(langue,
         'fr' = {labAn <- 'Ann\U{00E9}e'; labB85 <- "B85+ ('000t)"; labRecru <- "Recrues ('000)"; labRes <- 'R\U{00E9}sidus standardis\U{00E9}s'},
         'en' = {labAn <- 'Year'; labB85 <- "B85+ ('000t)"; labRecru <- "Recruits ('000)"; labRes <- 'Standardized residuals'},
         'bil' = {labAn <- 'Ann\U{00E9}e/Year'; labB85 <- "B85+ ('000t)"; labRecru <- "Recrues/Recruits ('000)"; labRes <- 'R\U{00E9}sidus/Residuals'})
  ##
  annees <- donnee$anneesBH[1]:donnee$anneesBH[2]
  if(1%in%lesquels){
    plot(utils::head(utils::tail(objReport$Bpred,-1),-(donnee$lagBH)), utils::tail(objReport$Rpred,-(donnee$lagBH)), type='n',
         xlim=c(0,max(objReport$Bpred,msyVal$B0)),
         ylim=c(0,max(objReport$Rpred)),
         ## subset(donnee$Robs,source==1 & annee%in%(annees+donnee$lagBH))[,'valeur'] / objReport$qRecru[1],
         ## subset(donnee$Robs,source==2 & annee%in%(annees+donnee$lagBH))[,'valeur'] / objReport$qRecru[2])),
         pch=21, bg=2, xlab=labB85, ylab=labRecru, axes=FALSE); graphics::box()
    graphics::axis(1, at=pretty(graphics::graphics::par('usr')[1:2]), labels=pretty(graphics::graphics::par('usr')[1:2])/1000000)
    graphics::axis(2, at=pretty(graphics::graphics::par('usr')[3:4]), labels=pretty(graphics::graphics::par('usr')[3:4])/1000)
    graphics::rug(utils::tail(objReport$Bpred,donnee$lagBH)); graphics::rug(utils::head(objReport$Rpred,donnee$lagBH), side=2)
    ## for(i in (donnee$lagBH+1):length(objReport$Rpred)){
    ##     graphics::points(objReport$Bpred[i-donnee$lagBH], objReport$Rbh[i-donnee$lagBH])
    ##     graphics::points(objReport$Bpred[i-donnee$lagBH], objReport$Rpred[i], col=2)
    ## }
    graphics::abline(h=0, col='grey70'); graphics::abline(v=0, col='grey70')
    ## if(donnee$switchBH){
    ##     graphics::curve(objReport$alphaBH*(x)/ (objReport$betaBH + x), add=TRUE, col=1, lty=2, lwd=2)
    ##     graphics::legend('bottomright', inset=0.03, legend=paste(c('alpha','beta'), '=', round(c(objReport$alphaBH, objReport$betaBH))))
    ## }else{
    if(is.null(msyVal)){
      msyVal <- calculerBH(donnee, objReport=objReport, steep=FALSE)
    }
    if(!steep){
      alphaBH <- msyVal$alphaBH
      betaBH <- msyVal$betaBH
    }else{
      alphaBH <- bh$alphaBH.steep
      betaBH <- bh$betaBH.steep
    }
    graphics::curve(alphaBH*(x)/ (betaBH + x), add=TRUE, col=1, lty=2, lwd=2)
    ## }
    ##
    ##
    ## eq <- substitute('R'==a%*%'B85+'/(b+'B85+'),list(a = round(alphaBH), b = round(betaBH)))
    ## graphics::text(0.5*graphics::graphics::par('usr')[2],0.9*graphics::graphics::par('usr')[4],labels=eq)
    ##
    graphics::points(utils::head(objReport$Bpred[annees+1],-(donnee$lagBH)), utils::tail(objReport$Rpred[annees],-(donnee$lagBH)), type='o', pch=21, bg=2)
    graphics::abline(a=0, b=alphaBH/betaBH, col=4, lty=3)
    graphics::lines(c(0,msyVal$B0), c(0,alphaBH*(msyVal$B0)/ (betaBH + msyVal$B0)), col=4, lty=3)
    CR <- (alphaBH/betaBH)/((alphaBH*(msyVal$B0)/ (betaBH + msyVal$B0))/msyVal$B0)
    ## print(paste('CR =', round(CR)))
    ## print(paste('h =', round(CR/(4+CR),2)))
    ## graphics::abline(v=msyVal$B0*c(1,0.2), lty=3, col=4)
    ## graphics::abline(h=alphaBH*(msyVal$B0*c(1,0.2))/ (betaBH + msyVal$B0*c(1,0.2)), lty=3, col=4)
    ## graphics::text(x=mean(graphics::graphics::par('usr')[c(1,2)]), y=diff(graphics::graphics::par('usr')[c(3,4)])*0.1 + graphics::graphics::par('usr')[3], label=paste('h=', round(alphaBH*(msyVal$B0*c(1))/ (betaBH + msyVal$B0*c(1)) / alphaBH*(msyVal$B0*c(0.2))/ (betaBH + msyVal$B0*c(0.2)),2)))
    ##
    noms <- utils::head(donnee$anneesFittees[annees],-(donnee$lagBH-1))
    noms[which(!noms %in% pretty(noms))] <- NA
    graphics::text(utils::head(objReport$Bpred[annees+1],-(donnee$lagBH)), utils::tail(objReport$Rpred[annees],-(donnee$lagBH)),
         labels=noms, pos=3)
    graphics::abline(v=msyVal$B, lty=2)
    graphics::abline(v=as.numeric(msyVal$B)*c(0.4,0.8), lty=2, col=c(2,3))
    graphics::legend('bottomright', inset=0.03, legend=paste(c('alpha','beta','CR','h'), '=', c(round(c(alphaBH, betaBH)),round(c(CR,CR/(4+CR)),2))))
    graphics::text(x=mean(graphics::graphics::par('usr')[c(1,2)]), y=diff(graphics::graphics::par('usr')[c(3,4)])*0.9 + graphics::graphics::par('usr')[3], label='A', cex=1.5)
  }
  ##
  ## agrandi
  if(2%in%lesquels){
    if(is.null(msyVal)){
      msyVal <- calculerBH(donnee, objReport, steep=FALSE)
    }
    if(!steep){
      alphaBH <- msyVal$alphaBH
      betaBH <- msyVal$betaBH
    }else{
      alphaBH <- bh$alphaBH.steep
      betaBH <- bh$betaBH.steep
    }
    plot(utils::head(utils::tail(objReport$Bpred,-1),-(donnee$lagBH)), utils::tail(objReport$Rpred,-(donnee$lagBH)), type='n',
         xlim=c(0,max(objReport$Bpred[donnee$anneesBH[1]:(donnee$anneesBH[2]-donnee$lagBH+1)])),
         ylim=c(0,max(objReport$Rpred, alphaBH*max(objReport$Bpred[donnee$anneesBH[1]:(donnee$anneesBH[2]-donnee$lagBH+1)])/
                        (betaBH + max(objReport$Bpred[donnee$anneesBH[1]:(donnee$anneesBH[2]-donnee$lagBH+1)])))),
         ## subset(donnee$Robs,source==1 & annee%in%(annees+donnee$lagBH))[,'valeur'] / objReport$qRecru[1],
         ## subset(donnee$Robs,source==2 & annee%in%(annees+donnee$lagBH))[,'valeur'] / objReport$qRecru[2],
         ## subset(donnee$Robs,source==3 & annee%in%(annees+donnee$lagBH))[,'valeur'] / objReport$qRecru[3])),
         pch=21, bg=2, xlab=labB85, ylab=labRecru, axes=FALSE); graphics::box()
    graphics::axis(1, at=pretty(graphics::graphics::par('usr')[1:2]), labels=pretty(graphics::graphics::par('usr')[1:2])/1000000)
    graphics::axis(2, at=pretty(graphics::graphics::par('usr')[3:4]), labels=pretty(graphics::graphics::par('usr')[3:4])/1000)
    graphics::rug(utils::tail(objReport$Bpred,donnee$lagBH)); graphics::rug(utils::head(objReport$Rpred,donnee$lagBH), side=2)
    ## for(i in (donnee$lagBH+1):length(objReport$Rpred)){
    ##     graphics::points(objReport$Bpred[i-donnee$lagBH], objReport$Rbh[i-donnee$lagBH])
    ##     graphics::points(objReport$Bpred[i-donnee$lagBH], objReport$Rpred[i], col=2)
    ## }
    graphics::abline(h=0, col='grey70'); graphics::abline(v=0, col='grey70')
    ## if(donnee$switchBH){
    ##     graphics::curve(objReport$alphaBH*(x)/ (objReport$betaBH + x), add=TRUE, col=1, lty=2, lwd=2)
    ##     graphics::legend('bottomright', inset=0.03, legend=paste(c('alpha','beta'), '=', round(c(objReport$alphaBH, objReport$betaBH))))
    ## }else{
    graphics::curve(alphaBH*(x)/ (betaBH + x), add=TRUE, col=1, lty=2, lwd=2)
    ##
    ##
    eq <- substitute('R'==a%*%'B85+'/(b+'B85+'),list(a = round(alphaBH), b = round(betaBH)))
    graphics::text(0.5*graphics::graphics::par('usr')[2],0.9*graphics::graphics::par('usr')[4],labels=eq)
    ##
    graphics::points(utils::head(objReport$Bpred[annees+1],-(donnee$lagBH)), utils::tail(objReport$Rpred[annees],-(donnee$lagBH)), type='o', pch=21, bg=2)
    graphics::abline(a=0,b=alphaBH/betaBH,col=2, lty=3)
    graphics::lines(c(0,msyVal$B0), c(0,alphaBH*(msyVal$B0)/ (betaBH + msyVal$B0)), col=2, lty=3)
    CR <- (alphaBH/betaBH)/((alphaBH*(msyVal$B0)/ (betaBH + msyVal$B0))/msyVal$B0)
    graphics::legend('bottomright', inset=0.03, legend=paste(c('alpha','beta','CR','h'), '=', c(round(c(alphaBH, betaBH)),round(c(CR,CR/(4+CR)),2))))
    ##
    noms <- utils::head(donnee$anneesFittees[annees],-(donnee$lagBH))
    noms[which(!noms %in% pretty(noms))] <- NA
    graphics::text(utils::head(objReport$Bpred[annees+1],-(donnee$lagBH)), utils::tail(objReport$Rpred[annees],-(donnee$lagBH)),
         labels=noms, pos=3)
    graphics::abline(v=msyVal$B, lty=2)
    graphics::abline(v=as.numeric(msyVal$B)*c(0.4,0.8), lty=3, col='grey70')
    graphics::text(x=mean(graphics::graphics::par('usr')[c(1,2)]), y=diff(graphics::graphics::par('usr')[c(3,4)])*0.9 + graphics::graphics::par('usr')[3], label='A', cex=1.5)
  }
  ##
  if(3%in%lesquels){
    ## SSR r\U{00E9}sidus
    annees <- donnee$anneesBH[1]:donnee$anneesBH[2]
    Rtemp <- utils::tail(objReport$Rpred[annees],-(donnee$lagBH))
    Btemp <- utils::head(objReport$Bpred[annees+1],-(donnee$lagBH))
    ## if(donnee$switchBH){
    ##     alpha <- objReport$alphaBH
    ##     beta <- objReport$betaBH
    ## }else{
    if(is.null(msyVal)){
      msyVal <- calculerBH(donnee, objReport, steep=FALSE)
    }
    if(!steep){
      alphaBH <- msyVal$alphaBH
      betaBH <- msyVal$betaBH
      sigmaBH <- msyVal$sigmaBH
    }else{
      alphaBH <- bh$alphaBH.steep
      betaBH <- bh$betaBH.steep
      sigmaBH <- bh$sigmaBH.steep
    }
    ## }
    plot(utils::head(donnee$anneesFittees[annees],-(donnee$lagBH)), (log(Rtemp) - log((alphaBH*Btemp / (betaBH + Btemp))))/sigmaBH,
         xlim=range(donnee$anneesFittees), ylim=ylimLog, xlab=labAn, ylab=labRes)
    graphics::abline(h=0, col='grey')
    for(i in seq_along(Rtemp)){
      graphics::lines(rep(donnee$anneesFittees[annees[i]], 2),
            c(0,log(Rtemp[i]) - log((alphaBH*Btemp[i] / (betaBH + Btemp[i]))))/sigmaBH, lwd=2, col=2)
      graphics::points(donnee$anneesFittees[annees[i]], (log(Rtemp[i]) - log((alphaBH*Btemp[i] / (betaBH + Btemp[i]))))/sigmaBH,
             pch=21, bg=2)
    }
    graphics::text(x=mean(graphics::graphics::par('usr')[c(1,2)]), y=diff(graphics::graphics::par('usr')[c(3,4)])*0.9 + graphics::graphics::par('usr')[3], label='B', cex=1.5)
    ##
  }
  ##
  if(4%in%lesquels){
    ## beverton-holt et R
    Rest <- msyVal$alphaBH*objReport$Bpred/(msyVal$betaBH+objReport$Bpred)
    plot(c(min(donnee$anneesFittees)-1,donnee$anneesFittees), Rest,
         type='o', pch=1, xlim=range(donnee$anneesFittees)-c(donnee$lagBH,0),
         ylim=c(0,max(Rest)),
         xlab=labAn, ylab=labRecru); graphics::abline(h=0, col='grey70')
    graphics::polygon(donnee$anneesFittees[donnee$anneesBH[c(1,1,2,2)]]-c(0,0,donnee$lagBH-1,donnee$lagBH-1)+c(-0.5,-0.5), graphics::graphics::par('usr')[c(3,4,4,3)], col='grey90', border=NA); graphics::box()
    graphics::lines(c(min(donnee$anneesFittees)-1,donnee$anneesFittees), Rest, type='o')
    graphics::lines(donnee$anneesFittees-donnee$lagBH, objReport$Rpred, type='o', pch=1, col=2)
    graphics::legend('topleft', inset=0.03, legend=c('BH','Recrues'), col=c(1,2), pch=1, lty=1)
    graphics::abline(v=donnee$anneesBH)
    ##
    ## plot(c(min(donnee$anneesFittees)-1,donnee$anneesFittees), log(Rest),
    ##      type='o', pch=1, xlim=range(donnee$anneesFittees)-c(donnee$lagBH,0),
    ##      ylim=c(min(log(objReport$Rpred)),max(log(Rest))),
    ##      xlab=labAn, ylab='Log-recrues'); graphics::abline(h=0, col='grey70')
    ## graphics::polygon(donnee$anneesFittees[donnee$anneesBH[c(1,1,2,2)]]-c(0,0,donnee$lagBH-1,donnee$lagBH-1)+c(-0.5,-0.5), graphics::graphics::par('usr')[c(3,4,4,3)], col='grey90', border=NA); graphics::box()
    ## graphics::lines(c(min(donnee$anneesFittees)-1,donnee$anneesFittees), log(Rest), type='o')
    ## graphics::lines(donnee$anneesFittees-donnee$lagBH, log(objReport$Rpred), type='o', pch=1, col=2)
    ## graphics::legend('topleft', inset=0.03, legend=c('BH','Recrues'), col=c(1,2), pch=1, lty=1)
    ## graphics::abline(v=donnee$anneesBH)
    ##
  }
}
##
graph.ageLong <- function(donnee, param, objReport=NULL, msyVal=NULL, langue=c('fr','en','bil')){
  switch(langue,
         'fr' = {labAge <- 'Age'; labTaille <- 'Taille (cm)'; labPoids <- 'Poids (kg)'; labPoidsA <- 'Poids (kg), \U{00E2}ge a'; labPoidsB <- 'Poids (kg), \U{00E2}ge a+1'},
         'en' = {labAge <- 'Age'; labTaille <- 'Size (cm)'; labPoids <- 'Weight (kg)'; labPoidsA <- 'Weight (kg), age a'; labPoidsB <- 'Weight (kg), \U{00E2}ge a+1'},
         'bil' = {labAge <- 'Age'; labTaille <- 'Taille/Size (cm)'; labPoids <- 'Poids/Weight (kg)'; labPoidsA <- 'Poids/Weight, age a'; labPoidsB <- 'Poids/Weight, age a+1'})
  ##
  ## age-longueur
  plot(donnee$croiss$age, donnee$croiss$longueur, xlim=c(0,max(donnee$croiss$age)), ylim=c(0,max(donnee$croiss$longueur)), xlab=labAge, ylab=labTaille)
  graphics::abline(v=0, col='grey70'); graphics::abline(h=0, col='grey70'); graphics::abline(h=c(81,85), col=3); graphics::abline(v=donnee$lagBH, col=3)
  graphics::axis(4, at=c(81,85), labels=c(81,85), col.ticks=3, col.axis=3, las=1, cex.axis=0.6)
  graphics::curve(objReport$linf * (1-exp(-objReport$K * (x-(objReport$t0)))), add=TRUE, col=2)
  graphics::legend('bottomright', inset=0.03, legend=paste(c('Linf=','K=','t0='),
                                                 c(round(objReport$linf), round(objReport$K,3), round(objReport$t0,2))))
  if(FALSE){ #graph des r\U{00E9}sidus
    plot(donnee$croiss$age, (donnee$croiss$longueur)-(objReport$linf * (1-exp(-objReport$K * (donnee$croiss$age-(objReport$t0))))),
         main='Residus croissance', xlab='Age', ylab='R\U{00E9}sidus'); graphics::abline(h=0)
    ## plot(donnee$croiss$age, log(donnee$croiss$longueur)-log(objReport$linf * (1-exp(-objReport$K * (donnee$croiss$age-(objReport$t0))))),
    ##      main='Residus croissance', xlab='Age', ylab='R\U{00E9}sidus'); graphics::abline(h=0)
  }
  graphics::text(x=mean(graphics::par('usr')[c(1,2)]), y=diff(graphics::par('usr')[c(3,4)])*0.9 + graphics::par('usr')[3], label='A', cex=1.5)
  ##
  ## age-poids
  plot(donnee$croiss$age, donnee$croiss$poids, xlim=c(0,max(donnee$croiss$age)), xlab=labAge, ylab=labPoids)
  graphics::abline(h=0, col='grey70'); graphics::abline(h=donnee$lpAlpha*c(81,85)^donnee$lpBeta, col=3); graphics::abline(v=donnee$lagBH, col=3)
  graphics::axis(4, at=donnee$lpAlpha*c(81,85)^donnee$lpBeta, labels=c(81,85), col.ticks=3, col.axis=3, las=1, cex.axis=0.6)
  graphics::curve(donnee$lpAlpha*(objReport$linf * (1-exp(-objReport$K * (x-(objReport$t0)))))^donnee$lpBeta, add=TRUE, col=2)
  ## graphics::legend('topleft', inset=0.03, legend=paste(c('Winf=','WK=','t0='),
  ##                                                c(round(donnee$lpAlpha*objReport$linf^donnee$lpBeta), round(objReport$K,3), round(objReport$t0,2))))
  graphics::text(x=mean(graphics::par('usr')[c(1,2)]), y=diff(graphics::par('usr')[c(3,4)])*0.9 + graphics::par('usr')[3], label='B', cex=1.5)
  ##
  ## poids+t vs poids_t+1
  taille.age <- objReport$linf * (1-exp(-objReport$K * (1:100-(objReport$t0))))
  poids.age <- donnee$lpAlpha*taille.age^donnee$lpBeta
  poids.min <- min(donnee$omega$valeur); poids.max <- max(donnee$omega$valeur)
  plot(utils::head(poids.age,-1), utils::tail(poids.age,-1), xlim=c(0,poids.age[utils::head(which(poids.age>poids.max),1)+2]), ylim=c(0,poids.age[utils::head(which(poids.age>poids.max),1)+3]),
       xlab=labPoidsA, ylab=labPoidsB)
  ## plot(utils::head(poids.age,-1), utils::tail(poids.age,-1), xlim=c(0,200), ylim=c(0,200),
  ##      xlab='Poids moyen, \U{00E2}ge a', ylab='Poids moyen, \U{00E2}ge a+1')
  graphics::abline(v=donnee$lpAlpha*85^donnee$lpBeta)
  graphics::axis(3, at=donnee$lpAlpha*85^donnee$lpBeta, label='85 cm')
  graphics::curve(objReport$alpha+x*objReport$rho, from=poids.min, to=poids.max, add=TRUE, col=4, lwd=3)
  graphics::abline(a=objReport$alpha, b=objReport$rho)
  graphics::points(utils::head(poids.age,-1), utils::tail(poids.age,-1))
  ## graphics::abline(a=0,b=1, col='grey70', lty=2)
  graphics::legend('bottomright', inset=0.03, legend=paste(c('alpha=','rho='), c(round(objReport$alpha,2), round(objReport$rho,3))))
  noms <- 1:20; noms[which(!noms %in% pretty(noms))] <- NA
  graphics::text(utils::head(poids.age,-1)[noms], utils::tail(poids.age,-1)[noms], labels=noms, pos=1)
  graphics::text(x=mean(graphics::par('usr')[c(1,2)]), y=diff(graphics::par('usr')[c(3,4)])*0.9 + graphics::par('usr')[3], label='C', cex=1.5)
}
##
graph.B.retro <- function(retro, langue=c('fr','en','bil')){
  switch(langue,
         'fr' = {labAn <- 'Ann\U{00E9}e'; labBiom <- "Biomasse ('000t)"; labBmsy <- 'Biomasse / Brms'},
         'en' = {labAn <- 'Year'; labBiom <- "Biomass ('000t)"; labBmsy <- 'Biomass / Bmsy'},
         'bil' = {labAn <- 'Ann\U{00E9}e/Year'; labBiom <- "Biomasse/Biomass ('000t)"; labBmsy <- 'Biomasse/Brms - Biomass/Bmsy'})
  ##
  ## biomasse absolue
  temp <- retro[[1]]
  plot(c(min(temp$donnee$anneesFittees)-1, temp$donnee$anneesFittees),
       temp$objReport$Bpred/1000/1000,
       type='l', xlim=range(temp$donnee$anneesFittees),
       ylim=c(0,max(temp$objReport$Bpred))/1000/1000,
       xlab=labAn, ylab=labBiom, col=grDevices::grDevices::rainbow(length(retro))[1]); graphics::abline(h=0, col='grey70')
  graphics::points(utils::tail(temp$donnee$anneesFittees,1), utils::tail(temp$objReport$Bpred/1000/1000,1), col=grDevices::grDevices::rainbow(length(retro))[1])
  qrel <- temp$objReport$qRel[1]
  for(i in 2:length(retro)){
    temp <- retro[[i]]
    graphics::lines(c(min(temp$donnee$anneesFittees)-1, temp$donnee$anneesFittees),
          temp$objReport$Bpred/1000/1000, type='l', col=grDevices::rainbow(length(retro))[i])
    graphics::points(utils::tail(temp$donnee$anneesFittees,1), utils::tail(temp$objReport$Bpred/1000/1000,1), col=grDevices::rainbow(length(retro))[i])
    qrel[i] <- temp$objReport$qRel[1]
  }
  ##
  graphics::legend('topleft', inset=0.03, legend=paste('qRel =',round(qrel,2)), lty=1, col=grDevices::rainbow(length(retro)))
  ##
  ##
  ## biomasse sur Bmsy
  temp <- retro[[1]]
  plot(c(min(temp$donnee$anneesFittees)-1, temp$donnee$anneesFittees),
       temp$objReport$Bpred/temp$msy$B,
       type='l', xlim=range(temp$donnee$anneesFittees),
       ylim=c(0,2),
       xlab=labAn, ylab=labBmsy, col=grDevices::rainbow(length(retro))[1]); graphics::abline(h=1); graphics::abline(h=c(0.4,0.8), lty=2)
  BsurBmsy <- utils::tail(temp$objReport$Bpred/temp$msy$B,1)
  graphics::points(utils::tail(temp$donnee$anneesFittees,1), BsurBmsy, col=grDevices::rainbow(length(retro))[1])
  for(i in 2:length(retro)){
    temp <- retro[[i]]
    graphics::lines(c(min(temp$donnee$anneesFittees)-1, temp$donnee$anneesFittees),
          temp$objReport$Bpred/temp$msy$B, type='l', col=grDevices::rainbow(length(retro))[i])
    BsurBmsy[i] <- utils::tail(temp$objReport$Bpred/temp$msy$B,1)
    graphics::points(utils::tail(temp$donnee$anneesFittees,1), BsurBmsy[i], col=grDevices::rainbow(length(retro))[i])
  }
  ##
  graphics::legend('topleft', inset=0.03, legend=paste('Bfin / Bmsy =',round(BsurBmsy,2)), lty=1, col=grDevices::rainbow(length(retro)))
}
