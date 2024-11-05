testfct <- function(x){print(x)}

fnll <- function(param, fit=TRUE){
  getAll(param, donnee)
  ##
  ## Ecart-type des fonctions vraissemblances
  sigma_Bobs <- exp(logSigma_Bobs)
  sigma_Bproc <- exp(logSigma_Bproc)
  sigma_oBar <- exp(logSigma_oBar)
  sigma_C <- exp(logSigma_C)
  sigma_Rrw <- exp(logSigma_Rrw)
  sigma_Robs <- exp(logSigma_Robs)
  sigma_longAge <- exp(logSigma_longAge)
  sigma_retourTag <- exp(logSigma_retourTag)
  ##
  Bobs <- OBS(Bobs)
  ##
  ## Parametres variables dans le temps
  Bpred <- exp(logBpred) #Bpred inclus B0, donc donc commence a l'an 0
  Rpred <- exp(logRpred)
  tauxExp <- 0.001 + 0.9*plogis(transTauxExp) # entre 0.001 et 0.9
  F <- -log(1-tauxExp)
  Z <- F+M
  s <- exp(-M) * (1-tauxExp) #equivalent a exp(-Z)
  ##
  ## Parametres invariables dans le temps
  N0 <- exp(logN0)
  qRel <- c(2*plogis(transQrelGSL), 2*plogis(transQrelAutre))
  qRecru <- 2*plogis(transRapportQrecru) * qRel[1] #par rapport au q des adultes
  ## qRecru <- rep(qRel[1],length(unique(Robs$source)))
  linf <- exp(logLinf)
  K <- exp(logK)
  t0 <- -5 + 10*plogis(transT0) #entre -5 et 5
  ##
  ## calcul des parametre du graphique Ford-Walford (rho et alpha), selon l'etendue des longueurs moyennes observes
  poidsMoyMin <- min(omega$valeur)
  ageMin <- log(1-((poidsMoyMin/lpAlpha)^(1/lpBeta))/linf) / -K + t0
  poidsMoyMin.plus1an <- lpAlpha * (linf * (1-exp(-K * ((ageMin+1)-t0))))^lpBeta
  ##
  poidsMoyMax <- max(omega$valeur)
  ageMax <- log(1-((poidsMoyMax/lpAlpha)^(1/lpBeta))/linf) / -K + t0
  poidsMoyMax.plus1an <- lpAlpha * (linf * (1-exp(-K * ((ageMax+1)-t0))))^lpBeta
  ##
  rho <- (poidsMoyMax.plus1an - poidsMoyMin.plus1an) / (poidsMoyMax - poidsMoyMin)
  alpha <- poidsMoyMin.plus1an - poidsMoyMin*rho
  ## winf <- lpAlpha * (linf)^lpBeta
  ## rho <- (poidsMoyMin.plus1an - winf) / (poidsMoyMin - winf)
  ## alpha <- winf*(1 - rho)
  ##
  ## Initialisation a l'an 1
  Bpred.proc <- s[1] * alpha * N0 +
    s[1] * rho * Bpred[1] + #Bpred commence a l'an 0
    omegaK[1,'valeur'] * Rpred[1]
  Npred <- s[1] * N0 + Rpred[1]
  Cpred <- tauxExp[1] * Bpred[1] * exp(-M) #peche apres mortalite naturelle
  ##
  ## Progression annuelle
  for(i in 2:(length(Bpred)-1)){
    Bpred.proc[i] <- s[i] * alpha * Npred[i-1] +
      s[i] * rho * Bpred[i] +
      omegaK[i,'valeur'] * Rpred[i]
    Npred[i] <- s[i] * Npred[i-1] + Rpred[i]
    Cpred[i] <- tauxExp[i] * Bpred[i] * exp(-M) #peche apres mortalite naturelle
    if(i==a2010){ # ajuster la biomasse pour le changement de taille legale en 2010
      Bpred.proc[i] <- Bpred.proc[i] * drop2010
      ## soustraction du nombre de poisson en moins, selon la biomasse soustraite et le poids moyen a cette taille
      Npred[i] <- Npred[i] - Bpred.proc[i] * (1-drop2010) / poidsMoy81a85
    }
  }
  omegaPred <- tail(Bpred,-1)/Npred
  ##
  ## suivi des tags presents dans l'eau et recaptures
  nTag <- matrix(nrow=length(Bpred)-1, ncol=length(Bpred)-1)
  nTagRetourPred <- matrix(NA, nrow=length(Bpred)-1, ncol=length(Bpred)-1) #pas de captures l'annee de marquage
  for(i in 1:(nrow(nTag)-1)){
    nTag[i,i] <- nTagsPoses$valeur[i] * sPostMarquage
    for(j in (i+1):ncol(nTag)){ #remplir le triangle superieur
      nTag.temp <- nTag[i,j-1] * (1-perteTag[j-1,'perteAnnuelle2tag'])
      nTag[i,j] <- nTag.temp * s[j]
      nTagRetourPred[i,j] <- nTag.temp * exp(-M) * tauxExp[j] * tauxRetour
    }
  }
  nTag[nrow(nTag),ncol(nTag)] <- nTagsPoses$valeur[nrow(nTag)] * sPostMarquage
  ##
  ##
  ## erreur de processus sur la biomasse, retirer le premier Bpred
  nll.Bproc <- -sum(dnorm(tail(logBpred,-1), log(Bpred.proc), sigma_Bproc, log=TRUE), na.rm=TRUE)
  ##
  ## marche aleatoire walk du recrutement
  nll.recru <- -sum(dnorm(log(Rpred[-length(Rpred)]), log(Rpred[-1]), sigma_Rrw, log=TRUE), na.rm=TRUE)
  ##
  ## erreur d'ajustement de la longueur a l'age
  nll.longAge <- -sum(dnorm(croiss$longueur, linf*(1-exp(-K*(croiss$age-t0))), sigma_longAge, log=TRUE), na.rm=TRUE)
  ##
  ## erreur d'observation sur les captures
  nll.Cobs <- -sum(dnorm(log(Cobs[,'valeur']), log(Cpred), sigma_C, log=TRUE), na.rm=TRUE)
  ##
  ## erreur d'observation sur les indices d'abondance
  nll.Bobs <- 0
  for(i in 1:nrow(Bobs)){
    if(is.finite(log(Bobs[i,2]))){
      annee <- Bobs[i,1]
      source <- Bobs[i,3]
      sigma <- Bobs[i,4]
      nll.Bobs <- nll.Bobs - dnorm(log(Bobs[i,2]), log(Bpred[annee+1]*qRel[source]), sigma_Bobs[sigma], log=TRUE)
    }
  }
  ##
  ## erreur d'observation sur les indices de recrutement
  nll.Robs <- 0
  for(i in 1:nrow(Robs)){
    if(is.finite(log(Robs[i,2]))){
      annee <- Robs[i,1]
      source <- Robs[i,3]
      sigma <- Robs[i,4]
      nll.Robs <- nll.Robs - dnorm(log(Robs[i,2]), log(Rpred[annee]*qRecru[source]), sigma_Robs[sigma], log=TRUE)
    }
  }
  ##
  ## erreur d'observation sur les poids moyens
  nll.oBar <- 0
  for(i in 1:nrow(omega)){
    if(is.finite(log(omega[i,2]))){
      annee <- omega[i,1]
      source <- omega[i,3]
      sigma <- omega[i,4]
      ## nll.oBar <- nll.oBar - dnorm(log(omega[i,2]), log(omegaPred[annee]), sigma_oBar[sigma], log=TRUE)
      nll.oBar <- nll.oBar - dnorm(omega[i,2], omegaPred[annee], sigma_oBar[sigma], log=TRUE)
    }
  }
  ##
  ## erreur d'observation sur les retours d'etiquettes
  nll.tag <- 0
  for(i in 1:nrow(nTagsRetourObs)){
    anPose <- nTagsRetourObs[i,'anneePose']
    anRecap <- nTagsRetourObs[i,'anneeRecap']
    source <- nTagsRetourObs[i,'source']
    nll.tag <- nll.tag - dnorm(nTagRetourPred[anPose,anRecap], nTagsRetourObs[i,'valeur'], sigma_retourTag[source], log=TRUE)
  }
  ##
  ## vraissemblance totale
  nll <- nll.Bproc +
    nll.recru +
    nll.longAge +
    nll.Cobs +
    nll.Bobs +
    nll.Robs +
    nll.oBar +
    nll.tag
  ##
  ##
  ## sorties
  REPORT(nll.Bproc)
  REPORT(nll.recru)
  REPORT(nll.Robs)
  REPORT(nll.Bobs)
  REPORT(nll.oBar)
  REPORT(nll.Cobs)
  REPORT(nll.longAge)
  REPORT(nll.tag)
  REPORT(N0)
  REPORT(Bpred)
  REPORT(Bpred.proc)
  REPORT(Npred)
  REPORT(Rpred)
  REPORT(s)
  REPORT(tauxExp)
  REPORT(qRel)
  REPORT(qRecru)
  REPORT(Cpred)
  REPORT(omegaPred)
  REPORT(rho)
  REPORT(alpha)
  REPORT(M)
  REPORT(F)
  REPORT(Z)
  REPORT(linf)
  REPORT(K)
  REPORT(t0)
  REPORT(nTag)
  REPORT(nTagRetourPred)
  REPORT(tauxRetour)
  REPORT(sPostMarquage)
  ##
  REPORT(nll)
  ##
  return(nll)
}

calculerDonnee <- function(donneeInit, annees, valM=0.15, valTR=0.8, Rmin=3000, anOmega.min=1998, doublerC=FALSE, RobsSigma.unique=FALSE,
                           omegaSigma.unique=FALSE) {
  d <- donneeInit
  d$anneesFittees <- donneeInit$anneesFittees[which(donneeInit$anneesFittees%in%annees)] #ann\U{00E9}es utilis\U{00E9}es
  d$anneesFitteesID <- seq_along(d$anneesFittees); names(d$anneesFitteesID) <- d$anneesFittees
  d$anneesBH <- range(d$anneesFitteesID[which(d$anneesFittees%in%1990:max(d$anneesFittees))]) #ann\U{00E9}es min et max sur lesquelles ajuster la SRR
  d$Cobs <- subset(donneeInit$Cobs, annee%in%d$anneesFittees); d$Cobs$annee <- d$anneesFitteesID[match(d$Cobs$annee,d$anneesFittees)] #obs des captures
  if(doublerC) d$Cobs <- rbind(d$Cobs,d$Cobs)
  d$Bobs <- subset(donneeInit$Bobs, annee%in%d$anneesFittees); d$Bobs$annee <- d$anneesFitteesID[match(d$Bobs$annee,d$anneesFittees)]  #obs des abondances
  d$Bobs <- subset(d$Bobs, source%in%1:3)
  d$Robs <- subset(donneeInit$Robs, annee%in%d$anneesFittees); d$Robs$annee <- d$anneesFitteesID[match(d$Robs$annee,d$anneesFittees)]  #obs des recrues
  d$Robs$sigma <- d$Robs$source; d$Robs[which(d$Robs$valeur<Rmin),'valeur'] <- Rmin
  d$omega <- subset(donneeInit$omega, annee%in%d$anneesFittees); d$omega$annee <- d$anneesFitteesID[match(d$omega$annee,d$anneesFittees)]  #obs poids moy
  d$omega <- subset(d$omega, annee%in%d$anneesFitteesID[which(d$anneesFittees>=anOmega.min)]); d$omega$sigma <- d$omega$source;
  d$omegaK <- subset(donneeInit$omegaK, annee%in%d$anneesFittees)
  d$omegaK$annee <- d$anneesFitteesID[match(d$omegaK$annee,d$anneesFittees)] #poids moyen au recrutement (variable selon la taille minimale l\U{00E9}gale)
  d$nTagsPoses <- subset(donneeInit$nTagsPoses, annee%in%d$anneesFittees)
  d$nTagsPoses$annee <- d$anneesFitteesID[match(d$nTagsPoses$annee,d$anneesFittees)]  #nombre de fl\U{00E9}tans \U{00E9}tiquett\U{00E9}s par ann\U{00E9}e
  d$nTagsRetourObs <- subset(donneeInit$nTagsRetourObs, anneePose%in%d$anneesFittees & anneeRecap%in%d$anneesFittees);
  d$nTagsRetourObs$anneePose <- d$anneesFitteesID[match(d$nTagsRetourObs$anneePose,d$anneesFittees)] #nb de fl\U{00E9}tan recap, an marquage et an recap
  d$nTagsRetourObs$anneeRecap <- d$anneesFitteesID[match(d$nTagsRetourObs$anneeRecap,d$anneesFittees)] #nb de fl\U{00E9}tan recap par an marquage et an recap
  ## d$nTagsRetour <- donneeInit$nTagsRetour[as.character(d$anneesFittees)] #nombre total de retour d'\U{00E9}tiquettes par ann\U{00E9}e
  ## ajuster les num\U{00E9}ro de source si certaines sont manquantes, utile surtout si certains relev\U{00E9} ne sont pas consid\U{00E9}r\U{00E9}s
  for(i.source in c('Bobs','Robs','omega')){
    if(length(unique(d[[i.source]]$source))<max(d[[i.source]]$source)){
      sourceInit <- sort(unique(d[[i.source]]$source))
      for(i in seq_along(sourceInit)){
        d[[i.source]][d[[i.source]]$source==sourceInit[i], 'source'] <- i
        d[[i.source]][d[[i.source]]$sigma==sourceInit[i], 'sigma'] <- i
      }
    }
  }
  if(RobsSigma.unique) d$Robs$sigma <- 1
  if(omegaSigma.unique) d$omega$sigma <- 1
  ##
  ## donn\U{00E9}es non influenc\U{00E9} par le choix d'ann\U{00E9}es:
  if(FALSE){
    d$a2010 #indice de l'ann\U{00E9}e du changement de taille minimale l\U{00E9}gale
    d$drop2010 #taux de diminution de la biomasse li\U{00E9}e au changement de taille l\U{00E9}gale
    d$poidsMoy81a85 #poids moyen des fl\U{00E9}tan entre 81 et 85 cm
    d$lagBH         #nombre d'ann\U{00E9}es entre la biomasse 85+ et le recrutement produit par celle-ci
    d$lpAlpha #parametres de la relation poids-longueur de forme p=a*x^b
    d$lpBeta
  }
  ## d$croiss <- subset(d$croiss, sexe%in%c('f','femelle'))    #longueur-age-poids
  d$M <- valM
  d$sPostMarquage <- 0.95
  d$tauxRetour <- valTR
  temp <- diff(d$perteTag$perteCummul2tag)
  d$perteTag$perteAnnuelle2tag <- c(temp,tail(temp,1))
  ##
  return(d)
}
calculerParam <- function(donnee, logSigma_C=NULL){
  ## calcul des parametres d'entr\U{00E9}e
  annees <- donnee$anneesFittees
  par <- list()
  ##
  par$logSigma_Bobs <- rep(0, length(unique(donnee$Bobs$sigma))) #la premi\U{00E8}re source est la biomasse minimale chalutable des RV-mpo
  par$logSigma_Bproc <- 0
  par$logSigma_oBar <- rep(0, length(unique(donnee$omega$sigma)))
  if(is.null(logSigma_C)){par$logSigma_C <- 0}else{par$logSigma_C <- logSigma_C}
  par$logSigma_Rrw <- 0
  par$logSigma_Robs <- rep(0, length(unique(donnee$Robs$sigma)))
  par$logSigma_longAge <- 0
  par$logSigma_retourTag <- 0
  ##
  par$transQrelGSL <- qlogis(0.7/2) #biomasse minimale chalutable des RV-mpo
  if(length(unique(donnee$Bobs$source)) > 1) par$transQrelAutre <- rep(qlogis(0.000005/2),length(unique(donnee$Bobs$source))-1) #autres indices utilis\U{00E9}s
  par$transRapportQrecru <- rep(qlogis(1/2), length(unique(donnee$Robs$source))) #utiliser le rapport avec le q utilis\U{00E9} pour la biomasse minimale chalutable
  par$transTauxExp <- rep(qlogis((0.1-0.001)/0.9), length(annees)) #captures divis\U{00E9} par biomasse
  par$logBpred <- rep(log(2e7), length(annees)+1)
  par$logRpred <- rep(log(1e6), length(annees))
  ## par$logB0 <- par$logBpred[1]
  par$logN0 <- log(exp(par$logBpred[1])/quantile(donnee$omega$valeur, probs=0.25))
  ## par$transOmega0 <- 0
  ##
  ## poids-longueur
  par$logLinf <- log(223)
  par$logK <- log(0.053)
  par$transT0 <- 0
  ##
  return(par)
}

