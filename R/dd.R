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
  d$perteTag$perteAnnuelle2tag <- c(temp,utils::tail(temp,1))
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
  par$transQrelGSL <- stats::qlogis  (0.7/2) #biomasse minimale chalutable des RV-mpo
  if(length(unique(donnee$Bobs$source)) > 1) par$transQrelAutre <- rep(stats::qlogis  (0.000005/2),length(unique(donnee$Bobs$source))-1) #autres indices utilis\U{00E9}s
  par$transRapportQrecru <- rep(stats::qlogis (1/2), length(unique(donnee$Robs$source))) #utiliser le rapport avec le q utilis\U{00E9} pour la biomasse minimale chalutable
  par$transTauxExp <- rep(stats::qlogis ((0.1-0.001)/0.9), length(annees)) #captures divis\U{00E9} par biomasse
  par$logBpred <- rep(log(2e7), length(annees)+1)
  par$logRpred <- rep(log(1e6), length(annees))
  ## par$logB0 <- par$logBpred[1]
  par$logN0 <- log(exp(par$logBpred[1])/stats::quantile(donnee$omega$valeur, probs=0.25))
  ## par$transOmega0 <- 0
  ##
  ## poids-longueur
  par$logLinf <- log(223)
  par$logK <- log(0.053)
  par$transT0 <- 0
  ##
  return(par)
}

