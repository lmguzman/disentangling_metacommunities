hmsc_fun <- function(run_one, landscape){
  
  if(nrow(run_one)  < 2) {
    all_res <- c(space = NA, time = NA, environment = NA)
  }else{
  
    run_one$Species <- paste("s",run_one$Species, sep = "")
    run_one$Patch <- paste("p",run_one$Patch, sep = "")
    
    run_one_wide <- run_one %>% 
      dplyr::select(N, Species, Patch, Time, env) %>% 
      arrange(Time, Patch, Species) %>% 
      spread(key = Species, value = N, fill = 0)
    }
  
  if(nrow(run_one_wide) <2){
    all_res <- c(space = NA, time = NA, environment = NA)
  }else{
    
    studyDesign = data.frame(Space = factor(run_one_wide$Patch), 
                             Time = factor(paste("t", run_one_wide$Time, sep = "")), 
                             Space_Time = factor(paste(run_one_wide$Patch, "_t",run_one_wide$Time, sep = "")))
    
    rownames(landscape) <- paste("p",1:100, sep = "")
    landscape2 <- landscape
    landscape2$Patch <- rownames(landscape)
    xyt <- data.frame(left_join(run_one_wide, landscape2) %>% 
                        select(x, y, Time))
    
    rownames(xyt) <- studyDesign$Space_Time
    space_time_rL <- HmscRandomLevel(sData = xyt, sMethod = "NNGP", nNeighbours = 10)
    space_time_rL <- setPriors(space_time_rL,nfMin=1,nfMax=1)
    
    space_rL <- HmscRandomLevel(sData = landscape, sMethod = "NNGP", nNeighbours = 10)
    space_rL <- setPriors(space_rL,nfMin=1,nfMax=1)
    
    time.df <- data.frame(Time = unique(run_one_wide$Time))
    rownames(time.df) <- paste("t", unique(run_one_wide$Time), sep = "")
    time_rL <- HmscRandomLevel(sData = time.df)
    time_rL <- setPriors(time_rL,nfMin=1,nfMax=1)     
    
    # Add species composition in the previous time (Y1)?
    # Y1<- dplyr::select(run_one_wide, -Time, -Patch, -env, -x, -y)[time!=max(time),]
    # names(Y1) <- paste("X.", names(Y1),sep="")
    # time0 <- unique(time)
    # time2 <- c()
    # for(i in 1:nrow(Y1)){
    #   wt0 <- which(time0==time[time!=max(time)][i])
    #   time2[i] <- time0[wt0+1]
    # }
    # run_one_wide$time.patch <- paste(run_one_wide$Time,run_one_wide$Patch,sep="")
    # Y1$time.patch <- paste(time2,patch[time!=max(time)],sep="")
    # run_one_wide2 <- left_join(run_one_wide[time!=min(time),], Y1, by="time.patch")
    # del.na <- unique(which(is.na(run_one_wide2), arr.ind=TRUE)[,1])
    # if(length(del.na) != 0){
    #   run_one_wide2 <- run_one_wide2[-del.na,]
    # }
    
    Y<- data.matrix(dplyr::select(run_one_wide, -Time, -Patch, -env))
    
    XData <- data.frame(env = run_one_wide$env)
    
    mod <- Hmsc(Y = Y, 
                XFormula = ~ poly(env, degree = 2, raw = TRUE), 
                XData = XData, 
                studyDesign = studyDesign,
                ranLevels = list(Space_Time = space_time_rL, Space = space_rL, Time = time_rL), 
                distr = "poisson")
    
    thin = 1
    samples = 1000
    nChains = 4
    set.seed(1)
    m <- sampleMcmc(mod, thin = thin, samples = samples, transient = 5000, 
                    nChains = nChains, nParallel = nChains, updater=list(GammaEta=FALSE))
    
    mpost = convertToCodaObject(m)     

    predY = computePredictedValues(m, expected=FALSE)
    MF = evaluateModelFit(hM=m, predY=predY)
    VP = computeVariancePartitioning(m,group=c(1,1), groupnames = c("env"))
    
    OmegaCor = computeAssociations(m)
  
    
    all_res <- list(MF, VP, OmegaCor, mpost$Alpha, m$nc, m$ns)  
  }
  
  return(all_res)
}
