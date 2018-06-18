

### Fix:
# - nyears != nyears Data
# - missing data in Catch or Index - linearly interpolate?
# - missing years of CAA and CAL
# - max age for CAA 
# - number of size bins 
# Turing(OM, Data)

# CAA and CAL should match Year - work on Data object
# recommended by Andre Punt 

Turing <- function(OM, Data, wait=TRUE) {
  if (class(OM) != "OM") stop("OM must be class 'OM'")
  if (class(Data) != "Data") stop("Data must be class 'Data'")

  # Generate Historical Data 
  message("Simulating Data")
  # if length data exists, make sure length bins are the same
  if (!all(is.na(Data@CAL[1,,]))) { 
    CAL_bins <- Data@CAL_bins
    OM@cpars$CAL_bins <- CAL_bins
  }
  # make sure maxage is the same
  OM@cpars$maxage <- Data@MaxAge
  
  nsamp <- 5
  OM@nsim <- nsamp +2 
  SimDat <- runMSE(OM, Hist=TRUE, silent = TRUE)$Data
  # if(!all(SimDat@CAL_bins == Data@CAL_bins)) stop("CAL_bins not correct length")
  
  YrInd <- 1:OM@nyears
  
  if (length(Data@Year) < OM@nyears) {
    message("Note: length Data@Year (", length(Data@Year), ") is < length OM@nyears (", OM@nyears, ") \nUsing last ",
            length(Data@Year), " years of simulations")
    YrInd <- (OM@nyears - length(Data@Year)+1): OM@nyears

  } else if (length(Data@Year) > OM@nyears) {
    stop("OM@nyears is < length(Data@Year). Should OM@nyears be increased to match Data?", call. = FALSE) 
  }
  set.seed(runif(1, 1, 10001))
  samps <- sample(1:OM@nsim, nsamp)
  message("Randomly sampling 5 iterations")
  
  # ---- Catch Data ---- 
  plotTSdata(Ylab="St. Catch", slot="Cat", message="Catch Data", Data, SimDat, 
             samps, YrInd, wait)
  
  # ---- Index Data ---- 
  plotTSdata(Ylab="St. Index", slot="Ind", message="Index Data", Data, SimDat, 
             samps, YrInd, wait)
  
  # ---- Recruitment Data ----
  plotTSdata(Ylab="St. Recruitment", slot="Rec", message="Recruitment Data", Data, 
             SimDat, samps, YrInd, wait)
  
  # ---- Mean Length Data ----
  plotTSdata(Ylab="Mean Length", slot="ML", message="Mean Length Data", Data, 
             SimDat, samps, YrInd, wait, standarise=FALSE)
  
  # ---- Mean Length Data ----
  plotTSdata(Ylab="Lbar", slot="Lbar", message="Lbar Data", Data, 
             SimDat, samps, YrInd, wait, standarise=FALSE)
  
  # ---- Catch-at-Age ----
  plotCAAdata(Ylab="Count", slot="CAA", message="Catch-at-Age Data", Data, 
              SimDat, samps, YrInd, wait)
  
  # ---- Catch-at-Length ----
  plotCAAdata(Ylab="Count", slot="CAL", message="Catch-at-Length Data", Data, 
              SimDat, samps, YrInd, wait)
}

    
plotCAAdata <- function(Ylab="Count", slot="CAA", message="Catch-at-Age Data", 
                        Data, SimDat, samps, YrInd, wait) {
  realDat <- slot(Data, slot)[1,,]
  sampDat <- slot(SimDat, slot)[samps,YrInd,]
  
  if (all(is.na(realDat))) {
    message("No ", message, " found in Data object")
    return()
  }
  # plot last 4 years
  nyr <- nrow(realDat)
  realDat <- realDat[(nyr-3):nyr,]
  sampDat <- sampDat[,(nyr-3):nyr,]
  
  zeros <- dim(sampDat)[3]-dim(realDat)[2]
  if (zeros>0) {
    zeromat <- matrix(0, nrow=nrow(realDat), ncol=zeros)
    message('The number of columns in `CAA` is less than `Data@MaxAge`. Filling with 0s')
    realDat <- cbind(realDat, zeromat)
  }
 
  if (!all(is.na(realDat))) {
    message("Plotting: ", message)
   
    if (slot == "CAA") {
      nyrs <- nrow(realDat); maxage <- ncol(realDat)
      dimnames(realDat) <- list(1:nyrs, 1:maxage)
      
      df1 <- as.data.frame.table(realDat, stringsAsFactors = FALSE)
      colnames(df1) <- c("Year", "Val", "Freq")
      
      dimnames(sampDat) <- list(1:nsamp, 1:nyrs, 1:maxage)
      df2 <- as.data.frame.table(sampDat, stringsAsFactors = FALSE)
      colnames(df2) <- c("Sim", "Year", "Val", "Freq")
      
      Xlab <- "Age"
    } else if (slot == "CAL") {
      nyrs <- nrow(realDat); nbins <- length(Data@CAL_bins) - 1
      By <- Data@CAL_bins[2] - Data@CAL_bins[1]
      BinsMid <- seq(Data@CAL_bins[1] + 0.5*By, by=By,length.out = nbins)
      dimnames(realDat) <- list(1:nyrs, BinsMid)
      
      df1 <- as.data.frame.table(realDat, stringsAsFactors = FALSE)
      colnames(df1) <- c("Year", "Val", "Freq")
      
      dimnames(sampDat) <- list(1:nsamp, 1:nyrs, BinsMid)
      df2 <- as.data.frame.table(sampDat, stringsAsFactors = FALSE)
      colnames(df2) <- c("Sim", "Year", "Val", "Freq")
      
      Xlab <- "Length"
    }

    
    df1$Sim <- as.character(dim(sampDat)[1]+1)
    df3 <- dplyr::bind_rows(df1, df2)
    
    df3$Sim <- as.factor(df3$Sim)
    levels(df3$Sim) <- sample(1:6, 6) # randomise
    df3$Sim <- factor(df3$Sim, levels=1:(nsamp+1), ordered = TRUE)
    
    df3$Val <- as.numeric(df3$Val)
    df3$Year <- as.numeric(df3$Year)
    df3$Freq <- as.numeric(df3$Freq)

    realInd <- unique(df3$Sim[unique(match(df3$Freq, df1$Freq))])
    realInd <- as.numeric(as.character(realInd[!is.na(realInd)] ))

    Years <- unique(df3$Year)
    SampYears <- (max(Years)-3):max(Years)
    
    df4 <- df3 %>% filter(Year %in% SampYears)
    
    # filter empty bins
    Bcount <- df4 %>% dplyr::group_by(Val) %>% dplyr::summarize(count=sum(Freq>0))
    indmin <- min(which(Bcount$count >0))
    indmax <- max(which(Bcount$count >0))
    Bmin <- Bcount$Val[max(indmin-1, 1)]
    Bmax <- Bcount$Val[min(indmax+1, length(Bcount$Val))]
    df4 <- df4 %>% dplyr::filter(Val >= Bmin & Val <= Bmax)
    
    
    P1 <- ggplot2::ggplot(df4, ggplot2::aes(x=Val, y=Freq)) + ggplot2::geom_bar(stat="identity") +
      ggplot2::facet_grid(Sim~Year) + ggplot2::theme_classic() +
      ggplot2::labs(x=Xlab, y=Ylab) + ggplot2::ggtitle(message)
    
    df5 <- df4 %>% dplyr::filter(Sim == realInd)
    P2 <- P1 + ggplot2::geom_bar(fill="blue", data=df5, stat="identity")
    print(P1)
    if (wait) tt <- readline(paste0("Press any key to show real ", message, "..."))
    print(P2)
    if (wait)  tt <- readline("Press any key to continue...")
  
  } else {
    message("No ", message, " found in Data object")
  }
  
}


plotTSdata <- function(Ylab, slot, message, Data, SimDat, samps, YrInd, 
                       wait=TRUE, standarise=TRUE) {
  realDat <- slot(Data, slot)[1,]
  # remove NAs
  YrInd2 <- YrInd[is.finite(realDat)] 
  YrInd3 <- Data@Year[is.finite(realDat)]
  realDat <- realDat[is.finite(realDat)]
  sampDat <- slot(SimDat, slot)[samps,YrInd2]
  
  if (!all(is.na(realDat))) {
    message("Plotting: ", message)
    CombDat <- rbind(sampDat, realDat)
    CombDat <- CombDat[sample(1:nrow(CombDat)),]
    realInd <- which(match(data.frame(t(CombDat)), data.frame(realDat)) == 1)
    if (standarise) CombDat <- CombDat/matrix(apply(CombDat, 1, mean, na.rm=TRUE), 
                                              nrow=nrow(CombDat), ncol=ncol(CombDat))
    rownames(CombDat) <- rep("", nrow(CombDat))
    tCombDat <- as.data.frame(t(CombDat))
    df <- tidyr::gather(tCombDat, factor_key=TRUE)
    df$Year <- factor(rep(YrInd3, nrow(CombDat)))
    
    P1 <- ggplot2::ggplot(df, ggplot2::aes(x=Year, y=value, group=1)) + ggplot2::geom_line(size=1.2) +
      ggplot2::facet_wrap(~key, nrow=2) + ggplot2::theme_classic() +
      ggplot2::labs(x="Year", y=Ylab) + ggplot2::ggtitle(message) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) 
    
    
    df2 <- df %>% dplyr::filter(key == paste0("V", realInd))
    P2 <- P1 + ggplot2::geom_line(size=1.2, data=df2, color='blue')
    print(P1)
    if (wait) tt <- readline(paste0("Press any key to show real ", message, "..."))
    print(P2)
    if (wait) tt <- readline("Press any key to continue...")
    
  } else {
    message("No ", message, " found in Data object")
  }
}
  

  
  
