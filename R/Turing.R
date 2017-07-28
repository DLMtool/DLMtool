
Turing <- function(Data, OM) {
  if (class(Data) != "Data") stop("Data must be class 'Data'")
  if (class(OM) != "OM") stop("OM must be class 'OM'")
  
  if (length(Data@Year) != OM@nyears) {
    message("Note: length Data@Year (", length(Data@Year), ") is not of length OM@nyears (", OM@nyears, ") \nUsing last ", length(Data@Year), " years of simulations")
  } # fix this for when Data is longer than OM 
  
  
  nyr <- length(Data@Year)
  nyears <- OM@nyears <- length(Data@Year)
  # Run historical simulations
  Hist <- runMSE(OM, Hist=TRUE)
  
  sims <- sample(1:OM@nsim, 5)
  
  # What Data is available
  Cat <- Hist$TSdata$Catch[(nyears-nyr+1):nyears,sims]
  meancat <- matrix(apply(Cat, 2, mean), nrow=nyr, ncol=length(sims), byrow=TRUE)
  Cat <- Cat/meancat
  Cat_d <- as.numeric(Data@Cat/mean(Data@Cat, na.rm=TRUE))
  
  Cat <- cbind(Cat, Cat_d)
  
  
  ind <- sample(1:ncol(Cat), ncol(Cat))
  Cat <- Cat[,ind]
  
  par(mfrow=c(3,2))
  for (X in 1:ncol(Cat)) plot(Cat[,X], type="l", ylim=c(0, max(Cat)))
 
  
  
  dim(Data@CAA)
  
  dim(Data@CAL)
  
  
  
  
  
}