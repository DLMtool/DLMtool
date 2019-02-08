library(DLMtool)
library(kableExtra)
OM <- testOM
OM@nsim <- 10
MSE <- runMSE(OM, MPs=NA)



PLimit1 <- P50
formals(PLimit1)$Yrs <- -10
class(PLimit1) <- "PM"

PLimit2 <- P50
formals(PLimit2)$Yrs <- c(11,50)
class(PLimit2) <- "PM"


PMlist <- list('PLimit1', 'PLimit2')
Pthresh <- 0.8

maketable <- function(MSE, ..., PMlist=NULL) {
  
  if (is.null(PMlist)) {
    PMlist <- unlist(list(...))
  } else {
    PMlist <- unlist(PMlist)
  }
  for (X in seq_along(PMlist))
    if (!PMlist[X] %in% avail("PM")) stop(PMlist[X], " is not a valid PM method")
  runPM <- vector("list", length(PMlist))
  for (X in 1:length(PMlist)) {
    runPM[[X]] <- eval(call(PMlist[X], MSE))
  }
  nMPs <- MSE@nMPs
  
  df <- data.frame(MP=lapply(runPM, function(x) x@MPs) %>% unlist(),
             Prob=lapply(runPM, function(x) x@Mean) %>% unlist(),
             Cap=rep(lapply(runPM, function(x) x@Caption) %>% unlist(), each=nMPs))
  df$Prob <- round(df$Prob,2)
  
  temp <- df %>% group_by(MP) %>% summarize(min=min(Prob))
  df <- left_join(df, temp, by='MP') %>% arrange(MP)

  df <- df %>% mutate(
    Prob = cell_spec(Prob, "html", color = ifelse(Prob < Pthresh, "red", "green")),
    MP = cell_spec(MP, "html", color = ifelse(min < Pthresh, "red", "green")),
  )
  df <- tidyr::spread(df, Cap, Prob) 
  ind <- which(df$min >= Pthresh)
  df$min <- NULL
  df %>% kable(format = "html", escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    row_spec(ind, bold = T)
    
}
