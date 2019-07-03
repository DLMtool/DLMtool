library(DLMtool)
DLMextra()
library(DLMextra)

OMs <- avail("OM")

OM <- get(OMs[1])


Stock <- SubOM(OM, "Stock")
Fleet <- SubOM(OM, "Fleet")
Obs <- SubOM(OM, "Obs")
Imp <- Overages

OM@nsim <- 5
plot(OM)



plot(Imp, dev=TRUE)


plot(Fleet, Stock=Stock)


plot("Effort", Fleet, Stock=Stock, dev=TRUE)

plot("Catchability", Fleet, Stock=Stock, dev=TRUE)


plot("Selectivity", Fleet, Stock=Stock, dev=TRUE)

plot("MPA", Fleet, Stock=Stock, dev=TRUE)



Hist <- runMSE(OM, Hist=TRUE)
plot(Hist, dev=TRUE)

slotNames(Hist)

Hist@TSdata %>% names()

Hist@AtAge %>% names()

plot(OM)


# add plot.Hist function













