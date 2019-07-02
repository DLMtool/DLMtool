library(DLMtool)
DLMextra()
library(DLMextra)

OMs <- avail("OM")

OM <- get(OMs[1])


Stock <- SubOM(OM, "Stock")
Fleet <- SubOM(OM, "Fleet")
Obs <- SubOM(OM, "Obs")
Imp <- SubOM(OM, "Imp")

plot(Obs, dev=TRUE)


plot(Fleet, Stock=Stock)


plot("Effort", Fleet, Stock=Stock, dev=TRUE)

plot("Catchability", Fleet, Stock=Stock, dev=TRUE)


plot("Selectivity", Fleet, Stock=Stock, dev=TRUE)

plot("MPA", Fleet, Stock=Stock, dev=TRUE)



# add plot.Hist function













