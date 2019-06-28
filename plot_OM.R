library(DLMtool)
DLMextra()
library(DLMextra)

OMs <- avail("OM")

OM <- get(OMs[1])


Stock <- SubOM(OM, "Stock")
Fleet <- SubOM(OM, "Fleet")
Obs <- SubOM(OM, "Obs")
Imp <- SubOM(OM, "Imp")



plot("Effort", Fleet, Stock=Stock, dev=TRUE)

plot("Catchability", Fleet, Stock=Stock, dev=TRUE)

















