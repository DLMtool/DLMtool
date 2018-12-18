# Solve the Baranov catch equation 
# http://api.admb-project.org/baranov_8cpp.html
# Catch_age_area=relyield * RelRec
# Vuln_age=VAge
# B_age_area = lx * WtAge * VAge * RelRec
# Catch_age_area = Catch_tot[sim,,]
# M_at_Age = M_ageArray[sim,,y]
# Vuln_age=V_P[sim,,y]
# B_age_area=CurrentB[sim,,]
# maxF=maxF
# byage=TRUE
CalculateF <- function(Catch_age_area, M_at_Age, Vuln_age, B_age_area, maxF=NULL, byage=FALSE) {
  if (is.null(dim(B_age_area))) {
    nareas <- 1
    B_age_area <- matrix(B_age_area, ncol=nareas, nrow=length(B_age_area))
    Catch_age_area <- matrix(Catch_age_area, ncol=nareas, nrow=length(Catch_age_area))
  } else {
    nareas <- dim(B_age_area)[2]
  }
  
  # initial guess based on Pope's approximation
  Far <- (Catch_age_area)/((B_age_area) * exp(-0.5*M_at_Age)) # F by age & area 
  Fi <- sum(Catch_age_area)/(sum(Vuln_age*(B_age_area * exp(-0.5*M_at_Age)))) # F overall
  Fi[Fi > 1] <- 0.99

  M_at_Age_area <- matrix(M_at_Age, ncol=nareas, nrow=length(M_at_Age),byrow=FALSE)
  
  for (x in 1:50) {
    Z_age_area <- M_at_Age_area + Far
    Z <- M_at_Age + Fi*Vuln_age
    
    # predicted catch 
    pC_age_area <- Far/Z_age_area * B_age_area * (1-exp(-Z_age_area))
    pC <- sum((Fi*Vuln_age)/Z * rowSums(B_age_area) * (1-exp(-Z)))
    
    # calculate derivative 
    deriv <- ((B_age_area * (1-exp(-(Z_age_area))))/Z_age_area +
                (B_age_area * Far*Vuln_age * exp(-(Z_age_area)))/Z_age_area -
                (B_age_area * Far*Vuln_age * (1-exp(-(Z_age_area))))/Z_age_area^2)
    deriv[!is.finite(deriv)] <- 0
    Far <- Far - ((pC_age_area - (Catch_age_area))/deriv)
    Fi <- Fi - ((pC - sum(Catch_age_area))/sum(deriv))
    diff <- max(abs(pC - sum(Catch_age_area)))
    if (diff < 0.01) break()
  }
  if (byage) {
    if(!is.null(maxF))
      Far[Far>maxF] <- maxF
    return(Far)
  }
  if (!byage) {
    if(!is.null(maxF)) Fi <- min(Fi, maxF)
    return(Fi)
  }
}
