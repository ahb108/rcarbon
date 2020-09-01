# normalise densities to 1
normalise_densities <- function(dens,eps) {
  dens <- dens/sum(dens)
  dens[dens < eps] <- 0
  dens <- dens/sum(dens)
  return(dens)
}

# calibrates in F14C space
F14C_calibration <- function(age, error, calf14,calf14error, eps) {
  # F14 <- exp(calcurve[,2]/-8033) 
  # F14Error <-  F14*calcurve[,3]/8033 
  # calf14 <- approx(calcurve[,1], F14, xout=calBP)$y 
  # calf14error <-  approx(calcurve[,1], F14Error, xout=calBP)$y 
  f14age <- exp(age/-8033) 
  f14err <- f14age*error/8033 
  p1 <- (f14age - calf14)^2 
  p2 <- 2 * (f14err^2 + calf14error^2) 
  p3 <- sqrt(f14err^2 + calf14error^2) 
  dens <- exp(-p1/p2)/p3 
  dens[dens < eps] <- 0
  return(dens)
}

# calibrates in 14C BP space
BP14C_calibration <- function(age, error, mu, tau2, eps) {
  tau <- error^2 + tau2
  dens <- dnorm(age, mean=mu, sd=sqrt(tau))
  dens[dens < eps] <- 0
  return(dens)
}

# reads a cal curve file from extdata
read_cal_curve_from_file <- function(calCurves) {
  calCurveFile <- paste(system.file("extdata", package="rcarbon"), "/", calCurves,".14c", sep="")
  options(warn=-1)
  calcurve <- readLines(calCurveFile, encoding="UTF-8")
  calcurve <- calcurve[!grepl("[#]",calcurve)]
  calcurve.con <- textConnection(calcurve)
  calcurve <- as.matrix(read.csv(calcurve.con, header=FALSE, stringsAsFactors=FALSE))[,1:3]
  close(calcurve.con)
  options(warn=0)
  colnames(calcurve) <- c("CALBP","C14BP","Error")
  return(calcurve)
}
