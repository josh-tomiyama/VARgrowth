### Create several distribution here
### For now, assume no covariate in observed process.

CreateVARGrowthDist <- function(){

}


NormalDist <- function(mean_struct = "arp", var_struct = "iid"){
  list(dist = "Normal",
       est_function = GompertzEstFunction_nlme,
       generator = rnorm()
       )
}

TruncNormalDist <- function(bounds = c("upper", "lower", "both")){
  list(dist = "Normal",
       est_function = BoundedGompertzEstFunction_nlme,
       generator = rtruncnorm()
  )
}

ARMean <- function(arp = 1, y){
  n <- length(y) - arp
}

ARMeanNormalDist <- function(ARp = 1){

}

ARVarNormalDist <- function(ARVarp = 1){
  list(dist = "Normal")
}

ARCHNormalDist <- function(ARCHp = 1){

}
