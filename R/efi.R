#' Ensemble Fractional Imputation

#' @importFrom dplyr filter "%>%" group_by n select summarise
#' @importFrom plyr alply adply
#' @import CVXR
#' @import cat
#' @import igraph
#' @import cvam
#'
#' @param Y categorical data with missingnes.Missing valuesa are encoded as NA. If Freq = T, the last column serves as the integer frequencies of the grouped observations.
#' @param dp a class doublep computed from \code{doublep}
#' @param freq an optional logical value that indicates whether frequency column exists in the data frame Y.
#' @return A class  \code{"efi"}
#' @examples
#' data(HairEyeColor)
#' p = 3
#' Y = do.call("rbind", apply(as.data.frame.table(HairEyeColor), 1, function(x) matrix(rep(x[1:p], each = x[p+1]), nc = p)))
#' Y = as.data.frame(Y)
#' for(k in 1:p){
#' Y[[k]] = factor(Y[[k]])
#' }
#' names(Y) <- names(dimnames(HairEyeColor))
#' # (n = nrow(Y)); sum(HairEyeColor)
#' delta = matrix(rbinom(n * p, 1, 0.5), nr = n, nc = p)
#' Y[delta == 0] = NA
#'
#' cand.edges = as.list(data.frame(combn(p, 2)))
#' dp = doublep(Y, cand.edges, freq = F)
#' EFI = efi(Y, dp, freq = F)
#' @export

efi = function(Y, dp, freq = F){
  if(class(dp) != "doublep") stop("dp must have class doublep")
  weightvec = dp$weightvec
  weightmat = dp$weightmat
  cand.edges = dp$cand.edges
  supp = dp$supp
  n = dp$n
  p = dp$p
  marginalProb = dp$marginalProb
  namesY = dp$namesY

  if(!freq) Y = cbind(Y, Freq = 1)
  else colnames(Y)[ncol(Y)] <- "Freq"

  namestmp = names(cbind(Y, id = 1:nrow(Y)))
  Y_FI0 = adply(cbind(Y, id = 1:nrow(Y)), 1, function(Z){
    if(any(is.na(Z))){
      tmpY = cbind(expand.grid(lapply(Z[,is.na(Z),drop = F], levels)), Z[,!is.na(Z)])[namestmp]
      return(tmpY)
    }else{
      Z
    }
  })

  marginalProbmat = sapply(Y[-ncol(Y)], function(x) cvamLik(~ x, cvam(~ x, data = Y, freq = Freq), data = Y)$likVal)
  marginalProbmat = as.data.frame(marginalProbmat)
  tmpmat = sapply(cand.edges[weightvec != 0], function(x){
    apply(select(marginalProbmat, -x), 1, prod) *
      cvamLik(formula(paste("~", paste(namesY[x], collapse = "+"))), cvam(formula(paste("~", paste(namesY[x], collapse = "*"))), data = Y, freq = Freq), data = Y)$likVal
  })
  tmpmat = as.data.frame(tmpmat)

  marginalProbmat2 =sapply(names(Y)[-ncol(Y)], function(x){
    form = as.formula(paste("~", x))
    cvamLik(form, cvam(form, data = Y, freq = Freq), data = Y_FI0)$likVal
  })
  marginalProbmat2 = as.data.frame(marginalProbmat2)

  tmpmat2 = sapply(cand.edges[weightvec != 0], function(x){
    apply(select(marginalProbmat2, -x), 1, prod) *
      cvamLik(formula(paste("~", paste(namesY[x], collapse = "+"))), cvam(formula(paste("~", paste(namesY[x], collapse = "*"))), data = Y, freq = Freq), data = Y_FI0)$likVal
  })
  tmpmat2 = as.data.frame(tmpmat2)

  timestmp = unlist(table(Y_FI0$id))
  restmp = mapply(function(x, y) return(x / rep(y, times = timestmp)), tmpmat2, tmpmat)

  w = apply(restmp, 1, function(x) sum(x * weightvec))
  Y_FI0 = cbind(Y_FI0, w = w * Y_FI0$Freq)

  edges = cand.edges[weightvec != 0]
  EFI = list(imp = Y_FI0, edges = edges, weightmat = weightmat, n = n)

  class(EFI) <- "efi"
  return(EFI)
}
