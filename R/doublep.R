#' Double project data from candidate models
#'
#' @param Y categorical data with missingnes.Missing valuesa are encoded as NA. If Freq = T, the last column serves as the integer frequencies of the grouped observations.
#' @param cand.edges candidate edges that are considered in the candidate graphical exponential model. Typicially, as.list(data.frame(combn(ncol(Y), 2))) is used.
#' @param freq an optional logical value that indicates whether frequency column exists in the data frame Y.
#' @return A class  \code{"doubledp"}
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
#' @export

doublep = function(Y, cand.edges, freq = F){
  if(!freq) Y = cbind(Y, Freq = 1)
  else colnames(Y)[ncol(Y)] <- "Freq"

  supp = lapply(Y[-ncol(Y)], . %>% levels)
  n = sum(Y$Freq)
  p = length(supp)
  supplen = sapply(supp, length)
  namesY = names(Y[-ncol(Y)])

  # grid_tmp = expand.grid(supp)
  # marginalProb =lapply(Y[-ncol(Y)], function(x) xtabs(Freq ~ x, Y, addNA = F) / sum(xtabs(Freq ~ x, Y, addNA = F)))
  # mat = sapply(cand.edges, function(x){
  #   thetahat <- array(get.fitted(cvam(formula(paste("~", paste(namesY[x], collapse = "*"))), data = Y, freq = Freq))$fit, dim = supplen[x])
  #   dimnames(thetahat) <- supp[x]
  #
  #   list_tmp = append(list(thetahat), marginalProb[-x])
  #
  #   apply(grid_tmp, 1, function(y)
  #     prod(c(list_tmp[[1]][t(y[x])], mapply(function(z, y) z[y], list_tmp[2:(p-1)], y[-(x)]))))
  # })
  # tmpmat = apply(Y[-ncol(Y)], 1, function(y){
  #   colSums(mat[apply(grid_tmp, 1, function(x) all(x[which(!is.na(y))] == y[which(!is.na(y))] )), ,drop = F])
  # })
  # tmpmat = t(tmpmat)

  marginalProb = lapply(Y[-ncol(Y)], function(x) get.fitted( cvam(~ x, data = Y, freq = Freq), mfTrue  = F))

  marginalProbmat =sapply(Y[-ncol(Y)], function(x) cvamLik(~ x, cvam(~ x, data = Y, freq = Freq), data = Y)$likVal)
  marginalProbmat = as.data.frame(marginalProbmat)

  tmpmat = sapply(cand.edges, function(x){
    # print(x)
    apply(select(marginalProbmat, -x), 1, prod) *
      cvamLik(formula(paste("~", paste(namesY[x], collapse = "+"))), cvam(formula(paste("~", paste(namesY[x], collapse = "*"))), data = Y, freq = Freq), data = Y)$likVal
  })

  w = CVXR::Variable(length(cand.edges))
  constraints <- list(sum(w) == 1, w >= 0)
  Phi_R <- Maximize(sum(Y$Freq * log(tmpmat %*% w)))
  prob <- Problem(Phi_R, constraints)
  res <- solve(prob, solver = "ECOS")

  if(res$status != "optimal" & res$status != "optimal_inaccurate") stop(paste("CVXR solver error: res$status =", res$status))

  weightvec = round(c(res$getValue(w)), 3)
  weightvec = weightvec / sum(weightvec)

  weightmat = matrix(0, nr = p, nc = p)
  weightmat[t(sapply(cand.edges, cbind))] <- weightvec # This should be corrected
  weightmat[lower.tri(weightmat)] <- t(weightmat)[lower.tri(weightmat)]
  dp <- list(weightvec = weightvec, weightmat = weightmat, #mat = mat, grid_tmp = grid_tmp,
             cand.edges = cand.edges, supp = supp, n = n, p = p, marginalProb = marginalProb,
             namesY = namesY)
  class(dp) <- "doublep"
  return(dp)
}
