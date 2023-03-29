# library(foreach)
# library(doParallel)
# library(CVXR)
# library(cat)
# # library(mice)
# library(igraph)
# library(plyr)
# library(dplyr)
# library(cvam)

doublep = function(Y, cand.edges, freq = F){
  if(!freq) Y = cbind(Y, Freq = 1)

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
    apply(select(marginalProbmat, -x), 1, prod) *
      cvamLik(formula(paste("~", paste(namesY[x], collapse = "+"))), cvam(formula(paste("~", paste(namesY[x], collapse = "*"))), data = Y, freq = Freq), data = Y)$likVal
  })

  w = Variable(length(cand.edges))
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

estimate = function(data, expr){
  call <- match.call()
  imp = data$imp
  n = data$n
  formula = formula(paste(expr, "~ 1"))
  # model.summ = suppressWarnings(summary(glm(as.integer((Var1 == 1) & (Var2 == 2)) ~ 1, family = binomial(link = "identity"), data = imp, weights = w)))
  model.summ = suppressWarnings(summary(glm(formula, family = binomial(link = "identity"), data = imp, weights = w)))
  # print(model.summ$terms)
  Estimate = model.summ$coefficients[1]
  U = model.summ$coefficients[2]^2
  B = eval(parse(text = paste("sum(imp %>% group_by(id) %>% dplyr::summarise(sum(w[", expr, "]), sum(w)) %>% {ifelse(.[[3]] == 0, 0, .[[2]] * (.[[3]] - .[[2]]) / .[[3]])}) / n^2")), envir = imp)

  Std.Error = sqrt(U + B)
  return(list(call = call, expr = expr, Estimate = Estimate, Std.Error = Std.Error, sqrtU = sqrt(U), sqrtB = sqrt(B)))
}

plot.doublep = function(dp){
  weightmat = dp$weightmat
  g <- graph_from_adjacency_matrix(weightmat, mode = "lower", weighted = "weight")
  # plot(g, edge.width = E(g)$weight, edge.label = E(g)$weight, label = coords)
  plot(g, edge.width = E(g)$weight, edge.label = E(g)$weight)
}

