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
  if(!freq) Y = cbind(Y, 1)

  colnames(Y)[ncol(Y)] = "Freq"
  supp = lapply(Y[-ncol(Y)], . %>% levels)
  grid_tmp = expand.grid(supp)
  n = sum(Y$Freq)
  p = length(supp)
  supplen = sapply(supp, length)

  marginalProb =lapply(Y[-ncol(Y)], function(x) xtabs(Freq ~ x, Y, addNA = F) / sum(xtabs(Freq ~ x, Y, addNA = F)))
  namesY = names(Y[-ncol(Y)])

  mat = sapply(cand.edges, function(x){
    thetahat <- array(get.fitted(cvam(formula(paste("~", paste(namesY[x], collapse = "*"))), data = Y, freq = Freq))$fit, dim = supplen[x])
    dimnames(thetahat) <- supp[x]

    list_tmp = append(list(thetahat), marginalProb[-x])

    apply(grid_tmp, 1, function(y)
      prod(c(list_tmp[[1]][t(y[x])], mapply(function(z, y) z[y], list_tmp[2:(p-1)], y[-(x)]))))
  })

  tmpmat = apply(Y[-ncol(Y)], 1, function(y){
    colSums(mat[apply(grid_tmp, 1, function(x) all(x[which(!is.na(y))] == y[which(!is.na(y))] )), ,drop = F])
  })
  tmpmat = t(tmpmat)

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
  dp <- list(weightvec = weightvec, weightmat = weightmat, mat = mat,
             cand.edges = cand.edges, supp = supp, grid_tmp = grid_tmp, n = n, p = p)
  class(dp) <- "doublep"
  return(dp)
}

efi = function(Y, dp, freq = F){
  if(class(dp) != "doublep") stop("dp must have class doublep")
  weightvec = dp$weightvec
  weightmat = dp$weightmat
  mat = dp$mat
  cand.edges = dp$cand.edges
  supp = dp$supp
  grid_tmp = dp$grid_tmp
  n = dp$n
  p = dp$p

  if(!freq) Y = cbind(Y, 1)


  Y_FI0 = adply(cbind(Y, 1:nrow(Y)), 1, function(Z){
    y = Z[-c(p+1, p+2)]
    Freq = Z[p+1]
    id = Z[p+2]
    mis_idx = which(is.na(y))
    if(length(mis_idx) == 0){
      return(cbind(y, Freq, id = id))
    }else{
      obs_idx0 = (1:p)[-mis_idx]
      obs_idx = which(colSums(weightmat[mis_idx,, drop = F]) != 0)
      obs_idx = obs_idx[!(obs_idx %in% mis_idx)]

      grid_mis = expand.grid(supp[mis_idx])
      mat_FI = adply(grid_mis, 1, function(z){
        # mis_val = unlist(z)
        # obs_val0 = unlist(y[obs_idx0])
        # obs_val = unlist(y[obs_idx])
        mis_val = z
        obs_val0 = y[obs_idx0]
        obs_val = y[obs_idx]
        y_tmp = y
        y_tmp[mis_idx] <- mis_val; y_tmp[obs_idx0] <- obs_val0;
        # print(y_tmp)

        sel_cand = sapply(cand.edges, function(x) all(x %in% c(mis_idx, obs_idx)))

        # print(grid_tmp[,c(mis_idx, obs_idx),drop = F])
        # print(z)
        # print(cbind(mis_val, obs_val))
        # print(apply(grid_tmp[,c(mis_idx, obs_idx),drop = F], 1, function(x) all(x == cbind(mis_val, obs_val))) )
        # condP = colSums(mat[apply(t(grid_tmp[,c(mis_idx, obs_idx),drop = F]) == c(mis_val, obs_val), 2, all), sel_cand ,drop = F]) / colSums(mat[apply(t(grid_tmp[,obs_idx,drop = F]) == obs_val, 2, all), sel_cand, drop = F]) # P(Y2 = 1, Y4 = 2 | Y1 = 1, Y2 = 2, Y3 = 2)
        condP = colSums(mat[apply(grid_tmp[,c(mis_idx, obs_idx),drop = F], 1, function(x) all(x == cbind(mis_val, obs_val))), sel_cand ,drop = F]) / colSums(mat[apply(grid_tmp[,obs_idx,drop = F], 1, function(x) all(x == obs_val)), sel_cand, drop = F]) # P(Y2 = 1, Y4 = 2 | Y1 = 1, Y2 = 2, Y3 = 2)

        # print(paste("P(", paste(paste("Y", mis_idx, sep = "") , mis_val, sep = "=", collapse = ", "), "|", paste(paste("Y", obs_idx, sep = "") , obs_val, sep = "=", collapse = ", "), ")"));
        sumtmp = sum(weightvec[sel_cand])
        if(sumtmp == 0){
          weightvec_tmp = rep(1, length(weightvec[sel_cand])) / length(weightvec[sel_cand])
        }else{
          weightvec_tmp = weightvec[sel_cand] / sumtmp
        }

        if(length(condP) != 0){
          fweight = sum(condP * weightvec_tmp)
        }else{
          # fweight = marginalProb[[mis_idx]][mis_val]
          fweight = sapply(marginalProb[mis_idx], function(prob_tmp) prob_tmp[mis_val])
        }

        # print(fweight)
        # print(cbind(y_tmp, fweight * Freq, id))
        cbind(y_tmp, Freq = fweight * Freq, id = id)
      })
      return(mat_FI)
    }
  })
  # print(Y_FI0)

  # Y_FI0 = do.call("rbind", Y_FI0)
  # Y_FI0 = as.data.frame(Y_FI0, stringsAsFactors = FALSE)
  # Y_FI0[,p+1] <- as.numeric(Y_FI0[,p+1])
  names(Y_FI0)[c(p+1, p+2)] <- c("w", "id")
  # names(Y_FI0)[1:p] <- names(Y)

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

