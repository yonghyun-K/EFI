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
#' (n = nrow(Y)); sum(HairEyeColor)
#' delta = matrix(rbinom(n * p, 1, 0.5), nr = n, nc = p)
#' Y[delta == 0] = NA
#'
#' cand.edges = as.list(data.frame(combn(p, 2)))
#' dp = doublep(Y, cand.edges, freq = F)
#' @export

doublep = function(Y, edges_list, freq = F, R = 5){
  if(!freq) Y = cbind(Y, Freq = 1)
  else colnames(Y)[ncol(Y)] <- "Freq"

  supp = lapply(Y[-ncol(Y)], . %>% levels)
  n = sum(Y$Freq)
  p = length(supp)
  supplen = sapply(supp, length)
  namesY = names(Y[-ncol(Y)])

  # sampleIdx = sample(1:n, round(n / 2))
  # Y1 = Y[sampleIdx, ]
  # Y2 = Y[-sampleIdx, ]

  # marginalProb = lapply(Y[-ncol(Y)], function(x) get.fitted( cvam(~ x, data = Y, freq = Freq), mfTrue  = F))

  marginalProbmat = sapply(colnames(Y)[-ncol(Y)], function(x) {
    formula = as.formula(paste("~", x))
    cvamLik(formula, cvam(formula, data = Y, freq = Freq), data = Y)$likVal
  })
  marginalProbmat = as.data.frame(marginalProbmat)

  edges_list1 = NULL
  edges_list2 = edges_list

  weightveclist = NULL
  for(r in 1:R){
    print(r)

    tmpmat = sapply(c(edges_list1, edges_list2), function(x){
      # print(formula(paste("~", paste(sapply(x, function(z) paste(namesY[z], collapse = "*")), collapse = "+"))))
      apply(select(marginalProbmat, -unique(unlist(x))), 1, prod) *
        cvamLik(formula(paste("~", paste(namesY[unique(unlist(x))], collapse = "+"))),
                cvam(formula(paste("~", paste(sapply(x, function(z) paste(namesY[z], collapse = "*")), collapse = "+"))),
                     data = Y, freq = Freq), data = Y)$likVal
    })

    w = CVXR::Variable(ncol(tmpmat))
    constraints <- list(sum(w) == 1, w >= 0)
    Phi_R <- Maximize(sum(Y$Freq * log(tmpmat %*% w)))
    prob <- Problem(Phi_R, constraints)
    res <- solve(prob, solver = "ECOS")
    # res <- solve(prob)
    # print(paste("res$solver", res$solver))
    # CVXR::installed_solvers()

    if(res$status != "optimal" & res$status != "optimal_inaccurate"){
      res <- solve(prob, solver = "ECOS_BB")
    }
    if(res$status != "optimal" & res$status != "optimal_inaccurate"){
      res <- solve(prob, solver = "SCS")
    }
    if(res$status != "optimal" & res$status != "optimal_inaccurate"){
      res <- solve(prob, solver = "OSQP")
    }
if(res$status != "optimal" & res$status != "optimal_inaccurate"){
  stop(paste("CVXR solver error: res$status =", res$status))
}


    # print(paste("res$value", res$value))

    weightvec = round(c(res$getValue(w)), 3)
    weightvec = weightvec / sum(weightvec)
    # print(paste("sum(Y$Freq * log(tmpmat %*% w))", sum(Y$Freq * log(tmpmat %*% weightvec))))
    tmpidx = setdiff(1:length(weightvec), sequence(length(edges_list1)))
    # print(tmpidx)
    weightveclist = append(weightveclist, list(weightvec))
    if(length(edges_list2[weightvec[tmpidx] != 0]) == 0) break
    # print(weightvec[tmpidx])

    edges_list1 = append(edges_list1[[1]], edges_list2[weightvec[tmpidx] != 0])
    # edges_list1 = list(append(edges_list1[[1]], lapply(edges_list2[weightvec[tmpidx] != 0], unlist)))
    print(edges_list1)

    print(lapply(edges_list1[[1]], function(x) names(Y)[x]))
    # print(edges_list2[weightvec[tmpidx] != 0])
    edges_list2 = edges_list2[weightvec[tmpidx] == 0]
  }

  weightmat = matrix(0, nr = p, nc = p)
  weightmat[do.call("rbind", edges_list1[[1]])] <- 1
  weightmat[lower.tri(weightmat)] <- t(weightmat)[lower.tri(weightmat)]

  dp <- list(weightveclist = weightveclist, weightmat = weightmat, #mat = mat, grid_tmp = grid_tmp,
             edges_list = edges_list, edges_list1 = edges_list1, edges_list2 = edges_list2,
             supp = supp, n = n, p = p, marginalProbmat = marginalProbmat,
             namesY = namesY)
  class(dp) <- "doublep"
  return(dp)
}
