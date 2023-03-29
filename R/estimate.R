#' Estimate the joint probabilities
#'
#' @param efi.obj An object produced by \code{efi}.
#' @param expr A string that contains formulas with the desired probabilities.
#' @return Estimates and Standard error of the parameters.
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
#' EFI = efi(Y, dp, freq = F)
#' estimate(EFI, "(Hair == \"Black\") & (Eye == \"Brown\")")
#' estimate(EFI, "(Hair == \"Black\") & (Sex == \"Male\")")
#' @export

estimate = function(efi.obj, expr){
  call <- match.call()
  imp = efi.obj$imp
  n = efi.obj$n
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
