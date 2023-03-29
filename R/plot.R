#' plot class doubledp
#'
#' @param dp a class doublep computed from \code{doublep}
#' @rdname plot
#' @method plot doublep
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
#' plot(dp)
#' @export
plot.doublep = function(dp){
  weightmat = dp$weightmat
  namesY = dp$namesY
  g <- graph_from_adjacency_matrix(weightmat, mode = "lower", weighted = "weight")
  # plot(g, edge.width = E(g)$weight, edge.label = E(g)$weight, label = coords)

  V(g)$color <- "#C83200"
  plot(g, edge.label = E(g)$weight, vertex.label = namesY, vertex.label.dist = 2.5)
}
