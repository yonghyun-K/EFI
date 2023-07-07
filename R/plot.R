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
  namesY = dp$namesY

  sample_el <- do.call("rbind", unlist(dp$edges_list1, recursive = F))
  sample_el2 <- do.call("rbind", unique(unlist(dp$edges_list, recursive = F)))
  sample_el2 <- as.matrix(dplyr::anti_join(data.frame(sample_el2), data.frame(sample_el), by = dplyr::join_by(X1, X2)))

  g <- igraph::graph_from_edgelist(rbind(sample_el, sample_el2), directed=F)
  edgesintree = sapply(dp$edges_list1, length)
  tmp <- c(rep(1:length(edgesintree), times = edgesintree), rep(0, nrow(sample_el2)))

  # g2 <- graph_from_edgelist(rbind(sample_el), directed=F)

  igraph::plot.igraph(g, layout = igraph::layout_nicely(g), vertex.label = namesY, vertex.shape = "none", edge.width = ifelse(tmp, 10, 1),
       edge.color = ifelse(tmp == 0, "grey",tmp), asp = 0.5, vertex.label.cex = 1.2)
}
