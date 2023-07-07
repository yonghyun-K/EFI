#' Construct Tree
#'
#' Constructs a tree based on the given parameters.
#'
#' @param p Number of variables.
#' @param high_dim Logical value indicating whether the tree should be constructed in high dimension or not.
#' @param ntree Number of trees to be constructed.
#'
#' @return A list of edges representing the constructed trees.
#'
#' @details This function constructs a tree based on the specified number of variables (\code{p}),
#' whether it should be constructed in high dimension (\code{high_dim}), and the number of trees to be constructed (\code{ntree}).
#' The function returns a list of edges representing the constructed trees.
#'
#' @examples
#' construct_tree(10, TRUE, 100)
#'
#' construct_tree(5)
#'
#' @seealso \code{\link{combn}}, \code{\link{plyr::alply}}, \code{\link{igraph::graph_from_edgelist}},
#' \code{\link{igraph::has_eulerian_cycle}}
#'
#' @importFrom plyr alply
#' @importFrom igraph graph_from_edgelist has_eulerian_cycle
#' @export

construct_tree = function(p, high_dim = ifelse(p > 6, T, F), ntree = 100){
  if(p == 2) return(list(list(c(1,2))))
  varidx =  combn(p, 2)
  edges_list = list()
  if(high_dim == F){
    tmpidx = combn(ncol(varidx), round(sqrt(p)))

    cnt = 0
    for (tmp in 1:ncol(tmpidx)){
      tmplist <- plyr::alply(varidx[,tmpidx[,tmp]],2,c)
      attributes(tmplist) <- NULL
      cnt = cnt + 1
      edges_list[[cnt]] <- tmplist
    }
  }else{
    itmp = 1; cnt = 0
    while(T){
      cnt = cnt + 1
      tmplist <- plyr::alply(varidx[,sample(1:ncol(varidx), round(sqrt(p)))],2,c)
      attributes(tmplist) <- NULL
      tmplist <- tmplist[order(sapply(tmplist, function(x) p * (x[1]- 1) + (x[2] - 1)))]
      if(!(list(tmplist) %in% edges_list) &
         !igraph::has_eulerian_cycle(igraph::graph_from_edgelist(
           do.call("rbind", tmplist), directed = F))){
        edges_list[[itmp]] <- tmplist
        itmp = itmp + 1
      }
      if(itmp > ntree) break
      if(cnt > 1000){
        warning("cnt > 1000: ntree would be smaller than the specified number."); break
      }
    }
  }
  return(edges_list)
}
