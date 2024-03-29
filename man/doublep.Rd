% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doublep.R
\name{doublep}
\alias{doublep}
\title{Double project data from candidate models}
\usage{
doublep(Y, edges_list, freq = F, R = 5)
}
\arguments{
\item{Y}{categorical data with missingnes.Missing valuesa are encoded as NA. If Freq = T, the last column serves as the integer frequencies of the grouped observations.}

\item{freq}{an optional logical value that indicates whether frequency column exists in the data frame Y.}

\item{cand.edges}{candidate edges that are considered in the candidate graphical exponential model. Typicially, as.list(data.frame(combn(ncol(Y), 2))) is used.}
}
\value{
A class  \code{"doubledp"}
}
\description{
Double project data from candidate models
}
\examples{
data(HairEyeColor)
p = 3
Y = do.call("rbind", apply(as.data.frame.table(HairEyeColor), 1, function(x) matrix(rep(x[1:p], each = x[p+1]), nc = p)))
Y = as.data.frame(Y)
for(k in 1:p){
Y[[k]] = factor(Y[[k]])
}
names(Y) <- names(dimnames(HairEyeColor))
(n = nrow(Y)); sum(HairEyeColor)
delta = matrix(rbinom(n * p, 1, 0.5), nr = n, nc = p)
Y[delta == 0] = NA

cand.edges = as.list(data.frame(combn(p, 2)))
dp = doublep(Y, cand.edges, freq = F)
}
