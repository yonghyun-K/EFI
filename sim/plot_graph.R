high_dim = T
MAR = F

data_name = "congressional"

eval(parse(text = paste("df_true =", data_name)))
df_true = df_true %>%
  mutate_if(is.ordered, factor, ordered = F)

df_true <- df_true[,sapply(df_true, nlevels) > 1]
if(high_dim == T) df_true <- df_true[,sapply(df_true, nlevels) < 6]
modes = sapply(df_true, function(x) names(which.max(table(x))))

if(high_dim == T) df_true <- df_true[,mapply(function(x, y) sum(x == y), df_true, modes) / nrow(df_true) < 0.9]

modes = sapply(df_true, function(x) names(which.max(table(x))))

# df_true = car
p = ncol(df_true)
n = nrow(df_true)
# for(k in 1:p){
#   df_true[[k]] <- as.factor(df_true[[k]])
# }
# modes = sapply(df_true, function(x) names(which.max(table(x))))

levelone = names(which.max(summary(df_true[[1]])))
trueval = mean(df_true[[1]] == levelone, na.rm = T)

mis_rate_vec = seq(from = 0.1, to = 0.9, length = 9)

mis_rate = 0
# mis_rate = 0.9
obs_rate = 1 - mis_rate
# print(mis_rate)

cnt = 0
while(T){
  cnt = cnt + 1
  print(cnt)
  df <- df_true

  if(MAR == F){
    delta = sapply(1:p, function(x) rbinom(n, 1, obs_rate) == 1)
  }else{
    # MAR
    p2 = round(p / 2)
    delta = apply(df_true, 1, function(y){
      deltav = integer(p)
      tmpval = sample(1:p, 1)
      tmpidx = (tmpval:(tmpval + p2 - 1)) %% p + 1
      deltav[tmpidx] = 1
      if(tmpval %% 2 == 0){
        deltav[-tmpidx] = rbinom(p2,1,ifelse(y[tmpidx] == modes[tmpidx], min(c(3 * obs_rate, 1)), max(c(obs_rate / 3, 0))))[1:(p - p2)]
      }else{
        deltav[-tmpidx] = rbinom(p2,1,ifelse(y[tmpidx] == modes[tmpidx], max(c(obs_rate / 3, 0)), min(c(3 * obs_rate, 1))))[1:(p - p2)]
      }
      # deltav[-tmpidx] = rbinom(p2,1,ifelse(y[tmpidx] == 1, 0.1, 0.9))

      deltav
    })
    delta = t(delta)
  }

  # if(high_dim == F & MAR == F){
  #   # MCAR
  #   delta = sapply(1:p, function(x) rbinom(n, 1, obs_rate) == 1)
  # }else if(high_dim == T & MAR == F){
  #   # MCAR - high dim
  #   delta = matrix(1, nr =  n, nc = p)
  #   delta[,1] = (rbinom(n, 1, obs_rate) == 1)
  # }else if(high_dim == F & MAR == T){
  #   # MAR
  #   p2 = round(p / 2)
  #   delta = apply(df_true, 1, function(y){
  #     deltav = integer(p)
  #     tmpval = sample(1:p, 1)
  #     tmpidx = (tmpval:(tmpval + p2 - 1)) %% p + 1
  #     deltav[tmpidx] = 1
  #     if(tmpval %% 2 == 0){
  #       deltav[-tmpidx] = rbinom(p2,1,ifelse(y[tmpidx] == modes[tmpidx], min(c(3 * obs_rate, 1)), max(c(obs_rate / 3, 0))))[1:(p - p2)]
  #     }else{
  #       deltav[-tmpidx] = rbinom(p2,1,ifelse(y[tmpidx] == modes[tmpidx], max(c(obs_rate / 3, 0)), min(c(3 * obs_rate, 1))))[1:(p - p2)]
  #     }
  #     # deltav[-tmpidx] = rbinom(p2,1,ifelse(y[tmpidx] == 1, 0.1, 0.9))
  #
  #     deltav
  #   })
  #   delta = t(delta)
  # }else if(high_dim == T & MAR == T){
  #   # MAR - high dim
  #   tmpidx = order(abs(cor(sweep(df_true, 2, modes, "=="))[1,]), decreasing = T)[2:5]
  #   delta = apply(df_true, 1, function(y){
  #     deltav = rep(1, p)
  #     deltav[1] = rbinom(1,1,ifelse(y[sample(tmpidx, 1)] == modes[sample(tmpidx, 1)], min(c(3 * obs_rate, 1)), max(c(obs_rate / 3, 0))))
  #     deltav
  #   })
  #   delta = t(delta)
  # }

  df[delta == 0] <- NA

  # print(sapply(df, function(x) sum(!is.na(unique(x)))))
  if(all(sapply(df, function(x) sum(!is.na(unique(x)))) != 1) | cnt > 100) break
}
nacols = sapply(df, function(x) any(is.na(x)))
print("Missing values are generated")

methods_mice = c("polyreg", "pmm", "cart", "lda")
vals = c()
accs = c()
tmes = c()
for(method in methods_mice){
  time_perf1 = Sys.time()
  imp_mice <- tryCatch(mice(df, m = 5, seed = 123, printFlag = F, method = method),
                       error = function(e) {
                         cat("An error occurred in MICE: ", conditionMessage(e))
                       })
  if(!is.null(imp_mice)){
    df_mice = complete(imp_mice,1:5)
    valmice = mean(df_mice[[1]] == levelone); names(valmice) <- method
    accmice = mean(complete(imp_mice)[is.na(df[,1]),1, drop = F] == df_true[is.na(df[,1]),1, drop = F]); names(accmice) <- method
    vals = c(vals, valmice)
    accs = c(accs, accmice)
  }else{
    vals = c(vals, NA)
    accs = c(accs, NA)
  }
  time_perf2 = Sys.time()
  timetmp = time_perf2 - time_perf1; names(timetmp) = method
  tmes = c(tmes, timetmp)
}
print("MICE successfully done")

time_perf1 = Sys.time()
imp_missF <- tryCatch(withTimeout({missForest(cbind(1, df))}, timeout = 60 * 10),
                      error = function(e) {
                        cat("An error occurred in missF: ", conditionMessage(e))
                      })
if(!is.null(imp_missF)) df_missF = imp_missF$ximp[-1]
time_perf2 = Sys.time()

vals = c(vals, missF = ifelse(!is.null(imp_missF), mean(imp_missF$ximp[-1][[1]] == levelone), NA))
accs = c(accs, missF = ifelse(!is.null(imp_missF), mean(imp_missF$ximp[-1][is.na(df[,1]), 1, drop = F] == df_true[is.na(df[,1]),1, drop = F]), NA))
tmes = c(tmes, time_perf2 - time_perf1)

# if(high_dim == F){
#   imp_Amelia <- tryCatch(withTimeout({amelia(df, noms = names(df), p2s = 0)}, timeout = 60 * 10),
#                          error = function(e) {
#                            cat("An error occurred in Amelia: ", conditionMessage(e))
#                          })
#   if(!is.null(imp_Amelia)) df_Amelia = do.call("rbind", imp_Amelia$imputations)
#
#   vals = c(vals, Amelia = ifelse(!is.null(imp_Amelia), mean(df_Amelia[[1]] == levelone), NA))
#   accs = c(accs, Amelia = ifelse(!is.null(imp_Amelia), mean(imp_Amelia$imputations[[1]][is.na(df[,1]), 1, drop = F] == df_true[is.na(df[,1]),1, drop = F]), NA))
# }
time_perf1 = Sys.time()
complete_CC = df
for (k in 1:p) {
  complete_CC[, k][is.na(complete_CC[, k])] <- modes[k]
}
time_perf2 = Sys.time()
tmes = c(tmes, mode = time_perf2 - time_perf1)

modeacc = mean(complete_CC[is.na(df[,1]),1, drop = F] == df_true[is.na(df[,1]),1, drop = F])
modeval = mean(df[[1]] == levelone, na.rm = T)


# if(high_dim){
#   df <- df %>% select(order(cor(sweep(df[complete.cases(df),], 2, modes, "=="))[1,], decreasing = T)[1:6])
#   nacols <- 1
# }


if(high_dim == T){
  varidx =  combn(p, 2)
  # varidx =  combn(p, 2)[,1:(p-1)]

  edges_list = list()
  itmp = 1
  while(T){
    tmplist <- plyr::alply(varidx[,sample(1:ncol(varidx), round(sqrt(p)))],2,c)
    attributes(tmplist) <- NULL
    if((1 %in% unique(unlist(tmplist)))){
      edges_list[[itmp]] <- tmplist
      itmp = itmp + 1
    }
    if(itmp > 300) break
  }
  edges_list <- edges_list[!duplicated(purrr::map(edges_list, function(x) unique(sort(unlist(x)))))]
}else{
  # edges_list <- apply(combn(ncol(df), 2), 2, list)

  varidx =  combn(p, 2)
  tmpidx = combn(ncol(varidx), round(sqrt(p)))
  edges_list = list()

  cnt = 0
  for (tmp in 1:ncol(tmpidx)){
    tmplist <- plyr::alply(varidx[,tmpidx[,tmp]],2,c)
    attributes(tmplist) <- NULL
    g <- igraph::graph_from_edgelist(do.call("rbind", tmplist), directed = F)
    if(!igraph::has_eulerian_cycle(g)){
      cnt = cnt + 1
      edges_list[[cnt]] <- tmplist
    }
  }
}



time_perf1 = Sys.time()
dp <- tryCatch(withTimeout({doublep(df, edges_list, R = 1)}, timeout = 60 * 10),
               error = function(e) {
                 cat("An error occurred in dp: ", conditionMessage(e))
               })
# dp$weightveclist
# dp$edges_list1

# order(abs(cor(df_true[,1] == 1, df_true == 1)), decreasing = T)

imp_EFI <- tryCatch(withTimeout({efi(df, dp)}, timeout = 60 * 10),
                    error = function(e) {
                      cat("An error occurred in efi: ", conditionMessage(e))
                    })
time_perf2 = Sys.time()

tmes = c(tmes, EFI = time_perf2 - time_perf1)

if(!is.null(imp_EFI)){
  complete_EFI = imp_EFI$imp %>% group_by(id) %>% summarize(maxw = max(w)) %>%
    left_join(imp_EFI$imp, by = c("id", "maxw" = "w"), multiple = "first") %>%
    select(-maxw, -Freq, -id)
  # complete_EFI = imp_EFI$imp %>% group_by(id, .[[1]]) %>% summarize(sumw = sum(w)) %>%
  #   filter(sumw == max(sumw))  %>% ungroup() %>%
  #   left_join(imp_EFI$imp, by = c("id"), multiple = "first") %>%
  #   select(-sumw, -Freq, -".[[1]]", -id, -w)
  EFIval = sum((imp_EFI$imp[[1]] == levelone) * ( imp_EFI$imp$w)) / n
  EFIacc = mean(complete_EFI[is.na(df[,1]),1, drop = F] == df_true[is.na(df[,1]),1, drop = F])
}else{
  EFIval = NA
  EFIacc = NA
}

result1 = c(EFI = EFIval,
            mode = modeval,
            # estimate(imp_EFI, "(Class_Values == \"unacc\")"),
            vals
)
result2 = c(EFI = EFIacc,
            mode = modeacc,
            accs)
sort(abs(result1 - trueval))
result2

# dp$edges_list
# dp$weightveclist
#
wmat = matrix(0, ncol = p, nrow = p)
wmat[do.call("rbind", unlist(dp$edges_list, recursive = F))] <- 1
wmat[do.call("rbind", unlist(dp$edges_list1, recursive = F))] <- 10

g <- igraph::graph_from_adjacency_matrix(wmat, mode = "upper", weighted = "weight")
E(g)$weight
# g <- igraph::graph_from_edgelist(do.call("rbind", unique(unlist(dp$edges_list1, recursive = F))), directed = F)
# g <- igraph::graph_from_edgelist(plot_data, directed = F)

igraph::V(g)$color <- "#C83200"
par(mar = c(0, 0, 0, 0))  # Set all margin values to 0
plot.igraph(g, vertex.label = names(df), vertex.shape = "none", layout = layout_nicely(g), edge.width = E(g)$weight,
            edge.color = ifelse(E(g)$weight == 1, "grey", "orange"), asp = 0.5, vertex.label.cex = 1)
# EFI:::plot.doublep(dp)



# A graph from sample data
sample_el <- do.call("rbind", unlist(dp$edges_list1, recursive = F))
sample_el2 <- do.call("rbind", unique(unlist(dp$edges_list, recursive = F)))
sample_el2 <- as.matrix(anti_join(data.frame(sample_el2), data.frame(sample_el)))

g <- graph_from_edgelist(rbind(sample_el, sample_el2), directed=F)
tmp <- c(rep(1:(nrow(sample_el) / 4), each = 4), rep(0, nrow(sample_el2)))
ifelse(tmp == 0, "grey",tmp)

g2 <- graph_from_edgelist(rbind(sample_el), directed=F)

par(mar = c(0, 3, 0, 3))  # Set all margin values to 0
# Always plot graphs with this same layout
plot(g, layout = layout_nicely(g), vertex.label = names(df), vertex.shape = "none", edge.width = ifelse(tmp, 10, 1),
     edge.color = ifelse(tmp == 0, "grey",tmp), asp = 0.5, vertex.label.cex = 1.2)

