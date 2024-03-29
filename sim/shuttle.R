require(R.utils)

library(igraph)
library(readr)
library(tidyr)
library(dplyr)
library(reshape2)
library(doParallel)

library(vcd)

library(mice)
library(missForest)
library(Amelia)
library(EFI)

library(ggplot2)
library(ggpubr)

set.seed(123)
SIMNUM = 100

high_dim = F
MAR = T


data_name = "shuttle"

if(high_dim == F){
  data_list = c("UCBAdmissions", "Titanic", "esoph", "PreSex", "HairEyeColor")
}else{
  data_list = c("spect_data", "mushroom_data", "promoters_data", "chess_data")
}

timenow1 = Sys.time()
timenow0 = gsub(' ', '_', gsub('[-:]', '', timenow1))
timenow = paste(timenow0, ".txt", sep = "")

dir.create(timenow0)
setwd(timenow0)

cores = min(detectCores() - 2, 100)
myCluster <- makeCluster(cores, # number of cores to use
                         type = "PSOCK",
                         outfile = timenow) # type of cluster
registerDoParallel(myCluster)

if(!interactive()) sink(timenow, append=TRUE)

# spect_names <- c("OVERALL_DIAGNOSIS", "F1R", "F1S", "F2R", "F2S", "F3R", "F3S", "F4R", "F4S", "F5R", "F5S", "F6R", "F6S", "F7R", "F7S", "F8R", "F8S", "F9R", "F9S", "F10R", "F10S", "F11R", "F11S", "F12R", "F12S", "F13R", "F13S", "F14R", "F14S", "F15R", "F15S", "F16R", "F16S", "F17R", "F17S", "F18R", "F18S")
# spect_data_train <- read_csv("https://archive.ics.uci.edu/ml/machine-learning-databases/spect/SPECT.train", col_names = spect_names,
#                              col_types = cols(.default = "f"))
# spect_data_test <- read_csv("https://archive.ics.uci.edu/ml/machine-learning-databases/spect/SPECT.test", col_names = spect_names,
#                             col_types = cols(.default = "f"))
# spect_data = rbind(spect_data_train, spect_data_test)
#
# promoters_data <- read_csv("https://archive.ics.uci.edu/ml/machine-learning-databases/molecular-biology/promoter-gene-sequences/promoters.data", col_names=FALSE)
# colnames(promoters_data) <- c("Class", "id", "Sequence")
# promoters_data <- separate(promoters_data, Sequence, into = paste0("pos_", 1:57), sep = "")
# promoters_data <- promoters_data %>% dplyr::select(-id, -pos_1) %>% mutate_all(factor)
# promoters_data <- promoters_data[,1:15]
#
# mushroom_names <- c("class", "cap.shape", "cap.surface", "cap.color", "bruises", "odor",
#                     "gill.attachment", "gill.spacing", "gill.size", "gill.color",
#                     "stalk.shape", "stalk.root", "stalk.surface.above.ring",
#                     "stalk.surface.below.ring", "stalk.color.above.ring",
#                     "stalk.color.below.ring", "veil.type", "veil.color", "ring.number",
#                     "ring.type", "spore.print.color", "population", "habitat")
# mushroom_data <- read_csv("http://archive.ics.uci.edu/ml/machine-learning-databases/mushroom/agaricus-lepiota.data", col_names = mushroom_names,
#                           col_types = cols(.default = "f"))
# # lapply(mushroom_data, levels)
#
# column_names <- c("age_gt_60", "air", "airBoneGap", "ar_c", "ar_u", "bone", "boneAbnormal", "bser",
#                   "history_buzzing", "history_dizziness", "history_fluctuating", "history_fullness",
#                   "history_heredity", "history_nausea", "history_noise", "history_recruitment",
#                   "history_ringing", "history_roaring", "history_vomiting", "late_wave_poor", "m_at_2k",
#                   "m_cond_lt_1k", "m_gt_1k", "m_m_gt_2k", "m_m_sn", "m_m_sn_gt_1k", "m_m_sn_gt_2k",
#                   "m_m_sn_gt_500", "m_p_sn_gt_2k", "m_s_gt_500", "m_s_sn", "m_s_sn_gt_1k", "m_s_sn_gt_2k",
#                   "m_s_sn_gt_3k", "m_s_sn_gt_4k", "m_sn_2_3k", "m_sn_gt_1k", "m_sn_gt_2k", "m_sn_gt_3k",
#                   "m_sn_gt_4k", "m_sn_gt_500", "m_sn_gt_6k", "m_sn_lt_1k", "m_sn_lt_2k", "m_sn_lt_3k",
#                   "middle_wave_poor", "mod_gt_4k", "mod_mixed", "mod_s_mixed", "mod_s_sn_gt_500",
#                   "mod_sn", "mod_sn_gt_1k", "mod_sn_gt_2k", "mod_sn_gt_3k", "mod_sn_gt_4k", "mod_sn_gt_500",
#                   "notch_4k", "notch_at_4k", "o_ar_c", "o_ar_u", "s_sn_gt_1k", "s_sn_gt_2k", "s_sn_gt_4k",
#                   "speech", "static_normal", "tymp", "viith_nerve_signs", "wave_V_delayed", "waveform_ItoV_prolonged",
#                   "identifier", "class")
# audiology <- read_csv("http://archive.ics.uci.edu/ml/machine-learning-databases/audiology/audiology.standardized.data",
#                       col_types = cols(.default = "f"), col_names = column_names)
# audiology <- (audiology %>% dplyr::select(-identifier))
# audiology <- audiology[c(ncol(audiology), 1:(ncol(audiology)-1))]
# # lapply(audiology, levels)
#
#
# chess_data <- read_csv("http://archive.ics.uci.edu/ml/machine-learning-databases/chess/king-rook-vs-king-pawn/kr-vs-kp.data",
#                        col_types = cols(.default = "f"))
# colnames(chess_data) <- c("bkblk", "bknwy", "bkon8", "bkona", "bkspr", "bkxbq", "bkxcr", "bkxwp", "blxwp", "bxqsq",
#                           "cntxt", "dsopp", "dwipd", "hdchk", "katri", "mulch", "qxmsq", "r2ar8", "reskd", "reskr",
#                           "rimmx", "rkxwp", "rxmsq", "simpl", "skach", "skewr", "skrxp", "spcop", "stlmt", "thrsk",
#                           "wkcti", "wkna8", "wknck", "wkovl", "wkpos", "wtoeg", "won")
# chess_data <- chess_data[c(ncol(chess_data), 1:(ncol(chess_data)-1))]
#
# car <- read_csv("https://archive.ics.uci.edu/ml/machine-learning-databases/car/car.data",
#                 col_names=c("buying", "maint", "doors", "persons", "lug_boot", "safety", "class_values"),
#                 col_types = cols(.default = "f"))
# colnames(car) <- c("Buying", "Maintenance", "Doors", "Persons", "Lug_Boot", "Safety", "Class_Values")
# car <- car[c(ncol(car), 1:(ncol(car)-1))]

data(HairEyeColor, Titanic, UCBAdmissions, esoph)
data(shuttle, package = "MASS")

HairEyeColor <- as.data.frame.table(HairEyeColor) %>% uncount(Freq)
Titanic <- as.data.frame.table(Titanic) %>% uncount(Freq)
Titanic <- Titanic[c(ncol(Titanic), 1:(ncol(Titanic)-1))]
UCBAdmissions <- as.data.frame.table(UCBAdmissions) %>% uncount(Freq)

esoph <- melt(esoph)
esoph <- esoph %>% uncount(value)
names(esoph)[ncol(esoph)] <- "cases"
esoph <- esoph[c(ncol(esoph), 1:(ncol(esoph)-1))]

data(PreSex, package = "vcd")
PreSex <- as.data.frame.table(PreSex) %>% uncount(Freq)

eval(parse(text = paste("df_true =", data_name)))
df_true = df_true %>%
  mutate_if(is.ordered, factor, ordered = F)

df_true <- df_true[,sapply(df_true, nlevels) > 1]
if(high_dim == T) df_true <- df_true[,sapply(df_true, nlevels) < 6]
modes = sapply(df_true, function(x) names(which.max(table(x))))

if(high_dim == T) df_true <- df_true[,mapply(function(x, y) sum(x == y), df_true, modes) / nrow(df_true) < 0.9]

# df_true = car
p = ncol(df_true)
n = nrow(df_true)
# for(k in 1:p){
#   df_true[[k]] <- as.factor(df_true[[k]])
# }
modes = sapply(df_true, function(x) names(which.max(table(x))))

levelone = names(which.max(summary(df_true[[1]])))
trueval = mean(df_true[[1]] == levelone, na.rm = T)

mis_rate_vec = seq(from = 0.1, to = 0.9, length = 9)
# mis_rate_vec = 0.5

oper <-
  foreach(mis_rate_idx= 1:length(mis_rate_vec), .packages = c("R.utils", "CVXR", "mice", "missForest", "Amelia", "EFI","dplyr", "igraph")) %:%
  foreach(simnum = 1:SIMNUM, .combine='rbind', .packages = c("R.utils", "CVXR", "mice", "missForest", "Amelia", "EFI", "dplyr", "igraph")) %dopar% {

mis_rate = mis_rate_vec[mis_rate_idx]
obs_rate = 1 - mis_rate
print(mis_rate)

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

# if(high_dim){
#   df <- df %>% select(order(cor(sweep(df[complete.cases(df),], 2, modes, "=="))[1,], decreasing = T)[1:6])
#   nacols <- 1
# }

if(high_dim == T){
  # varidx =  combn(p, 2)
  varidx =  combn(p, 2)[,1:(p-1)]

  edges_list = list()
  itmp = 1
  while(T){
    tmplist <- plyr::alply(varidx[,sample(1:ncol(varidx), round(sqrt(p)))],2,c)
    attributes(tmplist) <- NULL
    if((1 %in% unique(unlist(tmplist)))){
      edges_list[[itmp]] <- tmplist
      itmp = itmp + 1
    }
    if(itmp > 100) break
  }

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

dp <- tryCatch(withTimeout({doublep(df, edges_list, R = 1)}, timeout = 60 * 10),
               error = function(e) {
                 cat("An error occurred in dp: ", conditionMessage(e))
               })


do.call("rbind", unlist(dp$edges_list1, recursive = F))
}

names(oper) <- mis_rate_vec

save.image(paste(timenow0, ".RData", sep = ""))

for(misidx in 1:length(mis_rate_vec)){
  misrate =mis_rate_vec[misidx]
  matrix_data = oper[[misidx]]
  # matrix_data = do.call("rbind", unlist(dp$edges_list1, recursive = F))
  plot_data = data.frame(matrix_data) %>% count(X1, X2)
  plot_data$n = plot_data$n / (nrow(oper[[1]]) / 3)
  wmat = matrix(0, ncol = p, nrow = p)
  wmat[as.matrix(plot_data[,1:2])] <- plot_data[,3]
  g <- igraph::graph_from_adjacency_matrix(wmat, mode = "upper", weighted = "weight")
  E(g)$weight
  # g <- igraph::graph_from_edgelist(do.call("rbind", unique(unlist(dp$edges_list1, recursive = F))), directed = F)
  # g <- igraph::graph_from_edgelist(plot_data, directed = F)

  igraph::V(g)$color <- "#C83200"

  png(filename=paste(misrate, "_shuttle.png", sep = ""))
  plot.igraph(g, vertex.label = names(df), vertex.label.dist = 2.5, layout = layout_as_star(g), edge.width = E(g)$weight * 10, edge.label = round(E(g)$weight, 2))
  dev.off()
}
install.packages("igraph")
