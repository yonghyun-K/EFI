require(R.utils)

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
library(igraph)

library(ggplot2)
library(ggpubr)

set.seed(123)
SIMNUM = 100

high_dim = T
MAR = T

if(interactive()) data_name = "breast_cancer"

if(high_dim == F){
  data_list = c("esoph", "PreSex", "HairEyeColor", "UCBAdmissions", "Titanic")
}else{
  # data_list = c("spect_data", "mushroom_data", "promoters_data", "chess_data")
  data_list = c("breast_cancer", "credit", "car", "congressional")
  # data_list = c("endgame")
}

timenow1 = Sys.time()
timenow0 = gsub(' ', '_', gsub('[-:]', '', timenow1))
timenow = paste(timenow0, ".txt", sep = "")

dir.create(timenow0)
setwd(timenow0)

cores = min(detectCores() - 3, 100)
myCluster <- makeCluster(cores, # number of cores to use
                         type = "PSOCK") # type of cluster
                         # outfile = timenow)
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
#          "cntxt", "dsopp", "dwipd", "hdchk", "katri", "mulch", "qxmsq", "r2ar8", "reskd", "reskr",
#          "rimmx", "rkxwp", "rxmsq", "simpl", "skach", "skewr", "skrxp", "spcop", "stlmt", "thrsk",
#          "wkcti", "wkna8", "wknck", "wkovl", "wkpos", "wtoeg", "won")
# chess_data <- chess_data[c(ncol(chess_data), 1:(ncol(chess_data)-1))]

car <- read_csv("https://archive.ics.uci.edu/ml/machine-learning-databases/car/car.data",
                     col_names=c("buying", "maint", "doors", "persons", "lug_boot", "safety", "class_values"),
                     col_types = cols(.default = "f"))
colnames(car) <- c("Buying", "Maintenance", "Doors", "Persons", "Lug_Boot", "Safety", "Class_Values")
car <- car[c(ncol(car), 1:(ncol(car)-1))]
summary(car)

col_names <- c("Class", "Age", "Menopause", "Tumor_Size", "Inv_Nodes", "Node_Caps",
               "Deg_Malig", "Breast", "Breast_Quad", "Irradiat")

# Import the data from the URL and assign column names
breast_cancer <- read_csv("https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer/breast-cancer.data", col_names = col_names,
                               col_types = cols(.default = "f"), na = "?")
breast_cancer <- breast_cancer[complete.cases(breast_cancer),]

# Define the collapsing ranges and labels
ranges <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54")
labels <- c("0-9", "0-9", "10-19", "10-19", "20-29", "20-29", "30-39", "30-39", "40-54", "40-54", "40-54")

# Collapse the levels of the Tumor_Size column
breast_cancer <- breast_cancer %>%
  mutate(Tumor_Size = recode(Tumor_Size, !!!setNames(labels, ranges)))

# Define the collapsing ranges and labels
ranges <- c("0-2", "3-5", "6-8", "9-11", "12-14", "15-17", "18-20", "21-23", "24-26")
labels <- c("0-8", "0-8", "0-8", "9-26", "9-26", "9-26", "9-26", "9-26", "9-26")

# Collapse the levels of the Inv_Nodes variable
breast_cancer <- breast_cancer %>%
  mutate(Inv_Nodes = recode(Inv_Nodes, !!!setNames(labels, ranges)))

# Define the collapsing ranges and labels
ranges <- c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")
labels <- c("20-39", "20-39", "40-49", "50-59", "60-69", "70-79")

# Collapse the levels of the Inv_Nodes variable
breast_cancer <- breast_cancer %>%
  mutate(Age = recode(Age, !!!setNames(labels, ranges)))


summary(breast_cancer)


col_names <- c(
  "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10",
  "A11", "A12", "A13", "A14", "A15", "A16"
)

column_types <- cols(
  A1 = col_factor(),
  A2 = col_double(),
  A3 = col_double(),
  A4 = col_factor(),
  A5 = col_factor(),
  A6 = col_factor(),
  A7 = col_factor(),
  A8 = col_double(),
  A9 = col_factor(),
  A10 = col_factor(),
  A11 = col_double(),
  A12 = col_factor(),
  A13 = col_factor(),
  A14 = col_double(),
  A15 = col_double(),
  A16 = col_factor()
)

# Import the data from the URL and assign column names
credit <- read_csv("https://archive.ics.uci.edu/ml/machine-learning-databases/credit-screening/crx.data", col_names = col_names,
                        col_types = column_types, na = "?")
credit <- credit[complete.cases(credit),]

# Specify the columns to discretize
columns_to_discretize <- c("A2", "A3", "A8", "A11", "A14", "A15")

# Discretize and convert to factor
for (column in columns_to_discretize) {
  credit[[column]] <- log(credit[[column]] + 2)
  credit[[column]] <- cut(credit[[column]], breaks = 4)  # Adjust the number of breaks as needed
  credit[[column]] <- as.factor(credit[[column]])
}

credit = credit %>% filter(A4 != "l" & A5 != "gg" & A13 != "p")
credit = credit %>% select(-A6, -A7)
credit = droplevels(credit)
credit <- credit[c(ncol(credit), 1:(ncol(credit)-1))]

# Define the collapsing ranges and labels
ranges <- c("(0.69,1.58]", "(1.58,2.46]", "(2.46,3.35]", "(3.35,4.24]")
labels <- c("(0.69,1.58]", "(1.58,2.46]", "(2.46,4.24]", "(2.46,4.24]")

# Collapse the levels of the Inv_Nodes variable
credit <- credit %>%
  mutate(A11 = recode(A11, !!!setNames(labels, ranges)))
summary(credit)

col_names <- col_names <- c("party", "handicapped_infants", "water_project_cost_sharing", "adoption_of_the_budget_resolution",
                            "physician_fee_freeze", "el_salvador_aid", "religious_groups_in_schools", "anti_satellite_test_ban",
                            "aid_to_nicaraguan_contras", "mx_missile", "immigration", "synfuels_corporation_cutback",
                            "education_spending", "superfund_right_to_sue", "crime", "duty_free_exports", "export_administration_act_south_africa")

# Import the data from the URL and assign column names
congressional <- read_csv("https://archive.ics.uci.edu/ml/machine-learning-databases/voting-records/house-votes-84.data", col_names = col_names,
                               col_types = cols(.default = "f"), na = "?")
congressional <- congressional[complete.cases(congressional),]
summary(congressional)

col_names <- c("top_left_square", "top_middle_square", "top_right_square",
               "middle_left_square", "middle_middle_square", "middle_right_square",
               "bottom_left_square", "bottom_middle_square", "bottom_right_square",
               "class")

# Import the data from the URL and assign column names
endgame <- read_csv("https://archive.ics.uci.edu/ml/machine-learning-databases/tic-tac-toe/tic-tac-toe.data", col_names = col_names,
                    col_types = cols(.default = "f"), na = "?")
endgame <- endgame[complete.cases(endgame),]
endgame <- endgame[c(ncol(endgame), 1:(ncol(endgame)-1))]
summary(endgame)


data(HairEyeColor, Titanic, UCBAdmissions, esoph)

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

# data(UKSoccer, package = "vcd")
# as.data.frame.table(UKSoccer)

res_list <- NULL
res2_list <- NULL
res3_list <- NULL
res4_list <- NULL
for(data_name in data_list){
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

oper <-
  foreach(mis_rate_idx= 1:length(mis_rate_vec), .packages = c("R.utils", "CVXR", "mice", "missForest", "Amelia", "EFI","dplyr", "igraph")) %:%
  foreach(simnum = 1:SIMNUM, .combine='cbind', .packages = c("R.utils", "CVXR", "mice", "missForest", "Amelia", "EFI", "dplyr", "igraph")) %dopar% {

# mis_rate_idx = 4
#       for(simnum in 1:SIMNUM){

    print(paste("mis_rate_idx = ", mis_rate_idx))
    print(paste("simnum = ", simnum))
    mis_rate = mis_rate_vec[mis_rate_idx]
    # mis_rate = 0.2
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
    tmes = c(tmes, missF = time_perf2 - time_perf1)

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
    tmes = c(mode = time_perf2 - time_perf1, tmes)

    modeacc = mean(complete_CC[is.na(df[,1]),1, drop = F] == df_true[is.na(df[,1]),1, drop = F])
    modeval = mean(df[[1]] == levelone, na.rm = T)


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

    tmes = c(EFI = time_perf2 - time_perf1, tmes)

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
    c(result1, result2, tmes)

  }

if(MAR == T) mis_rate_vec = mis_rate_vec / 2

names(oper) <- mis_rate_vec

oper1 = lapply(oper, function(x) x[1:(nrow(x) / 3),])
oper2 = lapply(oper, function(x) x[(nrow(x) / 3 + 1):(nrow(x) / 3 * 2),])
oper3 = lapply(oper, function(x) x[(nrow(x) / 3 * 2 + 1):nrow(x),])

(res = do.call("rbind", lapply(oper1, function(x) sqrt(rowMeans((x - trueval)^2, na.rm = T)))))
res2 = do.call("rbind", lapply(oper2, rowMeans, na.rm = T))
res3 = do.call("rbind", lapply(oper3, rowMeans, na.rm = T))
res4 = do.call("rbind", lapply(oper2, apply, 1, sd, na.rm = T))

res_list = append(res_list, list(res))
res2_list = append(res2_list, list(res2))
res3_list = append(res3_list, list(res3))
res4_list = append(res4_list, list(res4))

(BIAS = do.call("rbind", lapply(oper1, function(x) rowMeans(x - trueval, na.rm = T)))) # BIAS
(SE = do.call("rbind", lapply(oper1, function(x) sqrt(apply(x, 1, function(y) mean((y - mean(y, na.rm = T))^2, na.rm = T) ))))) # SE
# apply(oper1[[1]], 1, function(x) mean(x) - trueval) # BIAS
# sqrt(apply(oper1[[1]], 1, function(x) mean((x - mean(x))^2))) # SE

##############
res_long = reshape2::melt(res)
res2_long = reshape2::melt(res2)
res3_long = reshape2::melt(res3)
names(res_long) = names(res2_long) = names(res3_long) = c("missrate", "method", "value")
# res_long <- res_long %>% filter(method != "pmm") # pmm not good

plot_res = ggplot(res_long, aes(x = missrate, y = value, group = method)) +
  geom_line(aes(color = method), linetype = 2) +
  geom_line(data = filter(res_long, method == "EFI"), aes(color = "EFI"), size = 2) +
  scale_x_continuous(breaks = mis_rate_vec, labels = (function(x) sprintf("%.2f", x)) ) +
  labs(x = NULL, y = NULL) +
  ggtitle(data_name) +
  # ggtitle(data_name, subtitle = paste("n = ", n, ", p = ", p)) +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = "top")

png(filename=paste(timenow0, "_", "RMSE", "_", data_name, ".png", sep = ""))
plot(plot_res)
dev.off()

plot_res2 = ggplot(res2_long, aes(x = missrate, y = value, group = method)) +
  geom_line(aes(color = method), linetype = 2) +
  geom_line(data = filter(res2_long, method == "EFI"), aes(color = "EFI"), size = 2) +
  scale_x_continuous(breaks = mis_rate_vec, labels = (function(x) sprintf("%.2f", x)) ) +
  labs(x = NULL, y = NULL) +
  ggtitle(data_name) +
  # ggtitle(data_name, subtitle = paste("n = ", n, ", p = ", p)) +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = "top")

png(filename=paste(timenow0, "_", "Accuracy", "_", data_name, ".png", sep = ""))
plot(plot_res2)
dev.off()

plot_res3 = ggplot(res3_long, aes(x = missrate, y = value, group = method)) +
  geom_line(aes(color = method), linetype = 2) +
  geom_line(data = filter(res3_long, method == "EFI"), aes(color = "EFI"), size = 2) +
  scale_x_continuous(breaks = mis_rate_vec, labels = (function(x) sprintf("%.2f", x)) ) +
  labs(x = NULL, y = NULL) +
  ggtitle(data_name) +
  # ggtitle(data_name, subtitle = paste("n = ", n, ", p = ", p)) +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"),
                 panel.background = element_rect(fill = "white"),
                 plot.title = element_text(size = 20, hjust = 0.5),
                 legend.title = element_text(size = 16),
                 legend.text = element_text(size = 14),
                 legend.position = "top")

png(filename=paste(timenow0, "_", "Time", "_", data_name, ".png", sep = ""))
plot(plot_res3)
dev.off()
#####################


save.image(paste(timenow0, data_name, ".RData", sep = ""))
}
stopCluster(myCluster)

timenow2 = Sys.time()
print(timenow2 - timenow1)

save.image(paste(timenow0, ".RData", sep = ""))

plot_list = NULL
plot_list2 = NULL
plot_list3 = NULL
for(idx in 1:length(res_list)){
  res = res_list[idx]
  res2 = res2_list[idx]
  res3 = res3_list[idx]
  data_name = data_list[idx]
  data_again = eval(parse(text = data_name))
  p = ncol(data_again)
  n = nrow(data_again)

  res_long = reshape2::melt(res)
  res2_long = reshape2::melt(res2)
  res3_long = reshape2::melt(res3)
  names(res_long)[1:3] = names(res2_long)[1:3] = names(res3_long)[1:3] = c("missrate", "method", "value")
  # res_long <- res_long %>% filter(method != "pmm") # pmm not good

  plot_res = ggplot(res_long, aes(x = missrate, y = log(value), group = method)) +
    geom_line(aes(color = method), linetype = 2) +
    geom_line(data = filter(res_long, method == "EFI"), aes(color = "EFI"), size = 2) +
    scale_x_continuous(breaks = mis_rate_vec, labels = (function(x) sprintf("%.2f", x)) ) +
    labs(x = NULL, y = NULL) +
    # ggtitle(data_name, subtitle = paste("n = ", n, ", p = ", p)) +
    ggtitle(data_name) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"),
          panel.background = element_rect(fill = "white"),
          plot.title = element_text(size = 20, hjust = 0.5),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.position = "top")

  plot_res2 = ggplot(res2_long, aes(x = missrate, y = log(value), group = method)) +
    geom_line(aes(color = method), linetype = 2) +
    geom_line(data = filter(res2_long, method == "EFI"), aes(color = "EFI"), size = 2) +
    scale_x_continuous(breaks = mis_rate_vec, labels = (function(x) sprintf("%.2f", x)) ) +
    labs(x = NULL, y = NULL) +
    # ggtitle(data_name, subtitle = paste("n = ", n, ", p = ", p)) +
    ggtitle(data_name) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"),
          panel.background = element_rect(fill = "white"),
          plot.title = element_text(size = 20, hjust = 0.5),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.position = "top")

  plot_res3 = ggplot(res3_long, aes(x = missrate, y = log(value), group = method)) +
    geom_line(aes(color = method), linetype = 2) +
    geom_line(data = filter(res3_long, method == "EFI"), aes(color = "EFI"), size = 2) +
    scale_x_continuous(breaks = mis_rate_vec, labels = (function(x) sprintf("%.2f", x)) ) +
    labs(x = NULL, y = NULL) +
    # ggtitle(data_name, subtitle = paste("n = ", n, ", p = ", p)) +
    ggtitle(data_name) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"),
                   panel.background = element_rect(fill = "white"),
                   plot.title = element_text(size = 20, hjust = 0.5),
                   legend.title = element_text(size = 16),
                   legend.text = element_text(size = 14),
                   legend.position = "top")

  plot_list = append(plot_list, list(plot_res))
  plot_list2 = append(plot_list2, list(plot_res2))
  plot_list3 = append(plot_list3, list(plot_res3))
}

plot_agg <- ggarrange(plotlist = plot_list, nrow = 1, common.legend = TRUE, legend="bottom")

png(filename=paste(ifelse(high_dim, "highd", "lowd"), "_", ifelse(MAR, "MAR", "MCAR"), "_","RMSE", ".png", sep = ""),  width = 300 * length(data_list),  height = 480)
annotate_figure(plot_agg, left = "log(RMSE)")
dev.off()

plot_agg2 <- ggarrange(plotlist = plot_list2, nrow = 1, common.legend = TRUE, legend="bottom")

png(filename=paste(ifelse(high_dim, "highd", "lowd"), "_", ifelse(MAR, "MAR", "MCAR"), "_", "Accuracy", ".png", sep = ""),  width = 300 * length(data_list),  height = 480)
annotate_figure(plot_agg2, left = "log(Accuracy)")
dev.off()

plot_agg3 <- ggarrange(plotlist = plot_list3, nrow = 1, common.legend = TRUE, legend="bottom")

png(filename=paste(ifelse(high_dim, "highd", "lowd"), "_", ifelse(MAR, "MAR", "MCAR"), "_", "Time", ".png", sep = ""),  width = 300 * length(data_list),  height = 480)
annotate_figure(plot_agg3, left = "log(Time)")
dev.off()

# load("~/GitHub/EFI/sim/ws1.RData")
