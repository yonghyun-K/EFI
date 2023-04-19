library(readr)
library(tidyr)
library(dplyr)
library(reshape2)
library(doParallel)

library(vcd)

library(mice)
library(missForest)
library(ggplot2)
library(Amelia)
library(EFI)

set.seed(123)

timenow1 = Sys.time()
timenow0 = gsub(' ', '_', gsub('[-:]', '', timenow1))
timenow = paste(timenow0, ".txt", sep = "")

# sink(timenow, append=TRUE)

spect_names <- c("OVERALL_DIAGNOSIS", "F1R", "F1S", "F2R", "F2S", "F3R", "F3S", "F4R", "F4S", "F5R", "F5S", "F6R", "F6S", "F7R", "F7S", "F8R", "F8S", "F9R", "F9S", "F10R", "F10S", "F11R", "F11S", "F12R", "F12S", "F13R", "F13S", "F14R", "F14S", "F15R", "F15S", "F16R", "F16S", "F17R", "F17S", "F18R", "F18S")
spect_data_train <- read_csv("https://archive.ics.uci.edu/ml/machine-learning-databases/spect/SPECT.train", col_names = spect_names,
                             col_types = cols(.default = "f"))
spect_data_test <- read_csv("https://archive.ics.uci.edu/ml/machine-learning-databases/spect/SPECT.test", col_names = spect_names,
                            col_types = cols(.default = "f"))
spect_data = rbind(spect_data_train, spect_data_test)

promoters_data <- read_csv("https://archive.ics.uci.edu/ml/machine-learning-databases/molecular-biology/promoter-gene-sequences/promoters.data", col_names=FALSE)
colnames(promoters_data) <- c("Class", "id", "Sequence")
promoters_data <- separate(promoters_data, Sequence, into = paste0("pos_", 1:57), sep = "")
promoters_data <- promoters_data %>% dplyr::select(-id, -pos_1) %>% mutate_all(factor)

mushroom_names <- c("class", "cap.shape", "cap.surface", "cap.color", "bruises", "odor",
                    "gill.attachment", "gill.spacing", "gill.size", "gill.color",
                    "stalk.shape", "stalk.root", "stalk.surface.above.ring",
                    "stalk.surface.below.ring", "stalk.color.above.ring",
                    "stalk.color.below.ring", "veil.type", "veil.color", "ring.number",
                    "ring.type", "spore.print.color", "population", "habitat")
mushroom_data <- read_csv("http://archive.ics.uci.edu/ml/machine-learning-databases/mushroom/agaricus-lepiota.data", col_names = mushroom_names,
                          col_types = cols(.default = "f"))
mushroom_data <- mushroom_data[,sapply(mushroom_data, nlevels) > 1]
# lapply(mushroom_data, levels)

column_names <- c("age_gt_60", "air", "airBoneGap", "ar_c", "ar_u", "bone", "boneAbnormal", "bser",
                  "history_buzzing", "history_dizziness", "history_fluctuating", "history_fullness",
                  "history_heredity", "history_nausea", "history_noise", "history_recruitment",
                  "history_ringing", "history_roaring", "history_vomiting", "late_wave_poor", "m_at_2k",
                  "m_cond_lt_1k", "m_gt_1k", "m_m_gt_2k", "m_m_sn", "m_m_sn_gt_1k", "m_m_sn_gt_2k",
                  "m_m_sn_gt_500", "m_p_sn_gt_2k", "m_s_gt_500", "m_s_sn", "m_s_sn_gt_1k", "m_s_sn_gt_2k",
                  "m_s_sn_gt_3k", "m_s_sn_gt_4k", "m_sn_2_3k", "m_sn_gt_1k", "m_sn_gt_2k", "m_sn_gt_3k",
                  "m_sn_gt_4k", "m_sn_gt_500", "m_sn_gt_6k", "m_sn_lt_1k", "m_sn_lt_2k", "m_sn_lt_3k",
                  "middle_wave_poor", "mod_gt_4k", "mod_mixed", "mod_s_mixed", "mod_s_sn_gt_500",
                  "mod_sn", "mod_sn_gt_1k", "mod_sn_gt_2k", "mod_sn_gt_3k", "mod_sn_gt_4k", "mod_sn_gt_500",
                  "notch_4k", "notch_at_4k", "o_ar_c", "o_ar_u", "s_sn_gt_1k", "s_sn_gt_2k", "s_sn_gt_4k",
                  "speech", "static_normal", "tymp", "viith_nerve_signs", "wave_V_delayed", "waveform_ItoV_prolonged",
                  "identifier", "class")
audiology <- read_csv("http://archive.ics.uci.edu/ml/machine-learning-databases/audiology/audiology.standardized.data",
                      col_types = cols(.default = "f"), col_names = column_names)
audiology <- (audiology %>% dplyr::select(-identifier))
audiology <- audiology[c(ncol(audiology), 1:(ncol(audiology)-1))]
# lapply(audiology, levels)


chess_data <- read_csv("http://archive.ics.uci.edu/ml/machine-learning-databases/chess/king-rook-vs-king-pawn/kr-vs-kp.data",
                       col_types = cols(.default = "f"))
colnames(chess_data) <- c("bkblk", "bknwy", "bkon8", "bkona", "bkspr", "bkxbq", "bkxcr", "bkxwp", "blxwp", "bxqsq",
         "cntxt", "dsopp", "dwipd", "hdchk", "katri", "mulch", "qxmsq", "r2ar8", "reskd", "reskr",
         "rimmx", "rkxwp", "rxmsq", "simpl", "skach", "skewr", "skrxp", "spcop", "stlmt", "thrsk",
         "wkcti", "wkna8", "wknck", "wkovl", "wkpos", "wtoeg", "won")
chess_data <- chess_data[c(ncol(chess_data), 1:(ncol(chess_data)-1))]


car_data <- read_csv("https://archive.ics.uci.edu/ml/machine-learning-databases/car/car.data",
                     col_names=c("buying", "maint", "doors", "persons", "lug_boot", "safety", "class_values"),
                     col_types = cols(.default = "f"))
colnames(car_data) <- c("Buying", "Maintenance", "Doors", "Persons", "Lug_Boot", "Safety", "Class_Values")
car_data <- car_data[c(ncol(car_data), 1:(ncol(car_data)-1))]

as.data.frame.table(HairEyeColor) %>% uncount(Freq)
as.data.frame.table(Titanic) %>% uncount(Freq)
as.data.frame.table(UCBAdmissions) %>% uncount(Freq)

esoph_data <- melt(esoph)
esoph_data <- esoph_data %>% uncount(value)
names(esoph_data)[ncol(esoph_data)] <- "Freq"

data(PreSex, package = "vcd")
as.data.frame.table(PreSex) %>% uncount(Freq)

# data(UKSoccer, package = "vcd")
# as.data.frame.table(UKSoccer)

df_true = spect_data
p = ncol(df_true)
n = nrow(df_true)
for(k in 1:p){
  df_true[[k]] <- as.factor(df_true[[k]])
}

modes = sapply(df_true, function(x) names(which.max(table(x))))

levelone = names(which.max(summary(df_true[[1]])))
trueval = mean(df_true[[1]] == levelone, na.rm = T)

mis_rate_vec = seq(from = 0.1, to = 0.9, length = 5)

cores = min(detectCores() - 3, 100)
myCluster <- makeCluster(cores, # number of cores to use
                         type = "PSOCK", setup_strategy = "sequential",
                         outfile = timenow) # type of cluster
registerDoParallel(myCluster)
SIMNUM = 200

oper <-
  foreach(mis_rate_idx= 1:length(mis_rate_vec), .packages = c("CVXR", "mice", "missForest", "Amelia", "EFI")) %:%
  foreach(simnum = 1:SIMNUM, .combine='cbind', .packages = c("CVXR", "mice", "missForest", "Amelia", "EFI")) %dopar% {
    print(paste("mis_rate_idx = ", mis_rate_idx))
    print(paste("simnum = ", simnum))
    mis_rate = mis_rate_vec[mis_rate_idx]
    obs_rate = 1 - mis_rate
    # print(mis_rate)
    df <- df_true

    # MCAR
    # for(k in 1:p){
    #   df[sample(1:n, round(n * mis_rate)),k] <- NA
    # }

    # MAR
    # p2 = round(p / 2)
    # delta = apply(df_true, 1, function(y){
    #   deltav = integer(p)
    #   tmpval = sample(1:p, 1)
    #   tmpidx = (tmpval:(tmpval + p2 - 1)) %% p + 1
    #   deltav[tmpidx] = 1
    #   if(tmpval %% 2 == 0){
    #     deltav[-tmpidx] = rbinom(p2,1,ifelse(y[tmpidx] == modes[tmpidx], min(c(3 * obs_rate, 1)), max(c(obs_rate / 3, 0))))[1:(p - p2)]
    #   }else{
    #     deltav[-tmpidx] = rbinom(p2,1,ifelse(y[tmpidx] == modes[tmpidx], max(c(obs_rate / 3, 0)), min(c(3 * obs_rate, 1))))[1:(p - p2)]
    #   }
    #   # deltav[-tmpidx] = rbinom(p2,1,ifelse(y[tmpidx] == 1, 0.1, 0.9))
    #
    #   deltav
    # })
    # delta = t(delta)

    tmpidx = order(abs(cor(sweep(df_true, 2, modes, "=="))[1,]), decreasing = T)[2]
    delta = apply(df_true, 1, function(y){
      deltav = rep(1, p)
      deltav[1] = rbinom(1,1,ifelse(y[tmpidx] == modes[tmpidx], min(c(3 * obs_rate, 1)), max(c(obs_rate / 3, 0))))
      deltav
    })
    delta = t(delta)
    mean(delta)

    df[delta == 0] <- NA
    print("Missing values are generated")

    # imp_mice <- mice(df, m = 5, seed = 123, printFlag = F)
    # df_mice = complete(imp_mice,1:5)

    # imp_missF <- missForest(cbind(1, df))
    # df_missF = imp_missF$ximp[-1]
    # print("missForest successfully done")

    # imp_Amelia <- amelia(df, noms = names(df), p2s = 0)
    # df_Amelia = do.call("rbind", imp_Amelia$imputations)
    # print("Amelia successfully done")

    edges_list = apply(rbind(1, 2:p), 2, list)
    # edges_list = apply(combn(p, 2), 2, list)
    dp = doublep(df, edges_list, R = 1)
    imp_EFI = efi(df, dp)
    print("EFI successfully done")

    result = c(CC = mean(df[[1]] == levelone, na.rm = T),
               EFI = sum((imp_EFI$imp[[1]] == levelone) * ( imp_EFI$imp$w)) / n
               # estimate(imp_EFI, "(Class_Values == \"unacc\")"),
               # MICE = mean(df_mice[[1]] == levelone),
               # missF = mean(df_missF[[1]] == levelone),
               # Amelia = mean(df_Amelia[[1]] == levelone, na.rm = T)
               )

    result

  }
stopCluster(myCluster)
names(oper) <- mis_rate_vec

res = do.call("rbind", lapply(oper, function(x) sqrt(rowMeans((x - trueval)^2))))

do.call("rbind", lapply(oper, function(x) rowMeans(x - trueval))) # BIAS
do.call("rbind", lapply(oper, function(x) sqrt(apply(x, 1, function(y) mean((y - mean(y))^2) )))) # SE
# apply(oper[[1]], 1, function(x) mean(x) - trueval) # BIAS
# sqrt(apply(oper[[1]], 1, function(x) mean((x - mean(x))^2))) # SE

res_long = reshape2::melt(res)
names(res_long) = c("missrate", "method", "value")

png(filename=paste(timenow0, "1.png", sep = ""))
ggplot(res_long, aes(x = missrate, y = value, color = method)) +
  geom_line() +
  scale_x_continuous(breaks = mis_rate_vec) +
  xlab("Missing rate") +
  ylab("RMSE") +
  # scale_color_manual(values = c("CC" = "red", "EFI" = "black", "MICE" = "blue", "missF" = "green", "Amelia" = "purple")) +
  # scale_color_manual(values = c("CC" = "#1b9e77", "EFI" = "#d95f02", "MICE" = "#7570b3", "missF" = "#e7298a", "Amelia" = "#66a61e")) +
  ggtitle("Comparison of imputation methods") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = "top",
        panel.grid.major = element_line(color = "gray", linetype = "dashed"))
dev.off()

timenow2 = Sys.time()
print(timenow2 - timenow1)




