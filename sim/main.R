library(readr)
library(tidyr)
library(dplyr)
library(doParallel)

library(mice)
library(missForest)
library(ggplot2)
library(Amelia)
library(EFI)

timenow1 = Sys.time()
timenow0 = gsub(' ', '_', gsub('[-:]', '', timenow1))
timenow = paste(timenow0, ".txt", sep = "")

# sink(timenow, append=TRUE)

spect_names <- c("OVERALL_DIAGNOSIS", "F1R", "F1S", "F2R", "F2S", "F3R", "F3S", "F4R", "F4S", "F5R", "F5S", "F6R", "F6S", "F7R", "F7S", "F8R", "F8S", "F9R", "F9S", "F10R", "F10S", "F11R", "F11S", "F12R", "F12S", "F13R", "F13S", "F14R", "F14S", "F15R", "F15S", "F16R", "F16S", "F17R", "F17S", "F18R", "F18S")
spect_data_train <- read_csv("https://archive.ics.uci.edu/ml/machine-learning-databases/spect/SPECT.train", col_names = spect_names)
spect_data_test <- read_csv("https://archive.ics.uci.edu/ml/machine-learning-databases/spect/SPECT.test", col_names = spect_names)
spect_data = rbind(spect_data_train, spect_data_test)

promoters_data <- read_csv("https://archive.ics.uci.edu/ml/machine-learning-databases/molecular-biology/promoter-gene-sequences/promoters.data", col_names=FALSE)
# Define the column names
colnames(promoters_data) <- c("Class", "id", "Sequence")
promoters_data <- separate(promoters_data, Sequence, into = paste0("pos_", 1:57), sep = "")
promoters_data <- select(promoters_data, -id, -pos_1)

car_data <- read_csv("https://archive.ics.uci.edu/ml/machine-learning-databases/car/car.data", col_names=c("buying", "maint", "doors", "persons", "lug_boot", "safety", "class_values"))
colnames(car_data) <- c("Buying", "Maintenance", "Doors", "Persons", "Lug_Boot", "Safety", "Class_Values")
car_data <- car_data[c(ncol(car_data), 1:(ncol(car_data)-1))]

df_true = spect_data
p = ncol(df_true)
n = nrow(df_true)
for(k in 1:p){
  df_true[[k]] <- as.factor(df_true[[k]])
}

levelone = names(which.max(summary(df_true[[1]])))
trueval = mean(df_true[[1]] == levelone, na.rm = T)

miss_rate_vec = seq(from = 0.1, to = 0.5, length = 4)

cores = min(detectCores() - 3, 100)
myCluster <- makeCluster(cores, # number of cores to use
                         type = "PSOCK",
                         outfile = timenow) # type of cluster
registerDoParallel(myCluster)
SIMNUM = 5

oper <-
  foreach(miss_rate_idx= 1:length(miss_rate_vec), .packages = c("CVXR", "mice", "missForest", "Amelia", "EFI")) %:%
  foreach(simnum = 1:SIMNUM, .combine='cbind', .packages = c("CVXR", "mice", "missForest", "Amelia", "EFI")) %dopar% {
    print(paste("miss_rate_idx = ", miss_rate_idx))
    print(paste("simnum = ", simnum))
    miss_rate = miss_rate_vec[miss_rate_idx]
    # print(miss_rate)
    df <- df_true

    for(k in 1:p){
      df[sample(1:n, round(n * miss_rate)),k] <- NA
    }

    # imp_mice <- mice(df, m = 5, seed = 123, printFlag = F)
    # df_mice = complete(imp_mice,1:5)

    imp_missF <- missForest(cbind(1, df))
    df_missF = imp_missF$ximp[-1]

    imp_Amelia <- amelia(df, noms = names(df), p2s = 0)
    df_Amelia = do.call("rbind", imp_Amelia$imputations)

    edges_list = apply(rbind(1, 2:p), 2, list)
    # edges_list = apply(combn(p, 2), 2, list)
    dp = doublep(df, edges_list)
    imp_EFI = efi(df, dp)

    result = c(CC = mean(df[[1]] == levelone, na.rm = T),
               EFI = sum((imp_EFI$imp[[1]] == levelone) * ( imp_EFI$imp$w)) / n,
               # estimate(imp_EFI, "(Class_Values == \"unacc\")"),
               # MICE = mean(df_mice[[1]] == levelone),
               missF = mean(df_missF[[1]] == levelone),
               Amelia = mean(df_Amelia[[1]] == levelone, na.rm = T))

    result

  }
stopCluster(myCluster)
names(oper) <- miss_rate_vec

res = do.call("rbind", lapply(oper, function(x) sqrt(rowMeans((x - trueval)^2))))

do.call("rbind", lapply(oper, function(x) rowMeans(x - trueval)))

res_long = reshape2::melt(res)
names(res_long) = c("missrate", "method", "value")

ggplot(res_long, aes(x = missrate, y = value, color = method)) +
  geom_line() +
  scale_x_continuous(breaks = miss_rate_vec) +
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

timenow2 = Sys.time()
print(timenow2 - timenow1)

ggplot(res_long, aes(x = missrate, y = value, color = method)) +
  geom_line() +
  scale_x_continuous(breaks = miss_rate_vec) +
  xlab("Missing rate") +
  ylab("RMSE") +
  scale_color_manual(values = c("CC" = "red", "EFI" = "black", "MICE" = "blue", "missF" = "green", "Amelia" = "purple")) +
  ggtitle("Comparison of imputation methods") +
  theme_bw()


