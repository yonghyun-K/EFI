args <- commandArgs(trailingOnly = TRUE)

# args <- c(500, 300, 100, 50, 0.2)

# Simulation  setup ####
n = as.numeric(args[1])
p = 2
n_B = as.numeric(args[2])
B = as.numeric(args[3])
SIMNUM = as.numeric(args[4])
res1 = NULL
res2 = NULL
for(simnum in 1:SIMNUM){
  print(simnum)
  
  # Generate X and Y ####
  p_x = c(0.3, 0.25, 0.25, 0.2)
  X = matrix(rep(1:4, n * p)[rmultinom(n * p, 1, p_x) == 1], nr = n, nc = p)
  colnames(X) <- paste("X", 1:p, sep = "")
  # table(X) / n
  
  p_Y1 = c(0.4, 0.3, 0.4, 0.8)
  p1 = p_Y1[X[,1]]
  Y1 = rbinom(n, 1, p1)
  theta1 = sum(p_x * p_Y1)
  
  # p2 = c((0.8 + Y1) / 2, 0.8 * Y1 + 0.1, 0.9 - 0.5 * Y1, 0.6)[X[,1]]
  
  p2 = rowSums(cbind((0.8 + Y1) / 2, 0.8 * Y1 + 0.1, 0.9 - 0.5 * Y1, 0.6) * model.matrix(~ 0 + as.factor(X[,1]))) 
  Y2 = rbinom(n, 1, p2)
  
  p_Y2 = c((0.8 + 0.4) / 2, 0.8 * 0.3 + 0.1, 0.9 - 0.5 * 0.4, 0.6)
  theta2 = sum(p_x * p_Y2)
  
  theta = c(theta1, theta2)
  
  delta1 = rbinom(n, 1, 0.8)
  delta2 = rbinom(n, 1, 0.7)
  
  Y_ogn = cbind(Y1, Y2)
  theta_full = colMeans(Y_ogn)
  delta = cbind(delta1, delta2)
  
  Y = Y_ogn
  Y[delta == 0] = NA
  
  lambda = as.numeric(args[5])
  
  train_idx = sample(1:n, round(n * 0.9), replace = FALSE)
  x = X[train_idx,]
  y = Y[train_idx,]
  
  w_B = NULL
  bstp_idx_B = NULL
  
  w_01 = array(0, dim = c(4,4,2,2))
  w_10 = array(0, dim = c(4,4,2,2))
  w_00 = array(0, dim = c(4,4,2,2))
  
  for(b in 1:B){
    bstp_idx = sample(1:length(train_idx), n_B, replace = TRUE)
    x_b = x[bstp_idx,]
    y_b = y[bstp_idx,]
    
    x_oob = x[!(1:length(train_idx) %in% bstp_idx),]
    y_oob = y[!(1:length(train_idx) %in% bstp_idx),]
    
    # unique(cbind(x_b, y_b))
    # unique(cbind(x_b, y_b))[complete.cases(unique(cbind(x_b, y_b))),]
    
    # EM algorithm to compute \pi_{ijkl} ####
    # for i, j = 1,2,3,4 and k, l = 0, 1 
    
    n_mat = table(data.frame(cbind(x_b, y_b)), useNA = "ifany")
    # apply(n_mat, (p+1) : (p+2), sum)
    # apply(n_mat, 1:p, sum)
    
    # p_mat = n_mat[,,-3,-3] / sum(n_mat[,,-3,-3])
    
    # (n_mat[,,-3,-3] / sum(n_mat[,,-3,-3]))[,,2,2]
    # n_mat[,,3,3]
    
    # n_mat[,,3,-3, drop = F] # n^(1101)
    # p_mat[,,1,] * n_mat[,,3,-3] / apply(p_mat, c(1,2,4), sum)
    
    p_mat = array(1, dim = dim(n_mat[,,-3,-3])) / length(n_mat[,,-3,-3])
    
    while(T){
      nhat_11 = n_mat[,,-3,-3] # n^(1111)
      nhat_01 = sweep(p_mat, MARGIN = c(1,2,4), n_mat[,,3,-3] / apply(p_mat, c(1,2,4), sum), FUN = "*") #n_hat^(1101)
      nhat_10 = sweep(p_mat, MARGIN = c(1,2,3), n_mat[,,-3,3] / apply(p_mat, c(1,2,3), sum), FUN = "*") #n_hat^(1110)
      nhat_00 = sweep(p_mat, MARGIN = c(1,2), n_mat[,,3,3] / apply(p_mat, c(1,2), sum), FUN = "*") #n_hat^(1100)
      
      nhat_01[is.nan(nhat_01)] <- 0
      nhat_10[is.nan(nhat_10)] <- 0
      nhat_00[is.nan(nhat_00)] <- 0
      
      p_mat2 = (nhat_11 + nhat_01 + nhat_10 + nhat_00) / n_B
      
      # print(norm(p_mat2 - p_mat, "2"))
      if(norm(p_mat2 - p_mat, "2") < 10^(-3)) break
      else p_mat = p_mat2
    }
    
    # plot(c(table(data.frame(cbind(X[train_idx,][bstp_idx,], Y_ogn[train_idx,][bstp_idx,]))) / n_B),
    #      c(p_mat)); abline(0,1)
    
    # plot(c(n_mat[,,-3,-3] / sum(n_mat[,,-3,-3])), c(p_mat2)); abline(0,1)
    # plot(c(table(data.frame(cbind(X[train_idx,][bstp_idx,], Y_ogn[train_idx,][bstp_idx,]))) / n_B),
    #      c(n_mat[,,-3,-3] / sum(n_mat[,,-3,-3])))
    
    # unique(cbind(x_oob, y_oob))
    
    # Compute observed likelihood l_{obs}^{(b)} for each bag b ####
    
    z_oob = cbind(x_oob, y_oob)
    z_oob_11 = z_oob[!is.na(z_oob[,3]) & !is.na(z_oob[,4]), ]
    z_oob_01 = z_oob[is.na(z_oob[,3]) & !is.na(z_oob[,4]), ]
    z_oob_10 = z_oob[!is.na(z_oob[,3]) & is.na(z_oob[,4]), ]
    z_oob_00 = z_oob[is.na(z_oob[,3]) & is.na(z_oob[,4]), ]
    
    l_11 = apply(z_oob_11, 1, function(x) log(p_mat[t(as.character(x))]))
    l_01 = apply(z_oob_01, 1, function(x) {
      tmp = as.character(x)
      tmp1 = tmp; tmp1[3] <- "0"
      tmp2 = tmp; tmp2[3] <- "1"
      log(p_mat[t(tmp1)] + p_mat[t(tmp2)])
    })
    l_10 = apply(z_oob_10, 1, function(x) {
      tmp = as.character(x)
      tmp1 = tmp; tmp1[4] <- "0"
      tmp2 = tmp; tmp2[4] <- "1"
      log(p_mat[t(tmp1)] + p_mat[t(tmp2)])
    })
    l_00 = apply(z_oob_00, 1, function(x) {
      tmp = as.character(x)
      tmp1 = tmp; tmp1[3:4] <- c("0", "0")
      tmp2 = tmp; tmp2[3:4] <- c("0", "1")
      tmp3 = tmp; tmp3[3:4] <- c("1", "0")
      tmp4 = tmp; tmp4[3:4] <- c("1", "1")
      log(p_mat[t(tmp1)] + p_mat[t(tmp2)] + p_mat[t(tmp3)]+ p_mat[t(tmp4)])
    })
    
    l_obs = sum(c(l_11, l_01, l_10, l_00)) / nrow(z_oob) # observed likelihood
    l_obs 
    
    expl_obs = exp(lambda * l_obs)
    
    # Find the conditional probabilities ####
    # P(y_mis | y_obs, X)
    phat_01 = sweep(p_mat, MARGIN = c(1,2,4), apply(p_mat, c(1,2,4), sum), FUN = "/") 
    phat_10 = sweep(p_mat, MARGIN = c(1,2,3), apply(p_mat, c(1,2,3), sum), FUN = "/") 
    phat_00 = sweep(p_mat, MARGIN = c(1,2), apply(p_mat, c(1,2), sum), FUN = "/") 
    
    if(any(is.nan(phat_01))){
      # print("is.nan(phat_01)")
      # print(expl_obs)
      next
    }else if(any(is.nan(phat_10))){
      # print("is.nan(phat_10)")
      # print(expl_obs)
      next
    }else if(any(is.nan(phat_00))){
      # print("is.nan(phat_00)")
      # print(expl_obs)
      next
    }
    
    bstp_idx_B = cbind(bstp_idx_B, bstp_idx)
    
    w_B = c(w_B, expl_obs)
    
    w_01 = w_01 + expl_obs * phat_01
    w_10 = w_10 + expl_obs * phat_10
    w_00 = w_00 + expl_obs * phat_00
  }
  
  w_01 = w_01 / sum(w_B); w_10 = w_10 / sum(w_B); w_00 = w_00 / sum(w_B)
  
  # Fractional Imputation step ####
  z_imp = cbind(x,y)
  # table(data.frame(z))
  
  z_imp_11 = z_imp[!is.na(z_imp[,3]) & !is.na(z_imp[,4]), ]
  z_imp_10 = z_imp[!is.na(z_imp[,3]) & is.na(z_imp[,4]), ]
  z_imp_01 = z_imp[is.na(z_imp[,3]) & !is.na(z_imp[,4]), ]
  z_imp_00 = z_imp[is.na(z_imp[,3]) & is.na(z_imp[,4]), ]
  
  z_imp_11 = cbind(z_imp_11, 1)
  
  # z_imp_10 = z_imp[-(1:450),]
  
  if(nrow(z_imp_10)!= 0){
    fw_10 = apply(z_imp_10, 1, function(x) {
      tmp = as.character(x)
      tmp1 = tmp; tmp1[4] <- "0"
      tmp2 = tmp; tmp2[4] <- "1"
      c(w_10[t(tmp1)], w_10[t(tmp2)])
    })
    z_imp_10 = z_imp_10[rep(1:nrow(z_imp_10), each = 2),]
    z_imp_10 = cbind(z_imp_10, c(fw_10))
    z_imp_10[,"Y2"] <- rep(c(0,1), ncol(fw_10))
  }else{
    z_imp_10 <- cbind(z_imp_10, z_imp_10)[,-(6:8)]
  }
  
  if(nrow(z_imp_01)!= 0){
    fw_01 = apply(z_imp_01, 1, function(x) {
      tmp = as.character(x)
      tmp1 = tmp; tmp1[3] <- "0"
      tmp2 = tmp; tmp2[3] <- "1"
      c(w_01[t(tmp1)], w_01[t(tmp2)])
    })
    z_imp_01 = z_imp_01[rep(1:nrow(z_imp_01), each = 2),]
    z_imp_01 = cbind(z_imp_01, c(fw_01))
    z_imp_01[,"Y1"] <- rep(c(0,1), ncol(fw_01))
  }else{
    z_imp_01 <- cbind(z_imp_01, z_imp_01)[,-(6:8)]
  }

  if(nrow(z_imp_00)!= 0){
    fw_00 = apply(z_imp_00, 1, function(x) {
      tmp = as.character(x)
      tmp1 = tmp; tmp1[3:4] <- c("0", "0")
      tmp2 = tmp; tmp2[3:4] <- c("0", "1")
      tmp3 = tmp; tmp3[3:4] <- c("1", "0")
      tmp4 = tmp; tmp4[3:4] <- c("1", "1")
      c(w_00[t(tmp1)], w_00[t(tmp2)], w_00[t(tmp3)], w_00[t(tmp4)])
    })
    z_imp_00 = z_imp_00[rep(1:nrow(z_imp_00), each = 4),]
    z_imp_00 = cbind(z_imp_00, c(fw_00))
    z_imp_00[,c("Y1", "Y2")] <- cbind(rep(c(0,0,1,1), ncol(fw_00)), rep(c(0,1,0,1), ncol(fw_00)))
  }else{
    z_imp_00 <- cbind(z_imp_00, z_imp_00)[,-(6:8)]
  }

  # Final estimates ####
  # sum(c(z_imp_10[,5]))
  theta_prop1 = sum(c(z_imp_11[,"Y1"] * z_imp_11[,5], z_imp_10[,"Y1"] * z_imp_10[,5], z_imp_01[,"Y1"] * z_imp_01[,5], z_imp_00[,"Y1"] * z_imp_00[,5])) /
    sum(c(z_imp_11[,5], z_imp_10[,5], z_imp_01[,5], z_imp_00[,5]))
  
  theta_prop2 = sum(c(z_imp_11[,"Y2"] * z_imp_11[,5], z_imp_10[,"Y2"] * z_imp_10[,5], z_imp_01[,"Y2"] * z_imp_01[,5], z_imp_00[,"Y2"] * z_imp_00[,5])) /
    sum(c(z_imp_11[,5], z_imp_10[,5], z_imp_01[,5], z_imp_00[,5]))
  # c(theta1, theta2)
  theta_full
  
  theta_prop = c(theta_prop1, theta_prop2)
  
  theta_cc = colMeans(y, na.rm = T)
  
  res1 = rbind(res1, c(theta_full[1], theta_cc[1], theta_prop[1]))
  res2 = rbind(res2, c(theta_full[2], theta_cc[2], theta_prop[2]))
}

# Summary table ####
BIAS = colMeans(res1 - theta1)
SE = apply(res1, 2, function(x) sqrt(var(x) * (length(x)-1)/length(x) ))
RMSE = apply(res1 - theta1, 2, function(x) sqrt(mean(x^2)))

tmp_tbl <- round(cbind(SB = BIAS / RMSE, SE = SE, RMSE = RMSE), 3)
tmp_tbl
xtable::xtable(tmp_tbl, digits = 3)

BIAS = colMeans(res2 - theta2)
SE = apply(res2, 2, function(x) sqrt(var(x) * (length(x)-1)/length(x) ))
RMSE = apply(res2 - theta2, 2, function(x) sqrt(mean(x^2)))

tmp_tbl <- round(cbind(SB = BIAS / RMSE, SE = SE, RMSE = RMSE), 3)
tmp_tbl
xtable::xtable(tmp_tbl, digits = 3)
