# args <- commandArgs(trailingOnly = TRUE)

# args <- c(n, n_B, B, SIMNUM, lambda)
args <- c(2000, 1500, 45, 100)

# Simulation  setup ####
n = as.numeric(args[1])
p = 10
q = 2
n_B = as.numeric(args[2])
B = as.numeric(args[3])
SIMNUM = as.numeric(args[4])
res = NULL
K = 4
set.seed(1)

lambda_vec = seq(from = 10, to = 100, len = 10)



# Generate X and Y ####
p_x = c(0.3, 0.25, 0.25, 0.2)
X = matrix(rep(1:4, n * p)[rmultinom(n * p, 1, p_x) == 1], nr = n, nc = p)
colnames(X) <- paste("X", 1:p, sep = "")
# table(X) / n

# p_Y1 = c(0.4, 0.3, 0.4, 0.8)
p_Y1 = c(0.2, 0.3, 0.6, 0.8)
p1 = p_Y1[X[, 1]]
Y1 = rbinom(n, 1, p1)
theta1 = sum(p_x * p_Y1)

# p2 = c((0.8 + Y1) / 2, 0.8 * Y1 + 0.1, 0.9 - 0.5 * Y1, 0.6)[X[,1]]

p2 = rowSums(cbind((0.8 + Y1) / 2, 0.8 * Y1 + 0.1, 0.9 - 0.5 * Y1, 0.6) * model.matrix( ~ 0 + as.factor(X[, 1])))
Y2 = rbinom(n, 1, p2)

cor(X[, 1], Y2) / cor(X[, 1], Y1)

p_Y2 = c((0.8 + p_Y1[1]) / 2, 0.8 * p_Y1[2] + 0.1, 0.9 - 0.5 * p_Y1[3], 0.6)
theta2 = sum(p_x * p_Y2)

theta = c(theta1, theta2)

# p_delta1 = X[,1] / 5 + 1 / 5
p_delta1 = X[, 1] / 8 + 0.5
# p_delta1 = X[,1] / 25 + 0.7
p_delta2 = X[, 1] / 6 + 1 / 3

# p_delta1 = 1 / (1 + exp(-(2 + x[,1] - x[,5]))) # 0.8
# p_delta2 = 1 / (1 + exp(-(-1 + x[,1] + x[,2] - x[,3]))) # 0.7

# p_delta1 = 0.8
# p_delta2 = 0.7

delta1 = rbinom(n, 1, p_delta1)
delta2 = rbinom(n, 1, p_delta2)

Y_ogn = cbind(Y1, Y2)
theta_full = colMeans(Y_ogn)
delta = cbind(delta1, delta2)

Y = Y_ogn
Y[delta == 0] = NA

bstp_idx_B = NULL
for (b in 1:B) {
  bstp_idx = sample(1:((1 - 1 / K) * n), n_B, replace = TRUE)
  bstp_idx_B = cbind(bstp_idx_B, bstp_idx)
}

for (lambda in lambda_vec) {
  print(lambda)
  res_K = rep(0, K)
  
  cuts = cut(sample(1:n, n), breaks = 10, labels = FALSE)
  for (k in 1:K) {
    # flds = createFolds(1:n, k = K, list = TRUE)
    # test_idx = flds[[k]]
    # train_idx = do.call(c, flds[-k])
    
    
    test_idx = which(cuts == k)
    train_idx = which(cuts != k)
    
    names(train_idx) = NULL
    
    # train_idx = sample(1:n, round(n * 0.9), replace = FALSE)
    # train_idx = 1:n
    x = X[train_idx, ]
    y = Y[train_idx, ]
    
    x_test = X[test_idx,]
    y_test = Y[test_idx,]
    
    w_B = NULL
    # bstp_idx_B = NULL
    
    p_mat_lambda = array(0, dim = c(4, 4, 2, 2))
    
    for (b in 1:B) {
      # select_x = sample(1:p, q, replace = F)
      select_x = combn(p, q)[, (b + 44) %% 45 + 1]
      # select_x = c(1, 2)
      # select_x = c(3, 4)
      
      bstp_idx = bstp_idx_B[,b]
      x_b = x[bstp_idx, select_x]
      y_b = y[bstp_idx, ]
      
      x_oob = x[!(1:length(train_idx) %in% bstp_idx), select_x]
      y_oob = y[!(1:length(train_idx) %in% bstp_idx), ]
      
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
      
      p_01 = array(1, dim = c(4, 4, 2, 2)) / 2 # P(Y1 | x1, x2, y2, delta1 = 0, delta2 = 1)
      p_10 = array(1, dim = c(4, 4, 2, 2)) / 2 # P(Y2 | x1, x2, y1, delta1 = 1, delta2 = 0)
      p_00 = array(1, dim = c(4, 4, 2, 2)) / 4 # P(Y1, Y2 | x1, x2, delta1 = 0, delta2 = 0)
      
      while (T) {
        n11_hat = n_mat[, , -3, -3]
        n01_hat = sweep(p_01, MARGIN = c(1, 2, 4), n_mat[, , 3, -3], "*") # sum P(Y1 | x, Y2, delta1 = 0, delta2 = 1)
        n10_hat = sweep(p_10, MARGIN = c(1, 2, 3), n_mat[, , -3, 3], "*") # sum P(Y2 | x, Y1, delta1 = 1, delta2 = 0)
        n00_hat = sweep(p_00, MARGIN = c(1, 2), n_mat[, , 3, 3], "*") # sum P(Y1, Y2 | x, delta1 = delta2 = 0)
        
        n_hat = n11_hat + n01_hat + n10_hat + n00_hat # \hat N(x1, x2, y1, y2)
        
        Py_x = sweep(n_hat, MARGIN = c(1, 2), apply(n_hat, c(1, 2), sum), "/") # P(Y1, Y2 | x)
        
        n1p_hat = n11_hat + n10_hat
        n0p_hat = n01_hat + n00_hat
        Pdel1_xy = apply(n1p_hat, c(1, 2, 3), sum) / apply(n_hat, c(1, 2, 3), sum) # P(delta1 = 1 | x, y1)
        
        np1_hat = n11_hat + n01_hat
        np0_hat = n10_hat + n00_hat
        Pdel2_xy = apply(np1_hat, c(1, 2, 4), sum) / apply(n_hat, c(1, 2, 4), sum) # P(delta2 = 1 | x, y2)
        
        # P(y1 | x, y2, delta1 = 0, delta2 = 1) \propto
        # P(del1 = 0 | x, y1) * P(del2 = 1 | x, y2) * p(y1, y2 | x)
        p_01_new = sweep(Py_x, MARGIN = c(1, 2, 3), (1 - Pdel1_xy), "*")
        p_01_new = sweep(p_01_new, MARGIN = c(1, 2, 4), Pdel2_xy, "*")
        p_01_new = sweep(p_01_new, MARGIN = c(1, 2, 4), apply(p_01_new, c(1, 2, 4), sum), "/")
        
        # P(y2 | x, y1, delta1 = 1, delta2 = 0) \propto
        # P(del1 = 1 | x, y1) * P(del2 = 0 | x, y2) * p(y1, y2 | x)
        p_10_new = sweep(Py_x, MARGIN = c(1, 2, 3), Pdel1_xy, "*")
        p_10_new = sweep(p_10_new, MARGIN = c(1, 2, 4), (1 - Pdel2_xy), "*")
        p_10_new = sweep(p_10_new, MARGIN = c(1, 2, 3), apply(p_10_new, c(1, 2, 3), sum), "/")
        
        # P(y1, y2 | x, delta1 = 0, delta2 = 0) \propto
        # P(del1 = 0 | x, y1) * P(del2 = 0 | x, y2) * p(y1, y2 | x)
        p_00_new = sweep(Py_x, MARGIN = c(1, 2, 3), (1 - Pdel1_xy), "*")
        p_00_new = sweep(p_00_new, MARGIN = c(1, 2, 4), (1 - Pdel2_xy), "*")
        p_00_new = sweep(p_00_new, MARGIN = c(1, 2), apply(p_00_new, c(1, 2), sum), "/")
        
        # If # of (x1, x2 & (y1 = NA | y2 = NA)) = 0, we ignore the associated probs.
        p_01_new[is.na(p_01_new)] = 0
        p_10_new[is.na(p_10_new)] = 0
        p_00_new[is.na(p_00_new)] = 0
        
        diff = norm(p_01_new - p_01, "2") + norm(p_10_new - p_10, "2") +
          norm(p_00_new - p_00, "2")
        # print(diff)
        # if(diff < 10^(-3)){
        if (diff < 10 ^ (-1)) {
          break
        }
        else{
          p_01 = p_01_new
          p_10 = p_10_new
          p_00 = p_00_new
        }
      }
      p_01 = p_01_new
      p_10 = p_10_new
      p_00 = p_00_new
      
      p_mat = n_hat / nrow(x_b)
      
      # Compute observed likelihood l_{obs}^{(b)} for each bag b ####
      
      z_oob = cbind(x_oob, y_oob)
      z_oob_11 = z_oob[!is.na(z_oob[, 3]) & !is.na(z_oob[, 4]),]
      z_oob_01 = z_oob[is.na(z_oob[, 3]) & !is.na(z_oob[, 4]),]
      z_oob_10 = z_oob[!is.na(z_oob[, 3]) & is.na(z_oob[, 4]),]
      z_oob_00 = z_oob[is.na(z_oob[, 3]) & is.na(z_oob[, 4]),]
      
      l_11 = apply(z_oob_11, 1, function(x)
        log(p_mat[t(as.character(x))] * Pdel1_xy[t(as.character(x[-4]))] * Pdel2_xy[t(as.character(x[-3]))]))
      l_01 = apply(z_oob_01, 1, function(x) {
        tmp = as.character(x)
        tmp1 = tmp
        tmp1[3] <- "0"
        tmp2 = tmp
        tmp2[3] <- "1"
        log(p_mat[t(tmp1)] * (1 - Pdel1_xy[t(tmp1[-4])]) * Pdel2_xy[t(tmp1[-3])] +
              p_mat[t(tmp2)] * (1 - Pdel1_xy[t(tmp2[-4])]) * Pdel2_xy[t(tmp2[-3])])
      })
      l_10 = apply(z_oob_10, 1, function(x) {
        tmp = as.character(x)
        tmp1 = tmp
        tmp1[4] <- "0"
        tmp2 = tmp
        tmp2[4] <- "1"
        log(p_mat[t(tmp1)] * Pdel1_xy[t(tmp1[-4])] * (1 - Pdel2_xy[t(tmp1[-3])]) +
              p_mat[t(tmp2)] * Pdel1_xy[t(tmp2[-4])] * (1 - Pdel2_xy[t(tmp2[-3])]))
      })
      l_00 = apply(z_oob_00, 1, function(x) {
        tmp = as.character(x)
        tmp1 = tmp
        tmp1[3:4] <- c("0", "0")
        tmp2 = tmp
        tmp2[3:4] <- c("0", "1")
        tmp3 = tmp
        tmp3[3:4] <- c("1", "0")
        tmp4 = tmp
        tmp4[3:4] <- c("1", "1")
        log(
          p_mat[t(tmp1)] * (1 - Pdel1_xy[t(tmp1[-4])]) * (1 - Pdel2_xy[t(tmp1[-3])]) +
            p_mat[t(tmp2)] * (1 - Pdel1_xy[t(tmp2[-4])]) * (1 - Pdel2_xy[t(tmp2[-3])]) +
            p_mat[t(tmp3)] * (1 - Pdel1_xy[t(tmp3[-4])]) * (1 - Pdel2_xy[t(tmp3[-3])]) +
            p_mat[t(tmp4)] * (1 - Pdel1_xy[t(tmp4[-4])]) * (1 - Pdel2_xy[t(tmp4[-3])])
        )
        
      })
      
      l_obs = sum(c(l_11, l_01, l_10, l_00)) / nrow(z_oob) # observed likelihood
      l_obs
  
      expl_obs = exp(lambda * l_obs)
      w_B = c(w_B, expl_obs)
      
      p_mat_lambda = p_mat_lambda + expl_obs * p_mat
      
    }
    
    if(sum(w_B) == 0){
      stop(sum(w_B) == 0)
    }
    
    p_mat_lambda = p_mat_lambda / sum(w_B)
    
    py = apply(p_mat_lambda, c(3, 4), sum) # P(y = a)
    
    tmp_ftn = function(y, py) {
      if (!is.na(y[1]) & !is.na(y[2])) {
        res = log(py[y[1] + 1, y[2] + 1])
      } else if (is.na(y[1]) & !is.na(y[2])) {
        res = log(sum(py[, y[2] + 1]))
      } else if (!is.na(y[1]) & is.na(y[2])) {
        res = log(sum(py[y[1] + 1,]))
      } else if (is.na(y[1]) & is.na(y[2])) {
        res = 0
      }
      return(res)
    }
    
    l_lambda = sum(apply(y_test, 1, tmp_ftn, py = py))
    
    res_K[k] = l_lambda
    
    
  }
  res = c(res, mean(res_K))
}
# Summary table ####
plot(lambda_vec, res)