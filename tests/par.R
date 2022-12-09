# Simulation  setup ####
n = 2000
print(paste("n(Sample size) =", n))

n_B = 1500
print(paste("n_B(Bootstrap sample size) =", n_B))

p = 100
print(paste("p =", p))

p_star = 2
print(paste("p_star =", p_star))

q = 2
print(paste("q =", q))

B = 30
print(paste("B(The number of bootstraps) =", B))

SIMNUM = 100
print(paste("SIMNUM(MC size) =", SIMNUM))

k_x = 2
print(paste("k_x(Number of categories of X) =", k_x))

lambda = 120
print(paste("lambda =", lambda))

cores=SIMNUM
print(paste("cores =", cores))
cl <- makeSOCKcluster(cores)
registerDoSNOW(cl)

progress <- function(n) cat(sprintf("task %d is complete\n", n))
opts <- list(progress=progress)

set.seed(1)

p_x = runif(k_x, min = 0.25, max = 0.75)
# p_x = rep(1, k_x)
p_x = p_x / sum(p_x)
print(paste("p_x =", paste(round(p_x,5), collapse = " ")))

p_Y = runif(20, min = 0.25, max = 0.75)
print(paste("p_Y[1:k_x] =", paste(round(p_Y[1:k_x],5), collapse = " ")))

theta = sum(p_x * p_Y[1:k_x])
# theta = 1/2
print(paste("theta =", round(theta, 5)))

# p_delta_ftn = function(k) function(k) rep(0.5, n); print("MCAR") # 0.5

#  p_delta_ftn = function(k) X_num[,k] / 5 + 0.2; print("MAR")# 0.5
#  p_delta_ftn = function(k) X_num[,k] / 3 + 1/3; print("MAR") # 0.84

# p_delta_ftn = function(k) 1 / (1 + exp(-(X_num[,k] + X_num[,k+1] - X_num[,k+2]))); print("MAR") # 0.8
# p_delta_ftn = function(k) 1 / (1 + exp(-(-1 + X_num[,k] + X_num[,k+1] - X_num[,k+2]))); print("MAR") # 0.61
p_delta_ftn = function(k) 1 / (1 + exp(-(-X_num[,k] / 4 + X_num[,k+2] / 4 - X_num[,k+3] / 4))); print("MAR") # 0.415
# p_delta_ftn = function(k) 1 / (1 + exp(-(-X_num[,k] / 4 + 1 / 4))); print("MAR") # 0.415

# p_delta_ftn = function(k) 1 / (1 + exp(-(1 - X_num[,k] / 4 + 1 / 4 - Y_num[,k] / 4))); print("NMAR") # 0.415
# p_delta_ftn = function(k) 1 / (1 + exp( -(1.5 + X_num[,k] + X_num[,k+1] - X_num[,k+2] - 5 * Y_num[,k]))); print("NMAR") # 0.561
# p_delta_ftn = function(k) 1 / (1 + exp( -(X_num[,k] + X_num[,k+1] - 5 * Y_num[,k]))); print("NMAR")
# p_delta_ftn = function(k) 1 / (1 + exp( -( X_num[,k] / 4 - Y_num[,k] / 4))); print("NMAR")
