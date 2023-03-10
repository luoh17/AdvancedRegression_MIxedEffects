library(geepack) 
library(nlme) 
library(MASS) 
library(broom.mixed)
# Set the regression coefficients
beta_0 <- 2
beta_1 <- 0.5 
beta_2 <- 1
m<-300 
n<-3
# Set the correlation parameter
rho <- 0.75
# Set the covariance matrix for the underlying normal

sigma <- matrix(c(1, rho, rho, rho, 1, rho, rho, rho, 1), nrow = n, ncol = n)
data_generation <- function() {
  time <- matrix(rep(1:n, times = m), nrow = m, ncol = n, byrow = TRUE)
  effect <- rnorm(100, 0, 1)
  mu <- beta_0 + beta_1 * effect + beta_2 * time 
  data <- matrix(rep(0, m*n), nrow = m, ncol = n)
  for (i in 1:m) {
    data[i, ] <- mvrnorm(n = 1, mu = mu[i, ], Sigma = sigma)
  }
  df = data.frame(id = rep(1:m, each = 3),
                  effects = rep(effect, each = 3),
                  times = rep(1:n, times = m),
                  y = matrix(t(data), nrow = m*n, ncol = 1))
  return(df) 
  }


run_one_times <- function(missfrac = 0.2) { # data generation
  df_sim <- data_generation()
  result <- matrix(rep(0, 36), nrow = 6, ncol = 6)
  #MCAR for gee
  missing_rows_MCAR <- sample(nrow(df_sim), size=round( missfrac*nrow(df_sim)), replace=FALSE) 
  df_MCAR <- df_sim
  df_MCAR[missing_rows_MCAR, "y"] <- NA
  gee_mod_MCAR <- geeglm(y ~ effects + times, data = df_MCAR, id = id, corstr = "exchange", family = gaussian(),
                         na.action = na.omit) 
  temp <- tidy(gee_mod_MCAR)
  result[1, 1:3] <- temp$estimate
  result[1, 4:6] <- temp$std.error
  #MAR for gee
  prob_MAR <- sapply(df_sim$effects, function(x) ifelse(x >= 0, 0.2, 0.5))
  missing_rows_MAR <- sample(nrow(df_sim), prob = prob_MAR, size = round(missfrac*nrow(df_sim))
                             , replace = FALSE)
  df_MAR <- df_sim
  df_MAR[missing_rows_MAR, "y"] <- NA
  gee_mod_MAR <- geeglm(y ~ effects + times, data = df_MAR, id = id,
                        corstr = "exchange", family = gaussian(),
                        na.action = na.omit) 
  temp <- tidy(gee_mod_MAR)
  result[2, 1:3] <- temp$estimate
    result[2, 4:6] <-temp$std.error
    #MNAR for gee
  missing_rows_MNAR <- sample(nrow(df_sim), prob = 0.05 * df_sim$times,
                              size = round(missfrac*nrow(df_sim)), replace = FALSE)
  df_MNAR <- df_sim
  df_MNAR[missing_rows_MNAR, "y"] <- NA
  gee_mod_MNAR <- geeglm(y ~ effects + times, data = df_MNAR,
                         id = id,
                         corstr = "exchange", family =gaussian(),
                         na.action = na.omit) 
  temp <- tidy(gee_mod_MNAR)
  result[3, 1:3] <- temp$estimate 
  result[3, 4:6] <- temp$std.error
  # lmm for MCAR
  lmm_mod_MCAR <- lme(y ~ effects + times, data = df_MCAR, random = ~1|id, na.action = na.omit)
  temp <- tidy(lmm_mod_MCAR)
  result[4, 1:3] <- temp$estimate[1:3] 
  result[4, 4:6] <- temp$std.error[1:3]
  # lmm for MAR
  lmm_mod_MAR <- lme(y ~ effects + times, data = df_MAR, random = ~1|id, na.action = na.omit)
  temp <- tidy(lmm_mod_MAR)
  result[5, 1:3] <- temp$estimate[1:3] 
  result[5, 4:6] <- temp$std.error[1:3]
  #lmmforMNAR
  lmm_mod_MNAR <- lme(y ~ effects + times, data = df_MNAR, random = ~1|id, na.action = na.omit)
  temp <- tidy(lmm_mod_MNAR)
  result[6, 1:3] <- temp$estimate[1:3] 
  result[6, 4:6] <- temp$std.error[1:3]
  result_df <- data.frame(result)
  colnames(result_df) <- c("beta_0 estimate", "beta_1
estimate", "beta_2 estimate", "beta_0 S.E.", "beta_1 S.E.", "beta_2 S.E.")
  rownames(result_df) <- c("gee_MCAR", "gee_MAR", "gee_MNAR", "lmm_MCAR", "lmm_MAR", "lmm_MNAR")
  return(result_df) 
  }


rep_time <- 300
beta_0_gee_MCAR <- c()
beta_1_gee_MCAR <- c()
beta_2_gee_MCAR <- c()
beta_0_gee_MAR <- c()
beta_1_gee_MAR <- c()
beta_2_gee_MAR <- c()
beta_0_gee_MNAR <- c()
beta_1_gee_MNAR <- c()
beta_2_gee_MNAR <- c()

beta_0_lmm_MCAR <- c()
beta_1_lmm_MCAR <- c()
beta_2_lmm_MCAR <- c()
beta_0_lmm_MAR <- c()
beta_1_lmm_MAR <- c()
beta_2_lmm_MAR <- c()
beta_0_lmm_MNAR <- c()
beta_1_lmm_MNAR <- c()
beta_2_lmm_MNAR <- c()

beta_0_gee_MCAR_se <- c()
beta_1_gee_MCAR_se <- c()
beta_2_gee_MCAR_se <- c()
beta_0_gee_MAR_se <- c()
beta_1_gee_MAR_se <- c()
beta_2_gee_MAR_se <- c()
beta_0_gee_MNAR_se <- c()
beta_1_gee_MNAR_se <- c()
beta_2_gee_MNAR_se <- c()

beta_0_lmm_MCAR_se <- c()
beta_1_lmm_MCAR_se <- c()
beta_2_lmm_MCAR_se <- c()
beta_0_lmm_MAR_se <- c()
beta_1_lmm_MAR_se <- c()
beta_2_lmm_MAR_se <- c()
beta_0_lmm_MNAR_se <- c()
beta_1_lmm_MNAR_se <- c()
beta_2_lmm_MNAR_se <- c()


set.seed(123)
for (i in 1:rep_time) {
  r <- run_one_times()
  beta_0_gee_MCAR <- c(beta_0_gee_MCAR, r[1, 1])
  beta_1_gee_MCAR <- c(beta_1_gee_MCAR, r[1, 2]) 
  beta_2_gee_MCAR <- c(beta_2_gee_MCAR, r[1, 3]) 
  beta_0_gee_MCAR_se <- c(beta_0_gee_MCAR_se, r[1, 4]) 
  beta_1_gee_MCAR_se <- c(beta_1_gee_MCAR_se, r[1, 5]) 
  beta_2_gee_MCAR_se <- c(beta_2_gee_MCAR_se, r[1, 6]) 
  beta_0_gee_MAR <- c(beta_0_gee_MAR, r[2, 1]) 
  beta_1_gee_MAR <- c(beta_1_gee_MAR, r[2, 2]) 
  beta_2_gee_MAR <- c(beta_2_gee_MAR, r[2, 3]) 
  beta_0_gee_MAR_se <- c(beta_0_gee_MAR_se, r[2, 4]) 
  beta_1_gee_MAR_se <- c(beta_1_gee_MAR_se, r[2, 5]) 
  beta_2_gee_MAR_se <- c(beta_2_gee_MAR_se, r[2, 6]) 
  beta_0_gee_MNAR <- c(beta_0_gee_MNAR, r[3, 1]) 
  beta_1_gee_MNAR <- c(beta_1_gee_MNAR, r[3, 2]) 
  beta_2_gee_MNAR <- c(beta_2_gee_MNAR, r[3, 3]) 
  beta_0_gee_MNAR_se <- c(beta_0_gee_MNAR_se, r[3, 4]) 
  beta_1_gee_MNAR_se <- c(beta_1_gee_MNAR_se, r[3, 5]) 
  beta_2_gee_MNAR_se <- c(beta_2_gee_MNAR_se, r[3, 6]) 
  beta_0_lmm_MCAR <- c(beta_0_lmm_MCAR, r[4, 1]) 
  beta_1_lmm_MCAR <- c(beta_1_lmm_MCAR, r[4, 2]) 
  beta_2_lmm_MCAR <- c(beta_2_lmm_MCAR, r[4, 3]) 
  beta_0_lmm_MCAR_se <- c(beta_0_lmm_MCAR_se, r[4, 4]) 
  beta_1_lmm_MCAR_se <- c(beta_1_lmm_MCAR_se, r[4, 5]) 
  beta_2_lmm_MCAR_se <- c(beta_2_lmm_MCAR_se, r[4, 6]) 
  beta_0_lmm_MAR <- c(beta_0_lmm_MAR, r[5, 1]) 
  beta_1_lmm_MAR <- c(beta_1_lmm_MAR, r[5, 2]) 
  beta_2_lmm_MAR <- c(beta_2_lmm_MAR, r[5, 3]) 
  beta_0_lmm_MAR_se <- c(beta_0_lmm_MAR_se, r[5, 4]) 
  beta_1_lmm_MAR_se <- c(beta_1_lmm_MAR_se, r[5, 5])
  beta_2_lmm_MAR_se <- c(beta_2_lmm_MAR_se, r[5, 6]) 
  beta_0_lmm_MNAR <- c(beta_0_lmm_MNAR, r[6, 1]) 
  beta_1_lmm_MNAR <- c(beta_1_lmm_MNAR, r[6, 2]) 
  beta_2_lmm_MNAR <- c(beta_2_lmm_MNAR, r[6, 3]) 
  beta_0_lmm_MNAR_se <- c(beta_0_lmm_MNAR_se, r[6, 4]) 
  beta_1_lmm_MNAR_se <- c(beta_1_lmm_MNAR_se, r[6, 5]) 
  beta_2_lmm_MNAR_se <- c(beta_2_lmm_MNAR_se, r[6, 6])
}



total_result <- data.frame(matrix(rep(0, 36), nrow = 6, ncol = 6))
colnames(total_result) <- c("beta_0 estimate", "beta_1 estimate", "beta_2 estimate", "beta_0 S.E.", "beta_1 S.E.",
                            "beta_2 S.E.")
rownames(total_result) <- c("gee_MCAR", "gee_MAR", "gee_MNAR"
                            , "lmm_MCAR", "lmm_MAR", "lmm_MNAR")


total_result[1, 1] <- mean(beta_0_gee_MCAR) 
total_result[1, 2] <- mean(beta_1_gee_MCAR) 
total_result[1, 3] <- mean(beta_2_gee_MCAR) 
total_result[1, 4] <- mean(beta_0_gee_MCAR_se) 
total_result[1, 5] <- mean(beta_1_gee_MCAR_se) 
total_result[1, 6] <- mean(beta_2_gee_MCAR_se) 
total_result[2, 1] <- mean(beta_0_gee_MAR) 
total_result[2, 2] <- mean(beta_1_gee_MAR) 
total_result[2, 3] <- mean(beta_2_gee_MAR) 
total_result[2, 4] <- mean(beta_0_gee_MAR_se) 
total_result[2, 5] <- mean(beta_1_gee_MAR_se) 
total_result[2, 6] <- mean(beta_2_gee_MAR_se) 
total_result[3, 1] <- mean(beta_0_gee_MNAR) 
total_result[3, 2] <- mean(beta_1_gee_MNAR) 
total_result[3, 3] <- mean(beta_2_gee_MNAR) 
total_result[3, 4] <- mean(beta_0_gee_MNAR_se) 
total_result[3, 5] <- mean(beta_1_gee_MNAR_se) 
total_result[3, 6] <- mean(beta_2_gee_MNAR_se) 
total_result[4, 1] <- mean(beta_0_lmm_MCAR) 
total_result[4, 2] <- mean(beta_1_lmm_MCAR) 
total_result[4, 3] <- mean(beta_2_lmm_MCAR) 
total_result[4, 4] <- mean(beta_0_lmm_MCAR_se) 
total_result[4, 5] <- mean(beta_1_lmm_MCAR_se) 
total_result[4, 6] <- mean(beta_2_lmm_MCAR_se) 
total_result[5, 1] <- mean(beta_0_lmm_MAR) 
total_result[5, 2] <- mean(beta_1_lmm_MAR) 
total_result[5, 3] <- mean(beta_2_lmm_MAR) 
total_result[5, 4] <- mean(beta_0_lmm_MAR_se) 
total_result[5, 5] <- mean(beta_1_lmm_MAR_se) 
total_result[5, 6] <- mean(beta_2_lmm_MAR_se) 
total_result[6, 1] <- mean(beta_0_lmm_MNAR) 
total_result[6, 2] <- mean(beta_1_lmm_MNAR)
total_result[6, 3] <- mean(beta_2_lmm_MNAR) 
total_result[6, 4] <- mean(beta_0_lmm_MNAR_se) 
total_result[6, 5] <- mean(beta_1_lmm_MNAR_se) 
total_result[6, 6] <- mean(beta_2_lmm_MNAR_se) 
total_result


