fitDABB <- function(Y, tot_read, bat_ind, bio_ind, gene_len, tech_para){
  G_num <- length(Y[ ,1])
  C_num <- length(Y[1, ])
  bat_num <- max(bat_ind)
  bio_num <- max(bio_ind)

  Y_binary <- ifelse(Y > 0, 1, 0)
  bat_mat <- c()
  for (i in 1:bat_num){
    bat_mat <- cbind(bat_mat, ifelse(bat_ind == i, 1, 0))
  }
  bat_cellnum <- colSums(bat_mat)
  bio_mat <- c()
  for (i in 1:bio_num){
    bio_mat <- cbind(bio_mat, ifelse(bio_ind == i, 1, 0))
  }
  
  RL_mat <- gene_len %*% t(tot_read) / 10^9
  
  lnR_vec <- c()
  for (i in 1:C_num){
    lnR_vec <- c(lnR_vec, rep(log(tot_read[i] / 10^6), G_num))
  }
  lnL_vec <- rep(log(gene_len / 10^3), C_num)
  
  Design_mat <- cbind(rep(1, G_num * C_num), lnR_vec, lnL_vec)
  Design_cov <- t(Design_mat) %*% Design_mat
  
  iter_num <- tech_para$iternum
  max_error <- tech_para$error
  MH_num <- tech_para$mhnum
  sig_jump <- tech_para$jump
  dec_rate <- tech_para$jump_rate
  burn_num <- tech_para$burnin
  sample_num <- MH_num - burn_num
  
  iter <- 0
  loglike_old <- 1
  nu_Bat_old <- 0
  error <- 1
  ##sigma2_Bat_old <- 1
  
  #initial value
  nu_Bat <- rep(0, bat_num)
  sigma2_Bat <- rep(1, bat_num)
  coef <- c(0, 0.1, 0.1, 1)
  mu_mat <- mu_update(Y, bio_mat, RL_mat, matrix(1, G_num, C_num), rep(1, C_num))
  
  while (iter < iter_num & error > max_error){
    iter <- iter + 1
    #E-step
    if (iter > 1 & iter < 3){
      sig_jump <- dec_rate * sig_jump
    }
    
    b_sample_mat <- c()
    loglike_mat <- c()
    expect_eta <- 0
    expect_eta_b <- 0
    p_weight <- 0
    Phi_mat <- 0
    ##loglike_all <- 0
    jump_mat <- c()
    
    b_old <- rnorm(C_num, as.vector(bat_mat %*% nu_Bat), 
                   as.vector(bat_mat %*% sigma2_Bat^0.5))
    lst_post_old <- b_post_loglike(Y, b_old, mu_mat, bio_mat, RL_mat, Y_binary,
                                   coef, Design_mat, bat_mat, nu_Bat, sigma2_Bat)
    log_post_old <- lst_post_old$loglike1
    
    log_post_bii <- lst_post_old$loglike2
    p_weightii <- lst_post_old$pmat
    expect_etaii <- lst_post_old$etaexpect
    Phi_matii <- lst_post_old$Phimat
    ##loglike_allii <- lst_post_old$loglikemat
      
    for (ii in 1:MH_num){
      #Metropolis-Hasting
      b_new <- rnorm(C_num, b_old, rep(sig_jump, C_num))
      lst_post_new <- b_post_loglike(Y, b_new, mu_mat, bio_mat, RL_mat, Y_binary,
                                     coef, Design_mat, bat_mat, nu_Bat, sigma2_Bat)
      log_post_new <- lst_post_new$loglike1
      log_post_diff <- log_post_new - log_post_old
      r <- ifelse(log_post_diff >= 0, 1, exp(log_post_diff))
      unif <- runif(C_num)
      b_old <- ifelse(r > unif, b_new, b_old)
      log_post_old <- ifelse(r > unif, log_post_new, log_post_old)
      jump_vec <- as.vector(ifelse(r > unif, 1, 0))
      jump_mat <- cbind(jump_mat, jump_vec)
      Phi_matii <- t((1 - jump_vec) * t(Phi_matii) + jump_vec * t(lst_post_new$Phimat))
      log_post_bii <- (1 - jump_vec) * log_post_bii + jump_vec * lst_post_new$loglike2
      p_weightii <- t((1 - jump_vec) * t(p_weightii) + jump_vec * t(lst_post_new$pmat))
      expect_etaii <- t((1 - jump_vec) * t(expect_etaii) +
                          jump_vec * t(lst_post_new$etaexpect))
      ##loglike_allii <- t((1 - jump_vec) * t(loglike_allii) +
        ##                   jump_vec * t(lst_post_new$loglikemat))
      
      #Monte Carlo
      if (ii > burn_num){
        b_sample_mat <- cbind(b_sample_mat, b_old)
        expect_eta <- expect_eta + expect_etaii
        expect_eta_b <- expect_eta_b + t(t(expect_etaii) * as.vector(b_old))
        p_weight <- p_weight + p_weightii
        loglike_mat <- cbind(loglike_mat, log_post_bii)
        Phi_mat <- Phi_mat + Phi_matii
        ##loglike_all <- loglike_all + loglike_allii
      }
    }
    expect_b <- rowMeans(b_sample_mat)
    expect_expb <- rowMeans(exp(b_sample_mat))
    expect_b2 <- rowMeans(b_sample_mat^2)
    expect_eta <- expect_eta / sample_num
    expect_eta_b <- expect_eta_b / sample_num
    p_weight <- p_weight / sample_num
    loglike_vec <- rowMeans(loglike_mat)
    loglike_new <- sum(loglike_vec)
    Phi_mat <- Phi_mat / sample_num
    ##loglike_all <- loglike_all / sample_num
    
    print(mean(jump_mat))
    
    #M-step
    mu_mat <- mu_update(Y, bio_mat, RL_mat, p_weight, expect_expb)
    A1_12 <- t(Design_mat) %*% as.vector(rep(1, G_num) %*% t(expect_b))
    A1_22 <- G_num * sum(expect_b2)
    A1 <- cbind(rbind(Design_cov, t(A1_12)), rbind(A1_12, A1_22))
    A2 <- rbind(t(Design_mat) %*% as.matrix(as.vector(expect_eta)), 
                sum(expect_eta_b)) 
    coef <- solve(A1) %*% A2
    expect_b_shift <- expect_b - mean(expect_b)

    nu_Bat <- (t(bat_mat) %*% expect_b_shift) / bat_cellnum
    sigma2_Bat <- (t(bat_mat) %*% ((expect_b - bat_mat %*% nu_Bat)^2)) / bat_cellnum
    
    print(c(iter, loglike_new))
    
    error1 <- abs((loglike_new - loglike_old) / (loglike_old + 1))
    loglike_old <- loglike_new
    error2 <- 1e3 * mean((nu_Bat - nu_Bat_old)^2 / (nu_Bat_old + 1))
    nu_Bat_old <- nu_Bat
    error <- min(error1, error2)
    ##print(c(error1, error2))
  }
  ##loglike_gene <- rowSums(loglike_all)
  return(list(nu = nu_Bat, sigma2 = sigma2_Bat, mu = mu_mat, bsample = b_sample_mat, 
              pweight = p_weight, coef = coef, eta = expect_eta, Phi = Phi_mat))
}

#update functions

mu_update <- function(Y, bio_mat, RL_mat, p_weight, exp_b){
  return((p_weight * Y) %*% bio_mat / (t(t(p_weight * RL_mat) * exp_b) %*% bio_mat))
}

b_post_loglike <- function(Y, b_vec, mu_mat, bio_mat, RL_mat, Y_binary,
                           coef, Design_mat, bat_mat, nu_Bat, sigma2_Bat){
  G_num <- length(Y_binary[,1])
  C_num <- length(Y_binary[1,])
  lambda_mat <- t(t(mu_mat %*% t(bio_mat) * RL_mat) * exp(b_vec))
  eta_mean_vec <- as.vector(cbind(Design_mat, b_vec) %*% coef)
  eta_mean <- matrix(eta_mean_vec, G_num, C_num)
  p_nodrop_mat <- matrix(pnorm(eta_mean_vec), G_num, C_num)
  zero_poisson <- p_nodrop_mat * (1 - Y_binary) * exp(-lambda_mat)
  zero_prob <- zero_poisson + (1 - p_nodrop_mat) * (1 - Y_binary)
  loglike_mat <- Y * log((1e-10) + lambda_mat) - lambda_mat * Y_binary + 
    log(p_nodrop_mat * Y_binary + zero_prob)
  loglike_cell <- colSums(loglike_mat)
  b_loglike_vec <- - 1 / (2 * bat_mat %*% sigma2_Bat) * (b_vec - bat_mat %*% nu_Bat)^2
  p_weight <- Y_binary + zero_poisson / (zero_prob + Y_binary) 
  phi_mat <- dnorm(eta_mean)
  expect_eta <- eta_mean + p_weight * phi_mat / p_nodrop_mat - 
    (1 - p_weight) * phi_mat / (1 - p_nodrop_mat)
  return(list(loglike1 = loglike_cell + b_loglike_vec, loglike2 = loglike_cell, 
              etaexpect = expect_eta, pmat = p_weight, Phimat = p_nodrop_mat,
              loglikemat = loglike_mat))
}

#Quality control
DABB_QC <- function(results, alternative = 'right', level = 0.05){
  b_sample_mat <- results$bsample
  sample_num <- length(b_sample_mat)
  C_num <- length(b_sample_mat[,1])
  S_num <- length(b_sample_mat[1,])
  pvalue_lst <- c()
  if (alternative == 'right'){
    level_ind <- as.integer((1 - level) * sample_num)
    up_bound <- sort(b_sample_mat)[level_ind]
    indicate_mat <- ifelse(b_sample_mat >= up_bound, 1, 0)
    sig_num <- rowSums(indicate_mat)
    for (i in 1:C_num){
      cont_table <- cbind(c(as.integer(sample_num * level), sample_num), 
                          c(sig_num[i], S_num))
      ##print(cont_table)
      p_value <- fisher.test(cont_table, alternative = 'less')$p.value
      pvalue_lst <- c(pvalue_lst, p_value)
    }
  }
  
  if (alternative == 'left'){
    level_ind <- as.integer(level * sample_num)
    low_bound <- sort(b_sample_mat)[level_ind]
    indicate_mat <- ifelse(b_sample_mat <= low_bound, 1, 0)
    sig_num <- rowSums(indicate_mat)
    for (i in 1:C_num){
      cont_table <- cbind(c(as.integer(sample_num * level), sample_num), 
                          c(sig_num[i], S_num))
      p_value <- fisher.test(cont_table, alternative = 'less')$p.value
      pvalue_lst <- c(pvalue_lst, p_value)
    }
  }
  
  if (alternative == 'two side'){
    level_up <- as.integer((1 - level / 2) * sample_num)
    level_low <- as.integer(level / 2 * sample_num)
    up_bound <- sort(b_sample_mat)[level_up]
    low_bound  <- sort(b_sample_mat)[level_low]
    indicate_mat <- ifelse(b_sample_mat < low_bound | b_sample_mat >= up_bound, 1, 0)
    sig_num <- rowSums(indicate_mat)
    for (i in 1:C_num){
      cont_table <- cbind(c(as.integer(sample_num * level), sample_num), 
                          c(sig_num[i], S_num))
      p_value <- fisher.test(cont_table, alternative = 'less')$p.value
      pvalue_lst <- c(pvalue_lst, p_value)
    }
  }
  return(pvalue_lst)
}


#Differential Expression
Random_assign <- function(bio_ind){
  num1 <- sum(ifelse(bio_ind == 1, 1, 0))
  num2 <- sum(ifelse(bio_ind == 2, 1, 0))
  C_num <- num1 + num2
  fake_group1 <- sample.int(C_num, size = num1)
  fake_mat <- matrix(0, C_num, 2)
  fake_mat[fake_group1, 1] <- 1
  fake_mat[as.vector(fake_mat[,1]) == 0, 2] <- 1
  return(fake_mat)
}

cal_loglike <- function(Y, Y_binary, RL_mat, bio_mat, p_mat, expect_expb, p_weight,
                 max_iter = 30, max_error = 1e-3){
  iter <- 0
  error <- 1
  G_num <- length(Y[,1])
  C_num <- length(Y[1,])
  loglike_old <- rep(-Inf, G_num)
  while (iter < max_iter & error > max_error){
    iter <- iter + 1
    #M-step
    mu_mat <- mu_update(Y, bio_mat, RL_mat, p_weight, expect_expb)
    lambda_mat <- t(t(mu_mat %*% t(bio_mat) * RL_mat) * expect_expb)
    #E-step
    zero_poisson <- p_mat * (1 - Y_binary) * exp(-lambda_mat)
    zero_prob <- zero_poisson + (1 - p_mat) * (1 - Y_binary)
    p_weight <- Y_binary + zero_poisson / (zero_prob + Y_binary) 
    loglike_mat <- Y * log((1e-10) + lambda_mat) - lambda_mat * Y_binary + 
      log(p_mat * Y_binary + zero_prob)
    loglike_gene <- rowSums(loglike_mat)
    error_vec <- ifelse(loglike_gene > loglike_old, 
                        abs(loglike_gene - loglike_old) / (1 + abs(loglike_gene)), 0)
    loglike_old <- ifelse(loglike_gene > loglike_old, loglike_gene, loglike_old)
    error <- max(error_vec)
  }
  ##print(iter)
  return(loglike_old)
}

DABB_DE <- function(Y, tot_read, bio_ind, gene_len, results, sample_num = 50){
  RL_mat <- gene_len %*% t(tot_read) / 10^9
  b_sample_mat <- results$bsample
  p_nodrop_mat <- results$Phi
  p_weight <- results$pweight
  expect_expb <- rowMeans(exp(b_sample_mat))
  C_num <- length(Y[1, ])
  Y_binary <- ifelse(Y > 0, 1, 0)
  bio_mat <- c()
  for (i in 1:max(bio_ind)){
    bio_mat <- cbind(bio_mat, ifelse(bio_ind == i, 1, 0))
  }
  loglike_null <- cal_loglike(Y, Y_binary, RL_mat, rep(1, C_num), p_nodrop_mat, expect_expb, p_weight,
                              max_iter = 30, max_error = 1e-3)
  diff_loglike_true <- 2 * 
    (cal_loglike(Y, Y_binary, RL_mat, bio_mat, p_nodrop_mat, expect_expb, p_weight,
                 max_iter = 40, max_error = 1e-3) - loglike_null)
  diff_loglike_gene_mat <- c()
  for (ii in 1:sample_num){
    fake_mat <- Random_assign(bio_ind)
    diff_loglike <- 2 *
      (cal_loglike(Y, Y_binary, RL_mat, fake_mat, p_nodrop_mat, expect_expb, p_weight,
                   max_iter = 40, max_error = 1e-3) - loglike_null)
    diff_loglike_gene_mat <- cbind(diff_loglike_gene_mat, diff_loglike)
    print(ii)
  }
  diff_loglike_mean <- rowMeans(diff_loglike_gene_mat)
  Test_Stat <- diff_loglike_true / diff_loglike_mean
  p_value <- pchisq(Test_Stat, df = 1)
  return(list(p.value = p_value, loglike = diff_loglike_gene_mat))
}

#Visualization
DABB_visualize <- function(Y, tot_read, gene_len, p_weight, b_sample_mat, 
                           mdim = 2, method = 'PCA', k = 5){
  RL_mat <- gene_len %*% t(tot_read) / 10^9
  expect_expb <- rowMeans(exp(b_sample_mat))
  ##Y_adj <- log(t(t(Y / RL_mat) / expect_expb) + 1)
  Y_adj <- t(t(Y / RL_mat) / expect_expb)
  ##Y_adj <- t(t(log(Y / RL_mat + 1)) / expect_expb)
  Y_shift <- t(t(Y_adj) - colSums(Y_adj * p_weight) / colSums(p_weight))
  Y_weight <- Y_shift * p_weight
  ##var_vec <- diag(Y_weight %*% t(Y_shift)) / rowSums(p_weight)
  cov_mat <- (t(Y_weight) %*% Y_weight) / (t(p_weight) %*% p_weight)
  if (method == 'PCA'){
    PCA <- svd(cov_mat)
    U <- PCA$u[,c(1:mdim)]
    return(U)
  }
    
  if (method == 'ISOmap'){
    C_num <- length(Y[1,])
    Dissimilar <- as.dist(diag(cov_mat) %*% t(rep(1, C_num)) + 
      t(diag(cov_mat) %*% t(rep(1, C_num))) -  2 * cov_mat)
    return(isomap(Dissimilar, mdim, k = k))
  }
  
}

