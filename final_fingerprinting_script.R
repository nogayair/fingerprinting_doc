library(mice)
library(gee) 

dat <- read.csv("/Volumes/homes/Noga/Papers/finngerprint/GEE/GEE_table4.csv", stringsAsFactors = TRUE)
dat$id <- factor(dat$id)
dat <- dat[order(dat$group, dat$id, dat$time),]
dat$groupname <- relevel(dat$groupname, ref = "Students")
dat$success <- factor(dat$success, levels = c(0, 1), labels = c("failure", "success"))

# 1. Split the data
dat_students <- subset(dat, groupname == "Students")
dat_soldiers <- subset(dat, groupname == "Soldiers")

# 2. Impute separately 
m <- 100
imputed_students <- mice(dat_students, m = m, maxit = 50, method = 'pmm', seed = 500)
imputed_soldiers <- mice(dat_soldiers, m = m, maxit = 50, method = 'logreg', seed = 500)

# 3. Fit GEE and extract coefficients/variance for pooling
fits <- lapply(1:m, function(i) {
  # Recombine the i-th imputed datasets
  combined_imputed_data <- rbind(complete(imputed_students, action = i), 
                                 complete(imputed_soldiers, action = i))
  
  combined_imputed_data$success_num <- as.numeric(combined_imputed_data$success) - 1  
  fit <- gee(success_num ~ groupname * time, 
             data = combined_imputed_data, 
             id = id, 
             family = binomial, 
             corstr = "exchangeable")
  return(fit)
})

# 4. POOLING (Applying Rubin's Rules)
all_coefs <- sapply(fits, function(x) coef(x))
all_vars <- lapply(fits, function(x) x$robust.variance)

pooled_coefs <- rowMeans(all_coefs)
W <- Reduce("+", all_vars) / m
B <- (1 / (m - 1)) * Reduce("+", lapply(1:m, function(i) (all_coefs[,i] - pooled_coefs) %*% t(all_coefs[,i] - pooled_coefs)))
total_var <- W + (1 + 1/m) * B
pooled_se <- sqrt(diag(total_var))

# Final GEE Results
summary_table <- data.frame(
  Estimate = pooled_coefs,
  Std.Error = pooled_se,
  z_value = pooled_coefs / pooled_se,
  p_value = 2 * (1 - pnorm(abs(pooled_coefs / pooled_se)))
)
print("--- Pooled GEE Results ---")
print(summary_table)

# 5. POOLED GROUP AVERAGES
avg_success_rates <- lapply(1:m, function(i) {
  # Recombine
  comp_data <- rbind(complete(imputed_students, action = i), 
                     complete(imputed_soldiers, action = i))
  comp_data$success_num <- as.numeric(comp_data$success) - 1  
  
  # Calculate mean per cell
  aggregate(success_num ~ groupname + time, data = comp_data, FUN = mean)
})

# Average the rates across the 100 sets
final_rates <- Reduce(function(x, y) {
  res <- x
  res$success_num <- x$success_num + y$success_num
  return(res)
}, avg_success_rates)

final_rates$success_num <- final_rates$success_num / m
print("--- Pooled Average Success Rates ---")
print(final_rates)
