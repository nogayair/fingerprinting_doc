library(mice)
library(gee) 

# Load the NEW file path you provided
dat <- read.csv("/Volumes/homes/Noga/Papers/finngerprint/stab_disc_self_for_gee_290126.csv", stringsAsFactors = TRUE)

# 1. Clean and Setup
dat$id <- factor(dat$id)
dat$group <- factor(dat$group) # Assuming 1 and 2 based on your snippet
dat <- dat[order(dat$group, dat$id, dat$time),]

# SELECT OUTCOME: iself, iothers or discriminability
target_outcome <- "iothers" 

# 2. Split the data by group
# Adjusting subset logic to match your new file columns ('group' instead of 'groupname')
dat_g1 <- subset(dat, group == 1)
dat_g2 <- subset(dat, group == 2)

# 3. Impute separately
m <- 100
# Both use 'pmm' because these are continuous measures (discriminability, iself, iothers)
imputed_g1 <- mice(dat_g1, m = m, maxit = 50, method = 'pmm', seed = 500, printFlag = FALSE)
imputed_g2 <- mice(dat_g2, m = m, maxit = 50, method = 'pmm', seed = 500, printFlag = FALSE)

# 4. Fit GEE and extract coefficients/variance
fits <- lapply(1:m, function(i) {
  # Recombine the i-th imputed datasets
  combined_imputed_data <- rbind(complete(imputed_g1, action = i), 
                                 complete(imputed_g2, action = i))
  
  # Prepare formula dynamically
  form <- as.formula(paste(target_outcome, "~ group * time"))
  
  # Fit the model - Family is now GAUSSIAN for continuous data
  fit <- gee(form, 
             data = combined_imputed_data, 
             id = id, 
             family = gaussian, 
             corstr = "exchangeable")
  return(fit)
})

# 5. POOLING (Rubin's Rules)
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
print(paste("--- Pooled GEE Results for:", target_outcome, "---"))
print(summary_table)

# 6. POOLED GROUP AVERAGES
avg_results <- lapply(1:m, function(i) {
  comp_data <- rbind(complete(imputed_g1, action = i), 
                     complete(imputed_g2, action = i))
  
  # Calculate mean per cell for the continuous outcome
  form_agg <- as.formula(paste(target_outcome, "~ group + time"))
  aggregate(form_agg, data = comp_data, FUN = mean)
})

final_means <- Reduce(function(x, y) {
  res <- x
  res[[target_outcome]] <- x[[target_outcome]] + y[[target_outcome]]
  return(res)
}, avg_results)

final_means[[target_outcome]] <- final_means[[target_outcome]] / m
print(paste("--- Pooled Average Levels for:", target_outcome, "---"))
print(final_means)
