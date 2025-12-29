library(spatialprobit)
library(Matrix)
library(dplyr)
library(splines)

happy <- read.csv("happying_english.csv", stringsAsFactors = FALSE)

# Y variable is numeric
happy$top100 <- as.numeric(happy$top100)

# Create IDs for counties and years, then sort (year, county)
counties <- sort(unique(happy$county_name))
years    <- sort(unique(happy$year))

happy <- happy %>%
  mutate(
    county_id = match(county_name, counties),
    year_id   = match(year, years)
  ) %>%
  arrange(year_id, county_id)

N  <- length(counties)   # number of counties
TT <- length(years)      # number of years
stopifnot(N * TT == nrow(happy))

# Cross-sectional W_N (N x N) sparse
W_N <- Matrix(0, nrow = N, ncol = N, sparse = TRUE)

for (i in 1:N) {
  city_i <- happy$city[happy$county_id == i][1]
  for (j in 1:N) {
    if (i != j) {
      city_j <- happy$city[happy$county_id == j][1]
      if (!is.na(city_i) && !is.na(city_j) && city_i == city_j) {
        W_N[i, j] <- 1
      }
    }
  }
}

# Zero diagonal
diag(W_N) <- 0

# Row-standardize: rows sum to 1 when neighbors exist
rs <- rowSums(W_N)
rs[rs == 0] <- 1
W_N <- Diagonal(x = 1 / rs) %*% W_N

# Panel W = I_T x W_N (block-diagonal, sparse)
I_T     <- Diagonal(TT)
W_panel <- kronecker(I_T, W_N)

selected_vars <- c("GDP", "registered_pop_10k", "cy", "jr")

Xvars <- selected_vars
# 1) Extract predictors as numeric matrix
X_no_int <- as.matrix(happy[, Xvars])

# Make sure everything is numeric
X_no_int <- apply(X_no_int, 2, as.numeric)

# 2) SCALE the predictors: center = TRUE, scale = TRUE
X_scaled <- scale(X_no_int, center = TRUE, scale = TRUE)

# 3) Add intercept AFTER scaling
X <- cbind(Intercept = 1, X_scaled)
colnames(X)[1] <- "Intercept"

# Response
y <- happy$top100

set.seed(123)

# Linear Spatial Probit Model
######################################
fit_sar <- sar_probit_mcmc(
  y      = y,
  X      = X,
  W      = W_panel,
  ndraw  = 10000,
  burn.in = 2000,
  thinning = 5,
  prior  = NULL,
  showProgress = TRUE
)

# 'fit_sar' is an object of class 'sarprobit'
summary(fit_sar)
# diagnostic plots: chains, densities, ACF, etc.
plot(fit_sar)

mean(fit_sar$rho) #p rho = 0.619

rho_mean <- 6.187e-01
rho_sd   <- 4.024e-02

rho_ci_norm <- rho_mean + qnorm(c(0.025, 0.975)) * rho_sd
rho_ci_norm #[0.540, 0.698]

intercept_ci_norm <- -0.26051 + qnorm(c(0.025, 0.975)) * 0.06785
gdp_ci_norm <- 0.31924 + qnorm(c(0.025, 0.975)) * 0.08762
pop_ci_norm <- 0.43169 + qnorm(c(0.025, 0.975)) * 0.08662
cy_ci_norm <- 0.35673 + qnorm(c(0.025, 0.975)) * 0.08540
jr_ci_norm <- 0.19793 + qnorm(c(0.025, 0.975)) * 0.07576

intercept_ci_norm #(-0.393, -0.128)
gdp_ci_norm #(0.148, 0.491)
pop_ci_norm #(0.262, 0.601)
cy_ci_norm #(0.189, 0.524)
jr_ci_norm #(0.049, 0.346)


#Varying Coefficient Spatial Probit Model
#################################################
# Response is numeric
happy$year <- as.numeric(happy$year)

# Scale time to [0, 1]
t_min <- min(happy$year)
t_max <- max(happy$year)

happy$t_scaled <- (happy$year - t_min) / (t_max - t_min)

# Number of splines (B1)
K <- 1
#1 Degree
B <- bs(happy$t_scaled, df = K, degree = 1)
B <- as.matrix(B)
colnames(B) <- paste0("B", 1:K)

# Attach to data frame
for (k in 1:K) {
  happy[[paste0("B", k)]] <- B[, k]
}

spline_terms <- paste0("B", 1:K)

vars_static <- c("GDP", "registered_pop_10k", "cy", "jr")

# Extract and coerce to numeric
X_static_raw <- happy[, vars_static, drop = FALSE]
X_static_raw[] <- lapply(X_static_raw, as.numeric)

# Scale for numerical stability
X_static_scaled <- scale(as.matrix(X_static_raw))
colnames(X_static_scaled) <- vars_static

# GDP as numeric
#happy$GDP <- as.numeric(happy$GDP)

# Scale GDP
X_vc <- cbind(
  Intercept = 1,
  X_static_scaled,                   # constant effects
  as.matrix(happy[, spline_terms])   # varying intercept basis
)

X_vc <- as.matrix(X_vc)
storage.mode(X_vc) <- "double"

fit_vc <- sar_probit_mcmc(
  y      = y,
  X      = X_vc,
  W      = W_panel,
  ndraw  = 10000,
  burn.in = 2000,
  thinning = 5,
  prior  = NULL,       # or your prior list
  showProgress = TRUE
)

summary(fit_vc)

plot(fit_vc)

fit_vc$bdraw
#p rho = 0.620
rho_mean_vc <- 6.19e-01
rho_sd_vc   <- 4.0128e-02

rho_ci_norm_vc <- rho_mean_vc + qnorm(c(0.025, 0.975)) * rho_sd_vc
rho_ci_norm_vc #[0.540, 0.698]

int_ci_norm_vc <- 0.14883 + qnorm(c(0.025, 0.975)) * 0.16784
int_ci_norm_vc #[-0.180, 0.478]

gdp_ci_norm_vc <- 0.47206 + qnorm(c(0.025, 0.975)) * 0.12366
gdp_ci_norm_vc #[0.230, 0.714]

pop_ci_norm_vc <- 0.46081 + qnorm(c(0.025, 0.975)) * 0.09051
pop_ci_norm_vc #[0.283, 0.638]

cy_ci_norm_vc <- 0.20094 + qnorm(c(0.025, 0.975)) * 0.09928
cy_ci_norm_vc #[0.006, 0.396]

jr_ci_norm_vc <- 0.28026 + qnorm(c(0.025, 0.975)) * 0.09023
jr_ci_norm_vc #[0.103, 0.457]

b1_ci_norm_vc <- -0.79404 + qnorm(c(0.025, 0.975)) * 0.32471
b1_ci_norm_vc #[-1.430, -0.158]

##################
library(pROC)

predict_prob_probit_approx <- function(fit, X) {
  beta_hat <- colMeans(fit$bdraw)
  eta_hat  <- as.vector(X %*% beta_hat)
  p_hat    <- pnorm(eta_hat)
  pmin(pmax(p_hat, 1e-12), 1 - 1e-12)
}

p_lin <- predict_prob_probit_approx(fit_sar, X)
roc_lin <- roc(response = y, predictor = p_lin)
auc_lin <- auc(roc_lin) #0.7833
plot(roc_lin, plot = TRUE, col = "blue", print.auc = TRUE, main = "Linear")

p_vc <- predict_prob_probit_approx(fit_vc, X_vc)
roc_vc <- roc(response = y, predictor = p_vc)
auc_vc <- auc(roc_vc) #0.791

plot(roc_vc, plot = TRUE, col = "blue", print.auc = TRUE, main = "Varying Coefficient")
#Follow-up
############################################
compute_dic_probit_approx <- function(fit, X, y, ndraw_eval = 1000) {
  # X: design matrix used in fitting (with intercept)
  # y: 0/1 outcome vector
  
  nd <- min(ndraw_eval, nrow(fit$bdraw))
  idx <- seq(1, nrow(fit$bdraw), length.out = nd)
  
  dev_vec <- numeric(nd)
  
  for (k in seq_along(idx)) {
    j <- idx[k]
    beta_j <- fit$bdraw[j, ]
    
    eta_j <- as.vector(X %*% beta_j)
    p_j   <- pnorm(eta_j)
    p_j   <- pmin(pmax(p_j, 1e-12), 1 - 1e-12)
    
    loglik_j <- sum(y * log(p_j) + (1 - y) * log(1 - p_j))
    dev_vec[k] <- -2 * loglik_j
  }
  
  # Posterior mean deviance
  D_bar <- mean(dev_vec)
  
  # Deviance at posterior mean beta
  beta_hat <- colMeans(fit$bdraw)
  eta_hat  <- as.vector(X %*% beta_hat)
  p_hat    <- pnorm(eta_hat)
  p_hat    <- pmin(pmax(p_hat, 1e-12), 1 - 1e-12)
  
  loglik_hat <- sum(y * log(p_hat) + (1 - y) * log(1 - p_hat))
  D_hat <- -2 * loglik_hat
  
  # DIC
  DIC <- 2 * D_bar - D_hat
  
  list(D_bar = D_bar, D_hat = D_hat, DIC = DIC)
}

DIC_lin <- compute_dic_probit_approx(fit_sar, X, y)
DIC_lin$DIC #691.08

DIC_vc <- compute_dic_probit_approx(fit_vc, X_vc, y)
DIC_vc$DIC #680.78

#Plots
##########################################################
#Trace plots for SAR and VC
par(mfrow=c(2,5), oma=c(0,0,2,0))
plot(fit_sar$bdraw[,1], type = "l",
     main = "GDP",
     xlab = "",
     ylab = expression(bold("SAR")))
plot(fit_sar$bdraw[,2], type = "l",
     main = "Population",
     xlab = "",
     ylab = "")
plot(fit_sar$bdraw[,3], type = "l",
     main = "cy",
     xlab = "",
     ylab = "")
plot(fit_sar$bdraw[,4], type = "l",
     main = "jr",
     xlab = "",
     ylab = "")
plot(fit_sar$pdraw, type = "l",
     main = "rho",
     xlab = "",
     ylab = "")
plot(fit_vc$bdraw[,1], type = "l",
     main = "",
     xlab = "covariate Value",
     ylab = expression(bold("VC")))
plot(fit_vc$bdraw[,2], type = "l",
     main = "",
     xlab = "covariate Value",
     ylab = "")
plot(fit_vc$bdraw[,3], type = "l",
     main = "",
     xlab = "covariate Value",
     ylab = "")
plot(fit_vc$bdraw[,4], type = "l",
     main = "",
     xlab = "covariate Value",
     ylab = "")
plot(fit_vc$pdraw, type = "l",
     main = "",
     xlab = "covariate Value",
     ylab = "")
mtext("Trace Plot Comparison of SAR vs. VC", 
      outer = TRUE, cex = 1.2)
#Density plots for SAR and VC
par(mfrow=c(2,5), oma=c(0,0,2,0))
plot(density(fit_sar$bdraw[,1]),
     main = "GDP",
     xlab = "",
     ylab = expression(bold("SAR")))
plot(density(fit_sar$bdraw[,2]),
     main = "Population",
     xlab = "",
     ylab = "")
plot(density(fit_sar$bdraw[,3]),
     main = "cy",
     xlab = "",
     ylab = "")
plot(density(fit_sar$bdraw[,4]),
     main = "jr",
     xlab = "",
     ylab = "")
plot(density(fit_sar$pdraw),
     main = "rho",
     xlab = "",
     ylab = "")
plot(density(fit_vc$bdraw[,1]),
     main = "",
     xlab = "covariate Value",
     ylab = expression(bold("VC")))
plot(density(fit_vc$bdraw[,2]),
     main = "",
     xlab = "covariate Value",
     ylab = "")
plot(density(fit_vc$bdraw[,3]),
     main = "",
     xlab = "covariate Value",
     ylab = "")
plot(density(fit_vc$bdraw[,4]),
     main = "",
     xlab = "covariate Value",
     ylab = "")
plot(density(fit_vc$pdraw),
     main = "",
     xlab = "covariate Value",
     ylab = "")
mtext("Density Plot Comparison of SAR vs. VC", 
      outer = TRUE, cex = 1.2)
#Diagnostics
##########################################################
X <- happy[, Xvars]
X[] <- lapply(X, as.numeric)

drop_vars <- c("top100", "county_name", "city", "province", "county_id", "year_id")
drop_vars <- intersect(drop_vars, names(happy))

Xvars     <- setdiff(names(happy), drop_vars)

Xvars # check predictors

round(cor(X_scaled), 2)
# correlation matrix
cor_X <- cor(X, use = "pairwise.complete.obs")
round(cor_X, 2)

library(car)

fm2 <- lm(
  as.numeric(top100) ~ .,
  data = cbind(top100 = happy$top100, happy[, selected_vars])
)

vif_values2 <- vif(fm2)

sort(vif_values2, decreasing = TRUE)

#Scaling the selected variables and checking for multicollinearity
sel_vars <- happy[, selected_vars]              # your chosen predictors
sel_vars[] <- lapply(sel_vars, as.numeric)     # make sure numeric
sel_vars_scaled <- scale(sel_vars)

sel_vars_scaled_df <- as.data.frame(sel_vars_scaled)

lm_scaled_2 <- lm(
  as.numeric(top100) ~ .,
  data = cbind(top100 = happy$top100, sel_vars_scaled_df)
)

vif_vals2 <- vif(lm_scaled_2)
vif_vals2