# =============================================================================
#
# Closed-form, vectorised two-way mixed ANOVA simulation for the
# donor-by-treatment interaction extension.
#
# Design per sim: a donors x 2 conditions x k organoids per cell, balanced.
# Model:
#   y_{dck} = alpha_c + u_d + v_{dc} + eps_{dck}
#     u_d ~ N(0, sigma_donor^2)   donor baseline
#     v_{dc} ~ N(0, sigma_int^2) donor-by-condition interaction
#     eps ~ N(0, sigma_res^2)
#     alpha  fixed: 0 (control), beta (treated)
#
# Expected mean squares (unrestricted mixed model):
#   E[MS_C] = sigma_res^2 + k*sigma_int^2 + a*k*Q(alpha)
#   E[MS_D] = sigma_res^2 + k*sigma_int^2 + 2*k*sigma_donor^2
#   E[MS_DC] = sigma_res^2 + k*sigma_int^2
#   E[MS_E]= sigma_res^2
#
# Correct F-tests:
#   condition fixed effect: F = MS_C  / MS_DC~ F(1,   a-1)
#   sigma_int^2 > 0:   F = MS_DC / MS_E ~ F(a-1, 2a(k-1))
#   sigma_donor^2 > 0:F = MS_D  / MS_DC ~ F(a-1, a-1)
#
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(20260416)
# -----------------------------------------------------------------------------
# Core simulator: vectorised across n_sims
# -----------------------------------------------------------------------------
twoway_anova_fast <- function(a, k, beta,
                              sigma_donor, sigma_int, sigma_res,
                              n_sims) {
  cc <- 2L
  N <- a * cc * k
  
  # --- Random draws ---
  u <- matrix(rnorm(a * n_sims, 0, sigma_donor), a, n_sims)
  v <- array(if (sigma_int > 0)
    rnorm(a * cc * n_sims, 0, sigma_int)
    else rep(0, a * cc * n_sims),
    c(a, cc, n_sims))
  eps <- array(rnorm(k * a * cc * n_sims, 0, sigma_res),
               c(k, a, cc, n_sims))
  
  # Fixed condition effect
  alpha_arr <- array(0, c(k, a, cc, n_sims))
  alpha_arr[, , 2L, ] <- beta
  
  # Broadcast u -> (k, a, cc, n_sims)
  u_arr <- aperm(array(u, c(a, n_sims, k, cc)), c(3, 1, 4, 2))
  
  # Broadcast v -> (k, a, cc, n_sims)
  v_arr <- aperm(array(v, c(a, cc, n_sims, k)), c(4, 1, 2, 3))
  
  y <- alpha_arr + u_arr + v_arr + eps # (k,a,cc,n_sims)
  
  # --- Means per sim ---
  grand_means <- colMeans(array(y, c(k * a * cc, n_sims)))  # n_sims
  
  cell_means <- array(colMeans(array(y, c(k, a * cc * n_sims))),
                       c(a, cc, n_sims))  # (a,cc,n_sims)
  
  y_kc <- aperm(y, c(1, 3, 2, 4))  # (k,cc,a,n_sims)
  donor_means <- matrix(colMeans(array(y_kc, c(k * cc, a * n_sims))),
                        a, n_sims)# (a,n_sims)
  
  cond_means <- matrix(colMeans(array(y, c(k * a, cc * n_sims))),
                        cc, n_sims) # (cc,n_sims)
  
  # --- Sums of squares ---
  cm_dev <- cond_means  - matrix(grand_means, cc, n_sims, byrow = TRUE)
  ss_c <- a * k  * colSums(cm_dev^2)
  
  dm_dev <- donor_means - matrix(grand_means, a,  n_sims, byrow = TRUE)
  ss_d <- cc * k * colSums(dm_dev^2)
  
  dm_bc <- aperm(array(donor_means, c(a, n_sims, cc)), c(1, 3, 2))
  cm_bc <- aperm(array(cond_means, c(cc, n_sims, a)), c(3, 1, 2))
  gm_bc <- aperm(array(grand_means, c(n_sims, a, cc)), c(2, 3, 1))
  int_resid <- cell_means - dm_bc - cm_bc + gm_bc
  ss_dc <- k * colSums(array(int_resid, c(a * cc, n_sims))^2)
  
  cell_means_bc <- aperm(array(cell_means, c(a, cc, n_sims, k)),
                         c(4, 1, 2, 3))
  ss_e <- colSums(array((y - cell_means_bc)^2, c(N, n_sims)))
  
  # --- Mean squares, F tests ---
  df_c <- cc - 1L
  df_d <- a  - 1L
  df_dc <- (a - 1L) * (cc - 1L)
  df_e <- a * cc * (k - 1L)
  
  ms_c <- ss_c / df_c
  ms_d <- ss_d / df_d
  ms_dc <- ss_dc / df_dc
  ms_e <- ss_e / df_e
  
  f_cond <- ms_c / ms_dc; p_cond <- pf(f_cond, df_c,  df_dc, lower.tail = FALSE)
  f_int <- ms_dc / ms_e;  p_int <- pf(f_int,  df_dc, df_e,  lower.tail = FALSE)
  f_don <- ms_d / ms_dc; p_don <- pf(f_don,  df_d,  df_dc, lower.tail = FALSE)
  
  # Method-of-moments variance components (truncated at zero)
  sigma2_int   <- pmax((ms_dc - ms_e) / k,0)
  sigma2_donor <- pmax((ms_d  - ms_dc) / (cc * k), 0)
  
  # Direction-concordance diagnostic on observed per-donor effects
  donor_effect <- cell_means[, 2L, , drop = TRUE] - cell_means[, 1L, , drop = TRUE]
  if (a == 1L) donor_effect <- matrix(donor_effect, nrow = 1L)
  n_pos <- colSums(donor_effect > 0)
  prop_same <- pmax(n_pos, a - n_pos) / a
  all_same_sign <- (n_pos == 0L) | (n_pos == a)
  
  list(
    p_cond= p_cond,
    p_interaction = p_int,
    p_donor = p_don,
    sigma_int_hat = sqrt(sigma2_int),
    sigma_donor_hat = sqrt(sigma2_donor),
    sigma_res_hat = sqrt(ms_e),
    prop_same_sign = prop_same,
    all_same_sign = all_same_sign,
    ms_c = ms_c, ms_d = ms_d, ms_dc = ms_dc, ms_e = ms_e
  )
}

# -----------------------------------------------------------------------------
# Null calibration check
# -----------------------------------------------------------------------------
null_check <- function(a = 6, k = 6, n_sims = 20000,
                       sigma_donor = 1, sigma_res = 1) {
  r <- twoway_anova_fast(a = a, k = k, beta = 0,
                         sigma_donor = sigma_donor,
                         sigma_int = 0,
                         sigma_res = sigma_res,
                         n_sims = n_sims)
  cat(sprintf("Null calibration (a=%d, k=%d, n_sims=%d)\n", a, k, n_sims))
  cat(sprintf("Type I (condition): %.4f  (expect ~0.05)\n",
              mean(r$p_cond        < 0.05)))
  cat(sprintf("Type I (interaction): %.4f  (expect ~0.05)\n",
              mean(r$p_interaction < 0.05)))
  cat(sprintf("Type I (donor): %.4f  (expect ~0.05)\n",
              mean(r$p_donor       < 0.05)))
  invisible(r)
}

cat("=== Null calibration ===\n")
null_check(sigma_donor = 0)
# -----------------------------------------------------------------------------
# Main parameter sweep
# -----------------------------------------------------------------------------
param_grid <- expand.grid(
  a = c(3, 4, 5, 6, 8, 12, 20, 50),
  k = 6,
  beta = c(0.5, 1.0),
  sigma_int = c(0, 0.3, 0.7, 1.0),
  sigma_donor = 1.0,
  sigma_res = 1.0
)
n_sims <- 5000

cat("=== Running", nrow(param_grid), "parameter cells ===\n")
rows <- vector("list", nrow(param_grid))

for (i in seq_len(nrow(param_grid))) {
  p <- param_grid[i, ]
  cat(sprintf("[%2d/%2d] a=%2d  k=%d  beta=%.2f  sigma_int=%.2f\n",
              i, nrow(param_grid), p$a, p$k, p$beta, p$sigma_int))
  
  r <- twoway_anova_fast(a = p$a, k = p$k, beta = p$beta,
                         sigma_donor = p$sigma_donor,
                         sigma_int = p$sigma_int,
                         sigma_res = p$sigma_res,
                         n_sims = n_sims)
  
  rows[[i]] <- data.frame(
    a = p$a,
    k = p$k,
    beta = p$beta,
    sigma_int = p$sigma_int,
    power_cond = mean(r$p_cond < 0.05),
    power_interact = mean(r$p_interaction < 0.05),
    power_donor = mean(r$p_donor < 0.05),
    median_sigma_int = median(r$sigma_int_hat),
    mean_sigma_int = mean(r$sigma_int_hat),
    pct_zero_int = mean(r$sigma_int_hat < 1e-8),
    pct_all_same = mean(r$all_same_sign)
  )
}
summary_tbl <- tibble::as_tibble(bind_rows(rows))
print(summary_tbl, n = Inf)

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
theme_set(theme_minimal(base_size = 11))

fig_main <- ggplot(summary_tbl,
                   aes(x = a, y = power_cond,
                       colour = factor(sigma_int),
                       linetype = factor(beta))) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 0.8, linetype = "dashed", alpha = 0.4) +
  scale_x_continuous(breaks = unique(summary_tbl$a)) +
  labs(title = "Power to detect the mean treatment effect",
       x = "Number of donors (a)", y = "Power",
       colour = expression(sigma[int]), linetype = expression(beta))

tiff("../pics/figI1.tiff", units = "in", res = 300, width = 6, height = 5)
print(fig_main)
dev.off()

fig_het <- summary_tbl %>%
  filter(sigma_int > 0,
         beta == 1) %>%
  ggplot(aes(x = a, y = power_interact,
             colour = factor(sigma_int),
             )) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 0.8, linetype = "dashed", alpha = 0.4) +
  scale_x_continuous(breaks = unique(summary_tbl$a)) +
  labs(title = "Power to detect donor-by-treatment heterogeneity",
       x = "Number of donors (a)", y = "Power (F)",
       colour = expression(sigma[int]), linetype = expression(beta))

tiff("../pics/figI2.tiff", units = "in", res = 300, width = 6, height = 5)
print(fig_het)
dev.off()


fig_zero <- summary_tbl %>%
  filter(sigma_int > 0,
         beta == 1) %>%
  ggplot(aes(x = a, y = pct_zero_int,
             colour = factor(sigma_int))) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = unique(summary_tbl$a)) +
  labs(title = "Rate at which donor x treatment effect is estimated at zero",
       x = "Number of donors (a)", y = "P(σ interaction=0)",
       colour = expression(sigma[int]), linetype = expression(beta))


tiff("../pics/figI3.tiff", units = "in", res = 300, width = 6, height = 5)
print(fig_zero)
dev.off()

fig_sign <- summary_tbl %>%
  filter(sigma_int > 0) %>%
  ggplot(aes(x = a, y = pct_all_same,
             colour = factor(sigma_int),
             linetype = factor(beta))) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = unique(summary_tbl$a)) +
  labs(title = "Probability all sampled donors respond in the same direction",
       x = "Number of donors (a)", y = "P(all a donors same sign)",
       colour = expression(sigma[int]), linetype = expression(beta))
tiff("../pics/figI4.tiff", units = "in", res = 300, width = 6, height = 5)
print(fig_sign)
dev.off()


ggsave("fig_I1_power_main_effect.png",  fig_main, width = 7.5, height = 4.8, dpi = 160)
ggsave("fig_I2_power_heterogeneity.png", fig_het,  width = 7.5, height = 4.8, dpi = 160)
ggsave("fig_I3_sigma_int_zero_rate.png", fig_zero, width = 7.5, height = 4.8, dpi = 160)
ggsave("fig_I4_direction_concordance.png", fig_sign, width = 7.5, height = 4.8, dpi = 160)



