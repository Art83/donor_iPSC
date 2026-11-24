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
# Direction of response (unanumity). The true treatment effect in donor d is
#   beta + v_{d,treated} - v_{d,control}  ~  N(beta, 2 * sigma_int^2)
# so the share of donors in the population who truly respond in the opposite
# direction to the mean is
#   q = pnorm(-beta / (sqrt(2) * sigma_int))
#
# simulator separates three events:
#   observed unanimity   all a observed per-donor effects share a sign
#   true unanimity       all a sampled donors truly share a sign
#   masked discordance   observed unanimity although a sampled donor truly
#                        responds the other way (replicate noise hides it)
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
  
  # The same diagnostic on the true per-donor effects of the sampled donors
  true_effect <- beta + v[, 2L, , drop = TRUE] - v[, 1L, , drop = TRUE]
  if (a == 1L) true_effect <- matrix(true_effect, nrow = 1L)
  n_pos_true <- colSums(true_effect > 0)
  all_same_true <- (n_pos_true == 0L) | (n_pos_true == a)
  
  list(
    all_same_true = all_same_true,
    masked_discordance = all_same_sign & !all_same_true,
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
    pct_all_same = mean(r$all_same_sign),
    # population share of donors truly responding against the mean effect
    frac_opposite = pnorm(-p$beta / (sqrt(2) * p$sigma_int)),
    pct_all_same_true = mean(r$all_same_true),
    pct_masked = mean(r$masked_discordance),
    # closed form for observed unanimity, used as a check on the simulation
    pct_all_same_expected = {
      q_obs <- pnorm(-p$beta / sqrt(2 * p$sigma_int^2 + 2 * p$sigma_res^2 / p$k))
      (1 - q_obs)^p$a + q_obs^p$a
    }
  )
}
summary_tbl <- tibble::as_tibble(bind_rows(rows))
print(summary_tbl, n = Inf, width = Inf)
write.csv(summary_tbl, "tables/unanimity_summary.csv")

# -----------------------------------------------------------------------------
# Direction of response: unanimity against the true share of opposite responders
# -----------------------------------------------------------------------------
cat("\n=== Unanimous panels and the true share of opposite responders ===\n")
cat(sprintf("Largest gap between simulated and closed-form unanimity: %.1f percentage points\n\n",
            100 * max(abs(summary_tbl$pct_all_same - summary_tbl$pct_all_same_expected))))

sign_tbl <- summary_tbl %>%
  filter(sigma_int > 0) %>%
  arrange(frac_opposite, a)

cat(sprintf("  %-5s %-9s %-22s %s\n", "beta", "sigma_int",
            "donors truly opposite", "P(all donors same sign), observed, at a = 3, 5, 10*, 20"))
for (key in unique(paste(sign_tbl$beta, sign_tbl$sigma_int))) {
  s <- sign_tbl[paste(sign_tbl$beta, sign_tbl$sigma_int) == key, ]
  q_obs <- pnorm(-s$beta[1] / sqrt(2 * s$sigma_int[1]^2 + 2 / s$k[1]))
  cat(sprintf("  %-5.1f %-9.1f %20.1f%%   %5.1f%%  %5.1f%%  %5.1f%%  %5.1f%%\n",
              s$beta[1], s$sigma_int[1], 100 * s$frac_opposite[1],
              100 * s$pct_all_same[s$a == 3], 100 * s$pct_all_same[s$a == 5],
              100 * ((1 - q_obs)^10 + q_obs^10),
              100 * s$pct_all_same[s$a == 20]))
}
cat("  * a = 10 is not in the simulation grid and is given in closed form\n")

cat("\nObserved unanimity split by what the sampled donors truly did (a = 3 and 5)\n")
# masked discordance: looks unanimous, a sampled donor truly responds the other way
# false discordance:  looks discordant, all sampled donors truly respond the same way
cat(sprintf("  %-5s %-9s %-3s %10s %18s %20s %19s\n", "beta", "sigma_int", "a",
            "observed", "truly unanimous", "masked discordance", "false discordance"))
for (i in which(sign_tbl$a %in% c(3, 5))) {
  s <- sign_tbl[i, ]
  cat(sprintf("  %-5.1f %-9.1f %-3d %9.1f%% %17.1f%% %19.1f%% %18.1f%%\n",
              s$beta, s$sigma_int, s$a, 100 * s$pct_all_same,
              100 * s$pct_all_same_true, 100 * s$pct_masked,
              100 * (s$pct_all_same_true - (s$pct_all_same - s$pct_masked))))
}

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


# Each condition is labelled by the share of donors who truly respond in the
# opposite direction, which is what makes a unanimous panel informative or not.
sign_fig_data <- summary_tbl %>%
  filter(sigma_int > 0) %>%
  mutate(opp_label = sprintf("%.0f%%", 100 * frac_opposite)) %>%
  arrange(frac_opposite)
sign_fig_data$opp_label <- factor(sign_fig_data$opp_label,
                                  levels = unique(sign_fig_data$opp_label))

fig_sign_new <- ggplot(sign_fig_data,
                   aes(x = a, y = pct_all_same, colour = opp_label, group = opp_label)) +
  annotate("rect", xmin = 3, xmax = 5, ymin = 0, ymax = 1,
           fill = "grey50", alpha = 0.12) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  scale_x_continuous(trans = "log10", breaks = unique(summary_tbl$a)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_colour_manual(values = c("#FDD0A2", "#FDAE6B", "#FD8D3C",
                                 "#E6550D", "#A63603", "#5C1F02")) +
  labs(title = "Probability all sampled donors respond in the same direction",
       x = "Number of donors (log scale)", y = "P(all donors have same sign)",
       colour = "Donors truly responding\nin the opposite direction") +
  theme(legend.position = "right")
tiff("pics/figI4_new.tiff", units = "in", res = 300, width = 6, height = 5)
print(fig_sign_new)
dev.off()




