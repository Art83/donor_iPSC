################################################################################
#
# Disease contrast with donors nested in group
#
# Scripts 01-03 ask whether a study can detect or estimate donor-level variance.
# This script asks:
# with n patient donors against n control donors, each measured k times, can the
# study detect a disease effect, and what does it report when it does?
#
# Model:
#   y_{gdi} = delta * [g = patient] + u_{d(g)} + eps_{gdi}
#     u_d ~ N(0, sigma_donor^2)    donor effect, nested in group
#     eps ~ N(0, sigma_res^2)      replicate (organoid / well) noise
#     ICC = sigma_donor^2 / (sigma_donor^2 + sigma_res^2)
#
# Scale. Total variance is held at 1 (sigma_donor^2 = ICC, sigma_res^2 = 1 - ICC)
# and delta is the disease effect in total-SD units
#
# The donor is the unit of replication. The test is a two-sample t-test on the
# 2n donor means, with 2n - 2 degrees of freedom:
#   Var(donor mean) = sigma_donor^2 + sigma_res^2 / k
#   Var(group difference) = (2 / n) * (sigma_donor^2 + sigma_res^2 / k)
#
# Link to the design effect in script 02 (Kish, 1965):
#   Var(group difference) = 2 * sigma_total^2 * DEFF / (n * k),
#   DEFF = 1 + (k - 1) * ICC
# so the effective number of independent observations per group is
#   n * k / DEFF, which tends to n / ICC as k grows. Replicates cannot buy
#   power past that ceiling; only donors can.
#
# Link to scripts 01-02 (ICC estimation): the denominator of this t-test is the
# between-donor mean square, the same quantity whose imprecision at small n is
# shown there. Part C below reports that imprecision directly as the spread of
# the estimated standard error around its true value.
#
# Part A: exact power from the noncentral t distribution (no simulation)
# Part B: minimum detectable effect and donors required
# Part C: precision of the error term
# Part D: replicates versus donors (the ICC ceiling)
# Part E: Monte Carlo check of Part A, plus exaggeration of significant effects
# Part F: false positive rate when replicates are analysed as independent
# Part G: the same results by fold change, with exaggeration of significant effects
#
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

set.seed(42)

# Same grids as script 01
TRUE_ICCS <- c(0.05, 0.10, 0.20, 0.30, 0.50)
N_PER_GROUP <- c(2, 3, 5, 10, 20, 50)
K_REPS <- 6
TOTAL_SD <- 1
ALPHA <- 0.05

# Disease effect sizes, in total-SD units
EFFECTS <- c(0.5, 1, 1.5, 2, 3)

# Monte Carlo iterations (Parts E and F)
N_SIMS <- 20000

palette_n <- c(
  "n=2" = "#E63946", "n=3" = "#F4A261", "n=5" = "#2A9D8F",
  "n=10" = "#264653", "n=20" = "#6A4C93", "n=50" = "#1982C4"
)

palette_icc <- c(
  "5%" = "#BFDBFE", "10%" = "#93C5FD", "20%" = "#3B82F6",
  "30%" = "#1D4ED8", "50%" = "#1E3A5F"
)

# One hue, light to dark, for the ordered effect sizes
palette_effect <- c(
  "0.5" = "#FDAE6B", "1" = "#FD8D3C", "1.5" = "#E6550D",
  "2" = "#A63603", "3" = "#5C1F02"
)

theme_pub <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    plot.subtitle = element_text(size = 11, color = "grey40", hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 12),
    plot.margin = margin(15, 15, 15, 15)
  )


################################################################################
# CORE FUNCTIONS
################################################################################

# Variance components with total variance fixed at TOTAL_SD^2
donor_sd_from_icc <- function(icc) sqrt(icc) * TOTAL_SD
resid_sd_from_icc <- function(icc) sqrt(1 - icc) * TOTAL_SD

# Standard error of the group difference, donor as the unit
se_difference <- function(n, k, icc) {
  sqrt((2 / n) * (donor_sd_from_icc(icc)^2 + resid_sd_from_icc(icc)^2 / k))
}

# Exact power of the two-sided two-sample t-test on donor means
power_disease <- function(n, k, icc, delta, alpha = ALPHA) {
  df <- 2 * n - 2
  ncp <- delta / se_difference(n, k, icc)
  tcrit <- qt(1 - alpha / 2, df)
  pt(tcrit, df, ncp = ncp, lower.tail = FALSE) + pt(-tcrit, df, ncp = ncp)
}

# Smallest effect detected with the target power
min_detectable_effect <- function(n, k, icc, target = 0.80) {
  uniroot(function(d) power_disease(n, k, icc, d) - target,
          interval = c(1e-6, 200))$root
}

# Smallest number of donors per group reaching the target power
donors_required <- function(k, icc, delta, target = 0.80, n_max = 5000) {
  n <- 2
  while (n <= n_max && power_disease(n, k, icc, delta) < target) n <- n + 1
  if (n > n_max) NA_integer_ else as.integer(n)
}


################################################################################
# PART A: EXACT POWER
################################################################################

cat("=== PART A: Exact power for the disease effect ===\n\n")

power_grid <- expand.grid(icc = TRUE_ICCS, n = N_PER_GROUP, delta = EFFECTS) %>%
  mutate(
    power = mapply(power_disease, n = n, k = K_REPS, icc = icc, delta = delta),
    # the same effect expressed against the residual (within-donor) SD
    delta_resid_sd = delta / sqrt(1 - icc),
    icc_label = factor(paste0("ICC = ", icc * 100, "%"),
                       levels = paste0("ICC = ", TRUE_ICCS * 100, "%")),
    effect_label = factor(as.character(delta), levels = as.character(EFFECTS))
  )

cat(sprintf("Power (%%), k = %d replicates per donor, alpha = %.2f, donor as unit\n",
            K_REPS, ALPHA))
cat("Effect in total-SD units\n\n")
for (ic in TRUE_ICCS) {
  cat(sprintf("ICC = %.0f%%\n", ic * 100))
  cat(sprintf("  %-12s %s\n", "donors/group",
              paste(sprintf("d=%-4s", EFFECTS), collapse = " ")))
  for (nn in N_PER_GROUP) {
    pw <- power_grid$power[power_grid$icc == ic & power_grid$n == nn]
    cat(sprintf("  %-12d %s\n", nn,
                paste(sprintf("%5.1f ", pw * 100), collapse = " ")))
  }
  cat("\n")
}

pD1 <- ggplot(power_grid, aes(x = n, y = power, color = effect_label,
                              group = effect_label)) +
  annotate("rect", xmin = 1, xmax = 5, ymin = 0, ymax = 1,
           fill = "grey50", alpha = 0.12) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  facet_wrap(~icc_label, nrow = 1) +
  scale_x_continuous(trans = "log10", breaks = N_PER_GROUP) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_color_manual(values = palette_effect) +
  labs(
    title = "Power to detect a disease effect",
    x = "Donors per group (log scale)",
    y = "Power",
    color = "Disease effect (SDs)"
  ) +
  theme_pub

tiff("pics/figD_power_disease_effect.tiff", units = "in", res = 300, width = 8, height = 5, compression = "lzw")
print(pD1)
dev.off()

################################################################################
# PART B: MINIMUM DETECTABLE EFFECT AND DONORS REQUIRED
################################################################################

cat("=== PART B: Minimum detectable effect at 80% power ===\n\n")

mde_grid <- expand.grid(icc = TRUE_ICCS, n = c(2, 3, 4, 5, 7, 10, 15, 20, 30, 50)) %>%
  mutate(
    mde = mapply(min_detectable_effect, n = n, k = K_REPS, icc = icc),
    mde_resid_sd = mde / sqrt(1 - icc),
    icc_label = factor(paste0(icc * 100, "%"), levels = paste0(TRUE_ICCS * 100, "%"))
  )

cat("Smallest effect detected with 80% power (total SDs; residual SDs in brackets)\n\n")
cat(sprintf("  %-12s %s\n", "donors/group",
            paste(sprintf("ICC=%-9s", paste0(TRUE_ICCS * 100, "%")), collapse = " ")))
for (nn in N_PER_GROUP) {
  r <- mde_grid[mde_grid$n == nn, ]
  cat(sprintf("  %-12d %s\n", nn,
              paste(sprintf("%4.2f (%4.2f) ", r$mde, r$mde_resid_sd), collapse = " ")))
}

cat("\nDonors per group required for 80% power\n\n")
req_grid <- expand.grid(icc = TRUE_ICCS, delta = EFFECTS) %>%
  mutate(n_required = mapply(donors_required, k = K_REPS, icc = icc, delta = delta))
cat(sprintf("  %-10s %s\n", "ICC",
            paste(sprintf("d=%-5s", EFFECTS), collapse = " ")))
for (ic in TRUE_ICCS) {
  r <- req_grid[req_grid$icc == ic, ]
  cat(sprintf("  %-10s %s\n", paste0(ic * 100, "%"),
              paste(sprintf("%-7d", r$n_required), collapse = " ")))
}

pD2 <- ggplot(mde_grid, aes(x = n, y = mde, color = icc_label, group = icc_label)) +
  annotate("rect", xmin = 1, xmax = 5, ymin = 0, ymax = Inf,
           fill = "grey50", alpha = 0.12) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  scale_x_continuous(trans = "log10", breaks = N_PER_GROUP) +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, 10, 1)) +
  scale_color_manual(values = palette_icc) +
  labs(
    title = "Smallest disease effect detected with 80% power",
    x = "Donors per group (log scale)",
    y = "Minimum detectable effect (SDs)",
    color = "Donor ICC"
  ) +
  theme_pub

tiff("pics/figD_smallest_effect.tiff", units = "in", res = 300, width = 7, height = 5, compression = "lzw")
print(pD2)
dev.off()

################################################################################
# PART C: PRECISION OF THE ERROR TERM
#
# The t-test divides the group difference by a standard error estimated from
# the between-donor mean square on 2n - 2 degrees of freedom. The ratio of the
# estimated to the true standard error is sqrt(chi^2_df / df).
################################################################################

cat("\n=== PART C: Precision of the error term ===\n\n")

se_precision <- data.frame(n = N_PER_GROUP) %>%
  mutate(
    df = 2 * n - 2,
    t_crit = qt(1 - ALPHA / 2, df),
    se_ratio_lo = sqrt(qchisq(0.025, df) / df),
    se_ratio_hi = sqrt(qchisq(0.975, df) / df)
  )

cat("Estimated SE as a multiple of the true SE (central 95%), and the critical t\n\n")
cat(sprintf("  %-12s %4s  %8s  %s\n", "donors/group", "df", "t crit", "estimated/true SE"))
for (i in seq_len(nrow(se_precision))) {
  r <- se_precision[i, ]
  cat(sprintf("  %-12d %4d  %8.2f  %.2f to %.2f\n",
              r$n, r$df, r$t_crit, r$se_ratio_lo, r$se_ratio_hi))
}
cat("\n  1 donor per group: df = 0. The between-donor variance cannot be estimated\n")
cat("  and no test of the disease effect at the donor level exists.\n")


################################################################################
# PART D: REPLICATES VERSUS DONORS
################################################################################

cat("\n=== PART D: Replicates versus donors ===\n\n")

K_GRID <- c(1, 2, 3, 6, 12, 25, 50, 100)
DELTA_D <- 1

rep_grid <- expand.grid(icc = c(0.10, 0.20, 0.50), n = c(3, 5, 10, 20), k = K_GRID) %>%
  mutate(
    power = mapply(power_disease, n = n, k = k, icc = icc, delta = DELTA_D),
    power_ceiling = mapply(function(n, icc) power_disease(n, 1e9, icc, DELTA_D), n, icc),
    deff = 1 + (k - 1) * icc,
    n_eff_per_group = n * k / deff,
    n_eff_ceiling = n / icc,
    icc_label = factor(paste0("ICC = ", icc * 100, "%"),
                       levels = paste0("ICC = ", c(10, 20, 50), "%")),
    n_label = factor(paste0("n=", n), levels = paste0("n=", c(3, 5, 10, 20)))
  )

cat(sprintf("Power (%%) for an effect of %g SD as replicates increase\n\n", DELTA_D))
for (ic in c(0.10, 0.20, 0.50)) {
  cat(sprintf("ICC = %.0f%%\n", ic * 100))
  cat(sprintf("  %-12s %s  %s\n", "donors/group",
              paste(sprintf("k=%-4d", K_GRID), collapse = " "), "ceiling"))
  for (nn in c(3, 5, 10, 20)) {
    r <- rep_grid[rep_grid$icc == ic & rep_grid$n == nn, ]
    cat(sprintf("  %-12d %s  %5.1f\n", nn,
                paste(sprintf("%5.1f ", r$power * 100), collapse = " "),
                r$power_ceiling[1] * 100))
  }
  cat("\n")
}

cat("Effective independent observations per donor = k / DEFF (ceiling 1 / ICC)\n\n")
for (ic in TRUE_ICCS) {
  cat(sprintf("  ICC = %2.0f%%: k = 6 gives %.2f, k = 100 gives %.2f, ceiling %.1f\n",
              ic * 100, 6 / (1 + 5 * ic), 100 / (1 + 99 * ic), 1 / ic))
}

# Numerical check of the identity Var(diff) = 2 * sigma_total^2 * DEFF / (n * k)
chk <- with(list(n = 3, k = 6, icc = 0.2), {
  lhs <- se_difference(n, k, icc)^2
  rhs <- 2 * TOTAL_SD^2 * (1 + (k - 1) * icc) / (n * k)
  c(lhs = lhs, rhs = rhs)
})
stopifnot(isTRUE(all.equal(chk[["lhs"]], chk[["rhs"]])))

pD3 <- ggplot(rep_grid, aes(x = k, y = power, color = n_label, group = n_label)) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  facet_wrap(~icc_label, nrow = 1) +
  scale_x_continuous(trans = "log10", breaks = c(1, 3, 6, 12, 25, 50, 100)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_color_manual(values = palette_n) +
  labs(
    title = "Adding replicates does not substitute for donors",
    subtitle = paste0("Power for a disease effect of ", DELTA_D, " SD"),
    x = "Replicates per donor (log scale)",
    y = "Power",
    color = "Donors per group"
  ) +
  theme_pub

tiff("pics/figD_adding_replicates.tiff", units = "in", res = 300, width = 7, height = 5, compression = "lzw")
print(pD3)
dev.off()

################################################################################
# PART E: MONTE CARLO CHECK AND EXAGGERATION OF SIGNIFICANT EFFECTS
#
# Simulates donor means directly. Among studies reaching p < 0.05, the
# exaggeration ratio (Gelman & Carlin, 2014) is |estimated effect| / true
# effect, and the sign error rate is the share pointing the wrong way.
################################################################################

cat("\n=== PART E: Monte Carlo check and exaggeration ===\n\n")

simulate_contrast <- function(n, k, icc, delta, n_sims = N_SIMS) {
  sd_mean <- sqrt(donor_sd_from_icc(icc)^2 + resid_sd_from_icc(icc)^2 / k)
  ctrl <- matrix(rnorm(n * n_sims, 0, sd_mean), n_sims, n)
  pat <- matrix(rnorm(n * n_sims, delta, sd_mean), n_sims, n)
  est <- rowMeans(pat) - rowMeans(ctrl)
  v <- (apply(ctrl, 1, var) + apply(pat, 1, var)) / 2
  tval <- est / sqrt(2 * v / n)
  pval <- 2 * pt(abs(tval), 2 * n - 2, lower.tail = FALSE)
  sig <- pval < ALPHA
  data.frame(
    icc = icc, n = n, delta = delta,
    power_mc = mean(sig),
    power_exact = power_disease(n, k, icc, delta),
    exaggeration = if (any(sig)) median(abs(est[sig])) / delta else NA_real_,
    sign_error = if (any(sig)) mean(est[sig] < 0) else NA_real_
  )
}

mc_grid <- expand.grid(icc = c(0.20, 0.50), n = N_PER_GROUP, delta = c(0.5, 1, 2))
mc_results <- do.call(rbind, lapply(seq_len(nrow(mc_grid)), function(i) {
  simulate_contrast(mc_grid$n[i], K_REPS, mc_grid$icc[i], mc_grid$delta[i])
}))

cat(sprintf("Largest gap between simulated and exact power: %.2f percentage points\n\n",
            max(abs(mc_results$power_mc - mc_results$power_exact)) * 100))

cat("Among significant results: median exaggeration of the effect, and sign errors\n\n")
cat(sprintf("  %-6s %-8s %-12s %8s %14s %12s\n",
            "ICC", "effect", "donors/group", "power", "exaggeration", "wrong sign"))
for (i in seq_len(nrow(mc_results))) {
  r <- mc_results[i, ]
  cat(sprintf("  %-6s %-8g %-12d %7.1f%% %13.2fx %11.1f%%\n",
              paste0(r$icc * 100, "%"), r$delta, r$n, r$power_mc * 100,
              r$exaggeration, r$sign_error * 100))
}


################################################################################
# PART F: FALSE POSITIVES WHEN REPLICATES ARE ANALYSED AS INDEPENDENT
#
# No disease effect. A pooled t-test on all n * k replicates per group treats
# organoids as independent. With one donor per group this is the only test
# available, and every donor difference is read as a disease effect.
################################################################################

cat("\n=== PART F: False positives, replicates analysed as independent ===\n\n")

naive_fpr <- function(n, k, icc, n_sims = N_SIMS) {
  dsd <- donor_sd_from_icc(icc)
  gen <- function() {
    d <- matrix(rnorm(n * n_sims, 0, dsd), n_sims, n)
    d[, rep(seq_len(n), each = k), drop = FALSE] +
      matrix(rnorm(n * k * n_sims, 0, resid_sd_from_icc(icc)), n_sims, n * k)
  }
  a <- gen(); b <- gen()
  m  <- n * k
  sp <- (apply(a, 1, var) + apply(b, 1, var)) / 2
  tv <- (rowMeans(a) - rowMeans(b)) / sqrt(2 * sp / m)
  mean(2 * pt(abs(tv), 2 * m - 2, lower.tail = FALSE) < ALPHA)
}

fpr_grid <- expand.grid(icc = TRUE_ICCS, n = c(1, 2, 3, 5, 10, 20)) %>%
  mutate(
    fpr = mapply(naive_fpr, n = n, k = K_REPS, icc = icc),
    icc_label = factor(paste0(icc * 100, "%"), levels = paste0(TRUE_ICCS * 100, "%"))
  )

cat(sprintf("False positive rate (%%), k = %d, nominal 5%%\n\n", K_REPS))
cat(sprintf("  %-12s %s\n", "donors/group",
            paste(sprintf("ICC=%-4s", paste0(TRUE_ICCS * 100, "%")), collapse = " ")))
for (nn in c(1, 2, 3, 5, 10, 20)) {
  r <- fpr_grid[fpr_grid$n == nn, ]
  cat(sprintf("  %-12d %s\n", nn,
              paste(sprintf("%6.1f   ", r$fpr * 100), collapse = " ")))
}

pD4 <- ggplot(fpr_grid, aes(x = n, y = fpr, color = icc_label, group = icc_label)) +
  geom_hline(yintercept = ALPHA, linetype = "dashed", color = "grey50") +
  annotate("text", x = 14, y = ALPHA + 0.02, label = "nominal 5%",
           color = "grey40", size = 3.5) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  scale_x_continuous(trans = "log10", breaks = c(1, 2, 3, 5, 10, 20)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  scale_color_manual(values = palette_icc) +
  labs(
    title = "False positives when replicates are treated as independent",
    x = "Donors per group (log scale)",
    y = "False positive rate",
    color = "Donor ICC"
  ) +
  theme_pub

tiff("pics/figD_FP.tiff", units = "in", res = 300, width = 7, height = 5, compression = "lzw")
print(pD4)
dev.off()

################################################################################
# PART G: FOLD-CHANGE ANCHORING
#
# Parts A-F express the disease effect in SD units. The two scales are linked
# by how much lines differ from each other:
#   effect in SD = log2(fold change) / SD of log2 abundance across lines
#
################################################################################

cat("\n=== PART G: Fold-change anchoring ===\n\n")

specs <- read.csv("tables/empirical_specs.csv", stringsAsFactors = FALSE)
sd_rows <- specs %>%
  filter(parameter == "line_sd_log2",
         dataset == "Mirauta 2020" |
           cell_type %in% c("macrophages, naive", "cardiomyocytes, untreated", "sensory neurons")) %>%
  arrange(value)
LINE_SDS <- setNames(sd_rows$value,
                     sprintf("%s, %s (SD %.2f)", sub(",.*$", "", sd_rows$cell_type),
                             ifelse(sd_rows$assay == "proteomics", "protein", "RNA"),
                             sd_rows$value))
FOLD_CHANGES <- c(1.05, 1.1, 1.2, 1.5)
ICCS_G <- c(0.20, 0.50)
N_GRID_G <- c(2, 3, 5, 10, 20, 50, 100)

fc_to_sd <- function(fc, line_sd) log2(fc) / line_sd
sd_to_fc <- function(d, line_sd) 2^(d * line_sd)
sd_label <- function(line_sd) factor(names(LINE_SDS)[match(line_sd, LINE_SDS)],
                                     levels = names(LINE_SDS))

palette_fc <- c("1.05" = "#FDAE6B", "1.1" = "#FD8D3C",
                "1.2" = "#D94801", "1.5" = "#5C1F02")

fc_grid <- expand.grid(icc = ICCS_G, n = N_GRID_G, fc = FOLD_CHANGES,
                       line_sd = unname(LINE_SDS)) %>%
  mutate(
    delta = fc_to_sd(fc, line_sd),
    power = mapply(power_disease, n = n, k = K_REPS, icc = icc, delta = delta),
    cells = sd_label(line_sd),
    icc_label = factor(paste0("ICC = ", icc * 100, "%"),
                       levels = paste0("ICC = ", ICCS_G * 100, "%")),
    fc_label = factor(as.character(fc), levels = as.character(FOLD_CHANGES))
  )

for (s in LINE_SDS) {
  cat(sprintf("%s: median SD across lines = %.3f log2 units, k = %d\n\n",
              as.character(sd_label(s)), s, K_REPS))
  cat(sprintf("  %-6s %-8s %-6s %s  %s\n", "FC", "SD", "ICC",
              paste(sprintf("n=%-4d", c(3, 5, 10, 20)), collapse = " "),
              "donors/group for 80%"))
  for (f in FOLD_CHANGES) {
    for (ic in ICCS_G) {
      r <- fc_grid[fc_grid$line_sd == s & fc_grid$fc == f & fc_grid$icc == ic &
                     fc_grid$n %in% c(3, 5, 10, 20), ]
      cat(sprintf("  %-6g %-8.2f %-6s %s  %d\n", f, fc_to_sd(f, s), paste0(ic * 100, "%"),
                  paste(sprintf("%5.1f ", r$power * 100), collapse = " "),
                  donors_required(K_REPS, ic, fc_to_sd(f, s))))
    }
  }
  cat("\n")
}

cat("\nSmallest fold change detected with 80% power\n\n")
fc_threshold <- expand.grid(icc = ICCS_G, n = c(3, 5, 10, 20),
                            line_sd = unname(LINE_SDS)) %>%
  mutate(mde = mapply(min_detectable_effect, n = n, k = K_REPS, icc = icc),
         fc_min = sd_to_fc(mde, line_sd),
         cells = sd_label(line_sd))
for (s in LINE_SDS) {
  cat(sprintf("%s\n", as.character(sd_label(s))))
  for (nn in c(3, 5, 10, 20)) {
    r <- fc_threshold[fc_threshold$line_sd == s & fc_threshold$n == nn, ]
    cat(sprintf("  %2d donors/group: %s\n", nn,
                paste(sprintf("FC %.2f (ICC %.0f%%)", r$fc_min, r$icc * 100), collapse = ", ")))
  }
}

# Exaggeration of significant results at these fold changes
fc_mc_grid <- expand.grid(icc = ICCS_G, n = c(2, 3, 5, 10, 20, 50), fc = FOLD_CHANGES,
                          line_sd = unname(LINE_SDS))
fc_exag <- do.call(rbind, lapply(seq_len(nrow(fc_mc_grid)), function(i) {
  r <- simulate_contrast(fc_mc_grid$n[i], K_REPS, fc_mc_grid$icc[i],
                         fc_to_sd(fc_mc_grid$fc[i], fc_mc_grid$line_sd[i]))
  r$fc <- fc_mc_grid$fc[i]
  r$line_sd <- fc_mc_grid$line_sd[i]
  r
})) %>%
  mutate(
    cells = sd_label(line_sd),
    # the fold change a significant study reports, at the median exaggeration
    fc_reported = 2^(exaggeration * log2(fc)),
    icc_label = factor(paste0("ICC = ", icc * 100, "%"),
                       levels = paste0("ICC = ", ICCS_G * 100, "%")),
    fc_label = factor(as.character(fc), levels = as.character(FOLD_CHANGES))
  )

cat("\nAmong significant results: exaggeration and the fold change reported\n\n")
cat(sprintf("  %-28s %-6s %-6s %-12s %8s %14s %14s\n",
            "cells", "FC", "ICC", "donors/group", "power", "exaggeration", "reported FC"))
for (i in seq_len(nrow(fc_exag))) {
  r <- fc_exag[i, ]
  if (r$n %in% c(3, 5, 10, 20)) {
    cat(sprintf("  %-28s %-6g %-6s %-12d %7.1f%% %13.2fx %14.2f\n",
                as.character(r$cells), r$fc, paste0(r$icc * 100, "%"), r$n, r$power_mc * 100,
                r$exaggeration, r$fc_reported))
  }
}

pD5 <- ggplot(fc_grid, aes(x = n, y = power, color = fc_label, group = fc_label)) +
  annotate("rect", xmin = 1, xmax = 5, ymin = 0, ymax = 1,
           fill = "grey50", alpha = 0.12) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  facet_grid(cells ~ icc_label, labeller = labeller(cells = label_wrap_gen(18))) +
  scale_x_continuous(trans = "log10", breaks = N_GRID_G) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_color_manual(values = palette_fc) +
  labs(
    title = "Power to detect a disease effect of a given fold change",
    x = "Donors per group (log scale)",
    y = "Power",
    color = "Fold change"
  ) +
  theme_pub

pD6 <- ggplot(fc_exag, aes(x = n, y = exaggeration, color = fc_label, group = fc_label)) +
  annotate("rect", xmin = 1, xmax = 5, ymin = -Inf, ymax = Inf,
           fill = "grey50", alpha = 0.12) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  facet_grid(cells ~ icc_label, labeller = labeller(cells = label_wrap_gen(18))) +
  scale_x_continuous(trans = "log10", breaks = c(2, 3, 5, 10, 20, 50)) +
  scale_y_continuous(limits = c(0.9, NA), breaks = seq(1, 10, 1)) +
  scale_color_manual(values = palette_fc) +
  labs(
    title = "Significant results overstate small disease effects",
    x = "Donors per group (log scale)",
    y = "Exaggeration ratio",
    color = "True fold change"
  ) +
  theme_pub


tiff("pics/figD_power_FC.tiff", units = "in", res = 300, width = 8, height = 11, compression = "lzw")
print(pD5)
dev.off()

tiff("pics/figD_exaggeration.tiff", units = "in", res = 300, width = 8, height = 11, compression = "lzw")
print(pD6)
dev.off()


################################################################################
# PART H: MEASURED AD EFFECT SIZES
################################################################################

cat("\n=== PART H: Measured AD effect sizes ===\n\n")

effects_tbl <- read.csv("tables/effects_summary.csv", stringsAsFactors = FALSE)
variance_tbl <- read.csv("tables/variance_specs.csv", stringsAsFactors = FALSE)
ad <- effects_tbl[effects_tbl$contrast == "human: AD vs control", ]

MEASURED <- data.frame(
  effect = c("AD, median protein", "AD, top tenth of proteins"),
  delta = sqrt(pmax(c(ad$median_abs_d^2 - ad$null_median^2, ad$q90_abs_d^2 - ad$null_q90^2), 0)),
  fc = c(ad$median_fold_change, ad$q90_fold_change),
  stringsAsFactors = FALSE
)
ICC_ORGANOID <- variance_tbl$value[variance_tbl$parameter == "line_icc" &
                                     variance_tbl$condition == "within diagnosis"]
ICCS_H <- c(ICCS_G, ICC_ORGANOID)
ICC_LABELS_H <- c(paste0("ICC = ", ICCS_G * 100, "%"),
                  sprintf("ICC = %.0f%% (organoids)", ICC_ORGANOID * 100))
N_GRID_H <- c(2, 3, 5, 10, 20, 50, 100, 200)

measured_grid <- expand.grid(icc = ICCS_H, n = N_GRID_H, e = seq_len(nrow(MEASURED))) %>%
  mutate(
    effect = factor(MEASURED$effect[e], levels = MEASURED$effect),
    delta = MEASURED$delta[e],
    fc = MEASURED$fc[e],
    power = mapply(power_disease, n = n, k = K_REPS, icc = icc, delta = delta),
    icc_label = factor(ICC_LABELS_H[match(icc, ICCS_H)], levels = ICC_LABELS_H)
  )

cat(sprintf("Human AD against control, %d against %d donors, %d proteins\n\n",
            ad$n_group1, ad$n_group0, ad$proteins))
cat(sprintf("  %-26s %-6s %-6s %-24s %s  %s\n", "effect", "SD", "FC", "ICC",
            paste(sprintf("n=%-4d", c(3, 5, 10, 20, 50)), collapse = " "),
            "donors/group for 80%"))
for (i in seq_len(nrow(MEASURED))) {
  for (ic in ICCS_H) {
    r <- measured_grid[measured_grid$e == i & measured_grid$icc == ic &
                         measured_grid$n %in% c(3, 5, 10, 20, 50), ]
    cat(sprintf("  %-26s %-6.2f %-6.2f %-24s %s  %d\n", MEASURED$effect[i], MEASURED$delta[i],
                MEASURED$fc[i], ICC_LABELS_H[match(ic, ICCS_H)],
                paste(sprintf("%5.1f ", r$power * 100), collapse = " "),
                donors_required(K_REPS, ic, MEASURED$delta[i])))
  }
}

measured_mc_grid <- expand.grid(icc = ICCS_H, n = c(3, 5, 10, 20, 50), e = seq_len(nrow(MEASURED)))
measured_exag <- do.call(rbind, lapply(seq_len(nrow(measured_mc_grid)), function(i) {
  r <- simulate_contrast(measured_mc_grid$n[i], K_REPS, measured_mc_grid$icc[i],
                         MEASURED$delta[measured_mc_grid$e[i]])
  r$effect <- MEASURED$effect[measured_mc_grid$e[i]]
  r$fc <- MEASURED$fc[measured_mc_grid$e[i]]
  r
})) %>%
  mutate(
    fc_reported = 2^(exaggeration * log2(fc)),
    icc_label = ICC_LABELS_H[match(icc, ICCS_H)]
  )

cat("\nAmong significant results: exaggeration and the fold change reported\n\n")
cat(sprintf("  %-26s %-24s %-12s %8s %14s %14s\n",
            "effect", "ICC", "donors/group", "power", "exaggeration", "reported FC"))
for (i in seq_len(nrow(measured_exag))) {
  r <- measured_exag[i, ]
  cat(sprintf("  %-26s %-24s %-12d %7.1f%% %13.2fx %14.2f\n",
              r$effect, r$icc_label, r$n, r$power_mc * 100, r$exaggeration, r$fc_reported))
}

pD7 <- ggplot(measured_grid, aes(x = n, y = power, color = effect, group = effect)) +
  annotate("rect", xmin = 3, xmax = 5, ymin = 0, ymax = 1, fill = "grey50", alpha = 0.12) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  facet_wrap(~icc_label, nrow = 1) +
  scale_x_continuous(trans = "log10", breaks = N_GRID_H) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_color_manual(values = c("#2166AC", "#B2182B")) +
  labs(
    title = "Power at the effect sizes measured in human AD cortex",
    x = "Donors per group (log scale)",
    y = "Power",
    color = NULL
  ) +
  theme_pub

tiff("pics/figD_power_measured.tiff", units = "in", res = 300, width = 11, height = 5, compression = "lzw")
print(pD7)
dev.off()

write.csv(measured_grid %>% select(effect, delta, fc, icc, n, power),
          "tables/tableE7_power_measured_effects.csv", row.names = FALSE)
write.csv(measured_exag %>% select(effect, delta, fc, icc, n, power_mc, exaggeration, sign_error, fc_reported),
          "tables/tableE8_exaggeration_measured_effects.csv", row.names = FALSE)


################################################################################
# FIGURES AND TABLES
################################################################################

write.csv(fc_grid %>% select(cells, line_sd, icc, n, fc, delta, power),
          "tables/tableE5_power_by_fold_change.csv", row.names = FALSE)
write.csv(fc_exag %>% select(cells, line_sd, icc, n, fc, power_mc, exaggeration, sign_error, fc_reported),
          "tables/tableE6_exaggeration_by_fold_change.csv", row.names = FALSE)


write.csv(power_grid %>% select(icc, n, delta, delta_resid_sd, power),
          "tables/tableE1_power_disease_contrast.csv", row.names = FALSE)
write.csv(mde_grid %>% select(icc, n, mde, mde_resid_sd),
          "tables/tableE2_min_detectable_effect.csv", row.names = FALSE)
write.csv(fpr_grid %>% select(icc, n, fpr),
          "tables/tableE4_naive_false_positives.csv", row.names = FALSE)



################################################################################
# KEY STATISTICS FOR MANUSCRIPT TEXT
################################################################################

cat("\n--- Key statistics for manuscript text ---\n\n")
for (ic in c(0.20, 0.50)) {
  cat(sprintf("ICC = %.0f%%, %d replicates per donor:\n", ic * 100, K_REPS))
  for (nn in c(3, 5, 10, 20)) {
    cat(sprintf("  %d v %d donors: power %.0f%% at 0.5 SD, %.0f%% at 1 SD, %.0f%% at 2 SD; smallest effect at 80%% power %.2f SD\n",
                nn, nn,
                power_disease(nn, K_REPS, ic, 0.5) * 100,
                power_disease(nn, K_REPS, ic, 1) * 100,
                power_disease(nn, K_REPS, ic, 2) * 100,
                min_detectable_effect(nn, K_REPS, ic)))
  }
  cat(sprintf("  donors per group for 80%% power: %d at 1 SD, %d at 0.5 SD\n",
              donors_required(K_REPS, ic, 1), donors_required(K_REPS, ic, 0.5)))
  cat(sprintf("  false positives with replicates as independent: %.0f%% (1 v 1), %.0f%% (3 v 3)\n\n",
              fpr_grid$fpr[fpr_grid$icc == ic & fpr_grid$n == 1] * 100,
              fpr_grid$fpr[fpr_grid$icc == ic & fpr_grid$n == 3] * 100))
}
