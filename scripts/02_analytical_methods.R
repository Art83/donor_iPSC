################################################################################
# Analytical Arguments: Design Effect, Population Genetics, and Exact CIs
#
# These three analyses require NO simulation, they are closed-form results
# that complement the Monte Carlo framework.
#
# Part A: Design Effect / Effective Sample Size
# Part B: Population Genetics, Variant Capture by n Individuals
# Part C: Exact Confidence Intervals for ICC from F-distribution
#
################################################################################


library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)


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
# PART A: DESIGN EFFECT / EFFECTIVE SAMPLE SIZE
#
# The design effect (Kish, 1965) quantifies the loss of information due to
# clustering. When observations within a cluster (donor) are correlated,
# the effective sample size is smaller than the nominal sample size.
#
#   DEFF = 1 + (m - 1) * ICC
#   n_effective = n_total / DEFF
#
# where m = number of observations per cluster (organoids per donor)
# and ICC = intraclass correlation (donor-level variance / total variance).
#
# Consequences:
#   - p-values from naive tests (ignoring clustering) are anti-conservative
#   - confidence intervals are too narrow
#   - the "large n" from many organoids is misleading
#
################################################################################

cat("=== PART A: Design Effect ===\n\n")

# Parameter grid
n_donors_vec <- c(1, 2, 3, 5, 10)
m_per_donor  <- c(3, 6, 10, 20, 50, 100)
icc_vec      <- c(0.05, 0.10, 0.20, 0.30, 0.50)

# Compute design effect table
deff_data <- expand.grid(
  n_donors = n_donors_vec,
  m = m_per_donor,
  icc = icc_vec
) %>%
  mutate(
    n_total = n_donors * m,
    deff = 1 + (m - 1) * icc,
    n_effective = n_total / deff,
    # Ratio: what fraction of your data is "real"?
    efficiency = n_effective / n_total,
    # Labels
    icc_label = paste0("ICC = ", icc * 100, "%"),
    config = paste0(n_donors, " donors \u00D7 ", m, " org.")
  )

# Print key examples
cat("Key examples (what your n really is):\n\n")
cat(sprintf("%-30s  %6s  %6s  %10s  %s\n",
            "Design", "Total", "ICC", "Effective", "Efficiency"))
cat(paste(rep("-", 75), collapse = ""), "\n")

examples <- deff_data %>%
  filter(
    (n_donors == 3 & m %in% c(6, 20, 50) & icc %in% c(0.10, 0.20, 0.30)) |
    (n_donors == 1 & m %in% c(50, 100) & icc %in% c(0.20, 0.30)) |
    (n_donors == 5 & m %in% c(10, 50) & icc %in% c(0.20, 0.30)) |
    (n_donors == 10 & m %in% c(10, 50) & icc %in% c(0.20, 0.30))
  ) %>%
  arrange(n_donors, m, icc)

for (i in 1:nrow(examples)) {
  row <- examples[i, ]
  cat(sprintf("%-30s  %6d  %5.0f%%  %10.1f  %5.0f%%\n",
              row$config, row$n_total, row$icc * 100,
              row$n_effective, row$efficiency * 100))
}

# ---------- FIGURE A1: Effective sample size heatmap ----------

heatmap_3d <- deff_data %>%
  filter(n_donors == 3) %>%
  mutate(
    m_label = factor(paste0(m, " org/donor"),
                     levels = paste0(m_per_donor, " org/donor")),
    icc_pct = factor(paste0(icc * 100, "%"),
                     levels = paste0(icc_vec * 100, "%")),
    eff_label = sprintf("%.0f\n(of %d)", n_effective, n_total)
  )

pA1 <- ggplot(heatmap_3d, aes(x = icc_pct, y = m_label, fill = n_effective)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = eff_label), size = 3.8, lineheight = 0.85) +
  scale_fill_gradient2(
    low = "#E63946", mid = "#F4A261", high = "#264653",
    midpoint = 10, limits = c(0, NA)
  ) +
  labs(
    title = "Effective sample size for a 3-donor study",
    x = "Donor ICC (genetic contribution to variance)",
    y = "Organoids per donor",
    fill = "Effective\nsample size"
  ) +
  theme_pub +
  theme(panel.grid = element_blank(), axis.ticks = element_blank())

# ---------- FIGURE A2: Effective n across designs ----------

# Compare: what do you gain from more donors vs more organoids?
comparison_data <- bind_rows(
  # Strategy 1: add organoids (keep 3 donors)
  tibble(
    strategy = "Add organoids\n(3 donors, vary org/donor)",
    m = c(3, 6, 10, 20, 50, 100, 200),
    n_donors = 3,
    icc = 0.20
  ),
  # Strategy 2: add donors (keep 6 organoids)
  tibble(
    strategy = "Add donors\n(6 org/donor, vary donors)",
    m = 6,
    n_donors = c(3, 5, 8, 10, 15, 20, 30),
    icc = 0.20
  )
) %>%
  mutate(
    n_total = n_donors * m,
    deff = 1 + (m - 1) * icc,
    n_effective = n_total / deff
  )

pA2 <- ggplot(comparison_data,
              aes(x = n_total, y = n_effective, color = strategy)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60") +
  annotate("text", x = 400, y = 380,
           label = "Ideal: n_eff = n_total\n(no clustering)",
           color = "grey50", size = 3.5, hjust = 0.5, angle = 25) +
  scale_color_manual(values = c("#E63946", "#264653")) +
  labs(
    title = "Adding donors is far more efficient than adding organoids",
    x = "Total number of measurements",
    y = "Effective sample size",
    color = NULL
  ) +
  theme_pub +
  theme(legend.position = c(0.3, 0.85),
        legend.background = element_rect(fill = "white", color = "grey80"))

pA <- pA1 + pA2 +
  plot_annotation(
    tag_levels = "A",
    title = paste0("The design effect: why more organoids cannot compensate ",
                   "for fewer donors"),
    theme = theme(plot.title = element_text(face = "bold", size = 17))
  )

ggsave("figA_design_effect.pdf", pA, width = 16, height = 8, dpi = 300)
ggsave("figA_design_effect.png", pA, width = 16, height = 8, dpi = 300)



# ---------- FIGURE A3: Anti-conservative p-values ----------
# Show what happens to a t-test when you ignore clustering

cat("\nSimulating anti-conservative p-values from naive t-test\n")

set.seed(42)
n_sims <- 10000

# Scenario: 3 donors per group, 20 organoids per donor, ICC=0.20
# Treatment has NO effect (null is true)
# Compare: naive t-test on all 60 obs vs proper test on donor means

naive_p <- numeric(n_sims)
proper_p <- numeric(n_sims)

n_per_group <- 3
m_per <- 20
icc_sim <- 0.20
donor_sd <- sqrt(icc_sim / (1 - icc_sim))  # residual SD = 1

for (sim in 1:n_sims) {
  # Group A: 3 donors, no treatment effect
  donors_a <- rnorm(n_per_group, 0, donor_sd)
  y_a <- rep(donors_a, each = m_per) + rnorm(n_per_group * m_per)
  donor_id_a <- rep(1:n_per_group, each = m_per)

  # Group B: 3 donors, no treatment effect (null is true)
  donors_b <- rnorm(n_per_group, 0, donor_sd)
  y_b <- rep(donors_b, each = m_per) + rnorm(n_per_group * m_per)
  donor_id_b <- rep(1:n_per_group, each = m_per)

  # Naive t-test: ignores clustering, treats all organoids as independent
  naive_p[sim] <- t.test(y_a, y_b)$p.value

  # Proper test: average within donor first, then test donor means
  means_a <- tapply(y_a, donor_id_a, mean)
  means_b <- tapply(y_b, donor_id_b, mean)
  proper_p[sim] <- t.test(means_a, means_b)$p.value
}

cat(sprintf("Naive t-test: false positive rate = %.1f%% (should be 5%%)\n",
            mean(naive_p < 0.05) * 100))
cat(sprintf("Proper test: false positive rate = %.1f%% (should be 5%%)\n",
            mean(proper_p < 0.05) * 100))

fpr_data <- data.frame(
  p_value = c(naive_p, proper_p),
  method = rep(c("Naive t-test\n(organoids as independent)",
                 "Proper test\n(donor means)"),
               each = n_sims)
)

pA3 <- ggplot(fpr_data, aes(x = p_value, fill = method)) +
  geom_histogram(bins = 50, alpha = 0.85, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "#E63946",
             linewidth = 0.8) +
  facet_wrap(~method, ncol = 2) +
  scale_fill_manual(values = c("#E63946", "#264653")) +
  labs(
    title = "Ignoring donor clustering inflates false positive rates",
    x = "p-value",
    y = "Count",
    fill = NULL
  ) +
  theme_pub +
  theme(legend.position = "none")

ggsave("figA3_anticonservative_p.pdf", pA3, width = 12, height = 5, dpi = 300)
ggsave("figA3_anticonservative_p.png", pA3, width = 12, height = 5, dpi = 300)


# Print the false positive rate table
cat("\n--- False positive rates by design (naive t-test, null is true) ---\n\n")
cat(sprintf("%-40s  %10s  %10s\n", "Design", "Naive FPR", "Nominal"))
cat(paste(rep("-", 65), collapse = ""), "\n")

for (nd in c(2, 3, 5, 10)) {
  for (mp in c(10, 20, 50)) {
    for (ic in c(0.10, 0.20, 0.30)) {
      ds <- sqrt(ic / (1 - ic))
      np <- numeric(2000)
      for (s in 1:2000) {
        ya <- rep(rnorm(nd, 0, ds), each = mp) + rnorm(nd * mp)
        yb <- rep(rnorm(nd, 0, ds), each = mp) + rnorm(nd * mp)
        np[s] <- t.test(ya, yb)$p.value
      }
      fpr <- mean(np < 0.05)
      if (fpr > 0.08) {  # only print problematic ones
        cat(sprintf("%2d donors \u00D7 %3d org, ICC=%2.0f%%   %9.1f%%  %10s\n",
                    nd, mp, ic * 100, fpr * 100, "5%"))
      }
    }
  }
}


################################################################################
# PART B: POPULATION GENETICS, VARIANT CAPTURE
#
# Completely independent of ICC estimation. This asks:
# How much common genetic variation is even REPRESENTED when you
# sample n individuals?
#
# For a variant with minor allele frequency p:
#   P(observed in at least one copy) = 1 - (1-p)^(2n)
#
# Integrate over the site frequency spectrum to get the expected
# fraction of all common variants captured.
################################################################################

cat("\n\n=== PART B: Population Genetics ===\n\n")

# Proportion of variants captured at a given MAF, for n individuals
p_capture <- function(maf, n) {
  1 - (1 - maf)^(2 * n)
}

# Expected proportion of variants captured across a MAF range
# Assuming neutral SFS: density proportional to 1/p (folded)
expected_capture <- function(n, maf_min = 0.01, maf_max = 0.50,
                              n_points = 1000) {
  mafs <- seq(maf_min, maf_max, length.out = n_points)
  # Neutral folded SFS: f(p) proportional to 1/p for p in (0, 0.5)
  weights <- 1 / mafs
  weights <- weights / sum(weights)

  capture_probs <- p_capture(mafs, n)
  weighted_capture <- sum(capture_probs * weights)

  return(weighted_capture)
}

# Compute for various n
n_individuals <- c(1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000)

# Three MAF categories
capture_results <- data.frame()
for (n in n_individuals) {
  # Common variants (MAF 5-50%)
  common <- expected_capture(n, maf_min = 0.05, maf_max = 0.50)
  # Low-frequency (MAF 1-5%)
  low_freq <- expected_capture(n, maf_min = 0.01, maf_max = 0.05)
  # Rare (MAF 0.1-1%)
  rare <- expected_capture(n, maf_min = 0.001, maf_max = 0.01)

  capture_results <- rbind(capture_results, data.frame(
    n = n,
    common = common,
    low_freq = low_freq,
    rare = rare
  ))
}

cat("Fraction of genetic variants captured by n individuals:\n\n")
cat(sprintf("%6s  %10s  %10s  %10s\n",
            "n", "Common", "Low-freq", "Rare"))
cat(sprintf("%6s  %10s  %10s  %10s\n",
            "", "(5-50%)", "(1-5%)", "(0.1-1%)"))
cat(paste(rep("-", 42), collapse = ""), "\n")
for (i in 1:nrow(capture_results)) {
  cat(sprintf("%6d  %9.1f%%  %9.1f%%  %9.1f%%\n",
              capture_results$n[i],
              capture_results$common[i] * 100,
              capture_results$low_freq[i] * 100,
              capture_results$rare[i] * 100))
}

# ---------- FIGURE B1: Variant capture curves ----------

capture_long <- capture_results %>%
  pivot_longer(cols = c(common, low_freq, rare),
               names_to = "category", values_to = "fraction") %>%
  mutate(
    category = factor(category,
                      levels = c("rare", "low_freq", "common"),
                      labels = c("Rare (MAF 0.1-1%)",
                                 "Low-frequency (MAF 1-5%)",
                                 "Common (MAF 5-50%)"))
  )

pB1 <- ggplot(capture_long, aes(x = n, y = fraction, color = category)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey50") +
  annotate("text", x = 1.5, y = 0.98, label = "95% capture",
           color = "grey40", size = 3.5) +
  geom_vline(xintercept = c(2, 3, 5), linetype = "dotted",
             color = "grey70", linewidth = 0.3) +
  annotate("rect", xmin = 0.8, xmax = 5, ymin = 0, ymax = 1,
           fill = "red", alpha = 0.04) +
  annotate("text", x = 2.2, y = 0.85,
           label = "Typical iPSC study range",
           color = "#E63946", size = 3, fontface = "italic") +
  scale_x_continuous(trans = "log10",
                     breaks = c(1, 2, 3, 5, 10, 20, 50, 100, 500, 1000)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_color_manual(values = c("#E63946", "#F4A261", "#264653")) +
  labs(
    title = "Genetic variant capture as a function of sample size",
    x = "Number of individuals sampled (log scale)",
    y = "Fraction of variants captured",
    color = "Variant frequency class"
  ) +
  theme_pub + 
  theme(
    plot.title = element_text(size=14),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.8, "lines")
  )


tiff("../pics/fig_genetic_background.tiff", units = "in", res = 300, width = 7, height = 5)
print(pB1)
dev.off()



# ---------- FIGURE B2: What you miss at specific n ----------

maf_range <- seq(0.005, 0.50, by = 0.005)
miss_data <- expand.grid(maf = maf_range, n = c(2, 3, 5, 10, 50)) %>%
  mutate(
    p_miss = (1 - maf)^(2 * n),
    n_label = factor(paste0("n=", n),
                     levels = paste0("n=", c(2, 3, 5, 10, 50)))
  )

pB2 <- ggplot(miss_data, aes(x = maf, y = p_miss, color = n_label)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.4, y = 0.08, label = "5% miss rate",
           color = "grey40", size = 3.5) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(values = c("#E63946", "#F4A261", "#2A9D8F",
                                "#264653", "#1982C4")) +
  labs(
    title = "Probability of missing a variant entirely",
    x = "Minor allele frequency",
    y = "Probability variant is absent",
    color = "Donors sampled"
  ) +
  theme_pub

pB <- pB1 / pB2 +
  plot_annotation(
    tag_levels = "A",
    title = paste0("Small donor panels do not even sample the genetic ",
                   "variation they claim to address"),
    theme = theme(plot.title = element_text(face = "bold", size = 17))
  )

ggsave("figB_popgen_capture.pdf", pB, width = 11, height = 12, dpi = 300)
ggsave("figB_popgen_capture.png", pB, width = 11, height = 12, dpi = 300)

cat("\nFigure B saved.\n")

# Key stats for manuscript
cat("\n--- Key statistics for manuscript text ---\n\n")

for (n_val in c(2, 3, 5, 10)) {
  # P(miss) for a 5% MAF variant
  p_miss_5pct <- (1 - 0.05)^(2 * n_val)
  # P(miss) for a 10% MAF variant
  p_miss_10pct <- (1 - 0.10)^(2 * n_val)

  common_captured <- capture_results$common[capture_results$n == n_val]
  lowfreq_captured <- capture_results$low_freq[capture_results$n == n_val]

  cat(sprintf("n=%2d: miss 5%%-MAF variant %.0f%% | miss 10%%-MAF variant %.0f%% | ",
              n_val, p_miss_5pct * 100, p_miss_10pct * 100))
  cat(sprintf("common captured %.1f%% | low-freq captured %.1f%%\n",
              common_captured * 100, lowfreq_captured * 100))
}


################################################################################
# PART C: EXACT CONFIDENCE INTERVALS FOR ICC
#
# For balanced one-way random effects (Shrout & Fleiss, 1979):
#   F = MS_between / MS_within follows F(a-1, a(k-1)) under H0
#
# The relationship between F and ICC:
#   F = 1 + k * ICC / (1 - ICC)
#   ICC = (F - 1) / (F + k - 1)
#
# Exact 95% CI for ICC (Thomas & Hultquist, 1978):
#   F_L = F / F_{1-alpha/2, df1, df2}   gives ICC lower bound
#   F_U = F * F_{1-alpha/2, df2, df1}   gives ICC upper bound
#
# This is not a simulation result, it's the exact, minimum-possible
# CI width under ideal (balanced, normal) conditions.
################################################################################

cat("\n\n=== PART C: Exact Confidence Intervals for ICC ===\n\n")

# Exact CI for ICC given a, k, and a target ICC
# Returns the expected CI at the population F value
exact_icc_ci <- function(a, k, true_icc, alpha = 0.05) {
  df1 <- a - 1
  df2 <- a * (k - 1)

  # Expected F under the true ICC
  F_expected <- 1 + k * true_icc / (1 - true_icc)

  # CI for F (using the Searle approach)
  F_crit_upper <- qf(1 - alpha / 2, df1, df2)
  F_crit_lower <- qf(alpha / 2, df1, df2)

  # Transform to ICC scale
  F_L <- F_expected / F_crit_upper
  F_U <- F_expected / F_crit_lower

  icc_lower <- max((F_L - 1) / (F_L + k - 1), 0)
  icc_upper <- min((F_U - 1) / (F_U + k - 1), 1)

  return(list(
    lower = icc_lower,
    upper = icc_upper,
    width = icc_upper - icc_lower,
    a = a,
    k = k,
    true_icc = true_icc,
    F_expected = F_expected
  ))
}

# Build table
ci_data <- data.frame()
for (a in c(2, 3, 5, 10, 20, 50, 100)) {
  for (k in c(6)) {
    for (icc in c(0.05, 0.10, 0.20, 0.30, 0.50)) {
      res <- exact_icc_ci(a, k, icc)
      ci_data <- rbind(ci_data, data.frame(
        n_donors  = a,
        k         = k,
        true_icc  = icc,
        ci_lower  = res$lower,
        ci_upper  = res$upper,
        ci_width  = res$width
      ))
    }
  }
}

cat("Exact 95% CI width for ICC under ideal conditions (k=6 reps/donor):\n\n")
cat(sprintf("%8s  %10s  %12s  %20s  %10s\n",
            "Donors", "True ICC", "CI Width",
            "95% CI", "Useful?"))
cat(paste(rep("-", 68), collapse = ""), "\n")

for (i in 1:nrow(ci_data)) {
  row <- ci_data[i, ]
  useful <- ifelse(row$ci_width < 0.20, "Possibly",
                   ifelse(row$ci_width < 0.40, "Marginal", "No"))
  cat(sprintf("%8d  %9.0f%%  %11.0f pp  [%5.1f%%, %5.1f%%]  %10s\n",
              row$n_donors, row$true_icc * 100, row$ci_width * 100,
              row$ci_lower * 100, row$ci_upper * 100, useful))
}

# ---------- FIGURE C1: CI width by sample size ----------

ci_plot_data <- ci_data %>%
  mutate(
    icc_label = factor(paste0(true_icc * 100, "%"),
                       levels = paste0(c(5, 10, 20, 30, 50), "%"))
  )

palette_icc <- c(
  "5%" = "#BFDBFE",
  "10%" = "#93C5FD",
  "20%" = "#3B82F6",
  "30%" = "#1D4ED8",
  "50%" = "#1E3A5F"
)

pC1 <- ggplot(ci_plot_data, aes(x = n_donors, y = ci_width,
                                 color = icc_label)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  # geom_hline(yintercept = 0.20, linetype = "dashed", color = "grey50") +
  # annotate("text", x = 80, y = 0.23,
  #          label = "20 pp width\n(minimally informative)",
  #          color = "grey40", size = 3, hjust = 1) +
  annotate("rect", xmin = 1.5, xmax = 5, ymin = 0, ymax = 1,
           fill = "red", alpha = 0.04) +
  scale_x_continuous(trans = "log10",
                     breaks = c(2, 3, 5, 10, 20, 50, 100)) +
  scale_y_continuous(labels = function(x) paste0(x * 100, " pp"),
                     limits = c(0, 1)) +
  scale_color_manual(values = palette_icc) +
  labs(
    title = "Exact 95% CI width for ICC under ideal conditions",
    x = "Number of donors (log scale)",
    y = "95% CI width (percentage points)",
    color = "True ICC"
  ) +
  theme_pub + 
  theme(
    plot.title = element_text(size=14)
  )


tiff("../pics/figC_CIs_1.tiff", units = "in", res = 300, width = 6, height = 5)
print(pC1)
dev.off()




# ---------- FIGURE C2: CI bands ----------

pC2 <- ggplot(ci_plot_data %>% filter(true_icc == 0.20),
              aes(x = n_donors)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              fill = "#3B82F6", alpha = 0.25) +
  geom_line(aes(y = true_icc), color = "#E63946", linetype = "dashed",
            linewidth = 0.8) +
  geom_line(aes(y = ci_lower), color = "#3B82F6", linewidth = 0.8) +
  geom_line(aes(y = ci_upper), color = "#3B82F6", linewidth = 0.8) +
  annotate("text", x = 80, y = 0.24, label = "True ICC = 20%",
           color = "#E63946", size = 3.5, fontface = "bold") +
  scale_x_continuous(trans = "log10",
                     breaks = c(2, 3, 5, 10, 20, 50, 100)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title = "Exact 95% confidence band for ICC = 20%",
    x = "Number of donors (log scale)",
    y = "ICC"
  ) +
  theme_pub


tiff("../pics/figC_CIs_2.tiff", units = "in", res = 300, width = 6, height = 5)
print(pC2)
dev.off()


pC <- pC1 / pC2 +
  plot_annotation(
    tag_levels = "A",
    title = paste0("Mathematical proof: ICC cannot be meaningfully estimated ",
                   "at small n"),
    subtitle = paste0("Exact confidence intervals from F-distribution | ",
                      "No simulation assumptions"),
    theme = theme(
      plot.title = element_text(face = "bold", size = 17),
      plot.subtitle = element_text(size = 12, color = "grey40")
    )
  )

ggsave("figC_exact_ci.pdf", pC, width = 11, height = 12, dpi = 300)
ggsave("figC_exact_ci.png", pC, width = 11, height = 12, dpi = 300)

cat("\nFigure C saved.\n")


################################################################################
# SUMMARY TABLE: All three analytical arguments combined
################################################################################

cat("\n\n")
cat("################################################################\n")
cat("# COMBINED SUMMARY FOR MANUSCRIPT                             #\n")
cat("################################################################\n\n")

cat("For a study with 3 donors, 6 organoids per donor, ICC = 20%:\n\n")

deff_3 <- 1 + (6 - 1) * 0.20
neff_3 <- 18 / deff_3
ci_3 <- exact_icc_ci(3, 6, 0.20)
common_3 <- capture_results$common[capture_results$n == 3]
lowfreq_3 <- capture_results$low_freq[capture_results$n == 3]
miss_5pct_3 <- (1 - 0.05)^6

cat(sprintf("DESIGN EFFECT:DEFF = %.1f | effective n = %.1f (of 18 nominal)\n",
            deff_3, neff_3))
cat(sprintf("EXACT CI: 95%% CI for ICC = [%.1f%%, %.1f%%] (width = %.0f pp)\n",
            ci_3$lower * 100, ci_3$upper * 100, ci_3$width * 100))
cat(sprintf("VARIANT CAPTURE: %.1f%% of common variants, %.1f%% of low-frequency\n",
            common_3 * 100, lowfreq_3 * 100))
cat(sprintf("%.0f%% chance of completely missing a 5%%-MAF variant\n",
            miss_5pct_3 * 100))


