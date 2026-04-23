################################################################################
#
# This script demonstrates through simulation that typical iPSC/organoid study
# designs (n = 2-10 donors) are critically underpowered to:
#   detect donor-level genetic effects,
#   estimate their magnitude with any useful precision,
#   disentangle donor variance from differentiation-batch variance.
# Thus, making widely propagated N(donors)=3-5 requirements completely random and no different than N=1 
#
# ICC values are grounded in published large-scale variance decomposition studies:
#   - Kilpinen et al., Nature, 2017 (HipSci: 711 lines, 301 donors):
#       donor explains 5–46% of variance depending on assay; 46% of expression
#       probes have donor as dominant factor; morphology ICC up to 23%.
#   - Carcamo-Orive et al., Cell Stem Cell, 2017 (317 lines, 101 individuals):
#       ~50% of genome-wide expression variability explained by inter-individual
#       variation.
#   - Mirauta et al., eLife, 2020 (HipSci proteomics, 202 lines, 151 donors):
#       donor variance component ranges broadly across proteins.
#   - Cuomo et al., Nat Commun, 2020 (125 donors, scRNA-seq endoderm):
#       dynamic eQTLs during differentiation, donor effects on cell composition.
#
# Empirically-grounded ICC range (based on Kilpinen) used: 0.05, 0.10, 0.20, 0.30, 0.50
#
# Uses one-way and two-way ANOVA decompositions (exact F-tests, no iterative
# optimisation, no singular fits, lme4 was attempted but it struggles due to singularity). For the batch-confounding analysis, 
# a nested two-factor ANOVA is employed.
#
# Gelman & Carlin (2014) Type M and Type S error framework is applied.
#
################################################################################


library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

set.seed(42)

# Empirically-grounded ICC values (see header citations)
TRUE_ICCS <- c(0.05, 0.10, 0.20, 0.30, 0.50)

# Donor sample sizes spanning typical to aspirational
N_DONORS <- c(2, 3, 5, 10, 20, 50)

# Replicates per donor (organoids per line, typical for functional assays)
K_REPS <- 6

# Number of Monte Carlo iterations
N_SIMS <- 5000

# Population size for sampling donor effects
N_POP <- 1000

# Residual SD (set to 1; ICC controls the ratio)
RESID_SD <- 1

# --- Color palette ---
palette_n <- c(
  "n=2" = "#E63946",
  "n=3" = "#F4A261",
  "n=5" = "#2A9D8F",
  "n=10" = "#264653",
  "n=20" = "#6A4C93",
  "n=50" = "#1982C4"
)

palette_icc <- c(
  "5%" = "#BFDBFE",
  "10%" = "#93C5FD",
  "20%" = "#3B82F6",
  "30%" = "#1D4ED8",
  "50%" = "#1E3A5F"
)

# --- Theme for publication ---
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
# CORE FUNCTION: Vectorised one-way random-effects ANOVA
#
# Shrout & Fleiss (1979) ICC(1,1):
#   ICC = (MS_between - MS_within) / (MS_between + (k-1) * MS_within)
#   Truncated at 0 when MS_between < MS_within.
#
# F = MS_between / MS_within,  df1 = a-1, df2 = a(k-1)
################################################################################

anova_icc_fast <- function(y_matrix, k, a) {
  N <- a * k
  n_sims <- ncol(y_matrix)
  
  grand_means <- colMeans(y_matrix)
  
  dim(y_matrix) <- c(k, a, n_sims)
  group_means <- colMeans(y_matrix, dims = 1)  # (a x n_sims)
  
  ss_between <- k * colSums((group_means - rep(grand_means, each = a))^2)
  
  group_means_expanded <- group_means[rep(1:a, each = k), ]
  dim(y_matrix) <- c(N, n_sims)
  ss_within <- colSums((y_matrix - group_means_expanded)^2)
  
  ms_between <- ss_between / (a - 1)
  ms_within <- ss_within  / (N - a)
  
  f_stat <- ms_between / ms_within
  p_value <- pf(f_stat, df1 = a - 1, df2 = N - a, lower.tail = FALSE)
  
  icc_est <- (ms_between - ms_within) / (ms_between + (k - 1) * ms_within)
  icc_est <- pmax(icc_est, 0)
  
  return(list(
    f_stat = f_stat,
    p_value = p_value,
    icc = icc_est,
    ms_between = ms_between,
    ms_within = ms_within
  ))
}


################################################################################
# SIMULATION 1: Power & ICC estimation across empirical ICC range
################################################################################

cat("=== SIMULATION 1: Power & ICC Estimation ===\n")

simulate_power <- function(n_donors_vec = N_DONORS,
                           true_iccs = TRUE_ICCS,
                           k = K_REPS,
                           n_sims = N_SIMS,
                           n_pop = N_POP,
                           resid_sd = RESID_SD) {
  results <- data.frame()
  
  for (true_icc in true_iccs) {
    donor_sd <- sqrt((true_icc * resid_sd^2) / (1 - true_icc))
    
    cat(sprintf("\nICC = %.0f%% (donor SD = %.3f, residual SD = %.3f)\n",
                true_icc * 100, donor_sd, resid_sd))
    
    pop_effects <- rnorm(n_pop, mean = 0, sd = donor_sd)
    
    for (a in n_donors_vec) {
      cat(sprintf("n = %d donors: ", a))
      
      y_matrix <- matrix(NA_real_, nrow = a * k, ncol = n_sims)
      
      for (sim in 1:n_sims) {
        sampled <- sample(pop_effects, a)
        donor_effect <- rep(sampled, each = k)
        noise <- rnorm(a * k, 0, resid_sd)
        y_matrix[, sim] <- donor_effect + noise
      }
      
      res <- anova_icc_fast(y_matrix, k = k, a = a)
      
      power <- mean(res$p_value < 0.05)
      cat(sprintf("power = %.1f%%\n", power * 100))
      
      results <- rbind(results, data.frame(
        true_icc       = true_icc,
        n_donors       = a,
        power          = power,
        mean_est_icc   = mean(res$icc),
        median_est_icc = median(res$icc),
        sd_est_icc     = sd(res$icc),
        q025_icc       = quantile(res$icc, 0.025),
        q975_icc       = quantile(res$icc, 0.975),
        pct_icc_zero   = mean(res$icc == 0),
        icc_values     = I(list(res$icc)),
        p_values       = I(list(res$p_value))
      ))
    }
  }
  
  return(results)
}

power_results <- simulate_power()


################################################################################
# SIMULATION 2: Type M and Type S Errors (Gelman & Carlin, 2014)
#
# Among simulations where p < 0.05:
#   Type M (exaggeration ratio) = |estimated ICC| / true ICC
#   Type S (sign error rate) = proportion where estimated ICC is in wrong
#                                  direction (for variance components, this means
#                                  zero due to truncation)
#
# The key insight: at small n, your *significant* results are not just noisy,
# they systematically overestimate the true effect."
################################################################################

cat("\n=== SIMULATION 2: Type M / Type S Errors ===\n")

compute_type_m_s <- function(power_results) {
  type_m_results <- data.frame()
  
  for (i in 1:nrow(power_results)) {
    row <- power_results[i, ]
    icc_vals <- row$icc_values[[1]]
    p_vals <- row$p_values[[1]]
    
    # Only consider significant results
    sig_mask <- p_vals < 0.05
    
    if (sum(sig_mask) > 10) {
      sig_icc <- icc_vals[sig_mask]
      
      # Type M: exaggeration ratio (estimated / true)
      exaggeration <- sig_icc / row$true_icc
      type_m_median <- median(exaggeration)
      type_m_mean <- mean(exaggeration)
      
      # Distribution of significant estimates
      q25_sig <- quantile(sig_icc, 0.25)
      q75_sig <- quantile(sig_icc, 0.75)
      median_sig <- median(sig_icc)
      
      type_m_results <- rbind(type_m_results, data.frame(
        true_icc = row$true_icc,
        n_donors = row$n_donors,
        power = row$power,
        n_significant = sum(sig_mask),
        median_sig_icc = median_sig,
        q25_sig_icc = q25_sig,
        q75_sig_icc = q75_sig,
        type_m_median = type_m_median,
        type_m_mean = type_m_mean,
        pct_over_2x = mean(exaggeration > 2),
        pct_over_3x   = mean(exaggeration > 3)
      ))
    }
  }
  
  return(type_m_results)
}

type_m_results <- compute_type_m_s(power_results)

cat("\nType M exaggeration ratios (median estimated ICC / true ICC among sig results):\n")
print(type_m_results %>%
        select(true_icc, n_donors, power, type_m_median, pct_over_2x) %>%
        mutate(across(where(is.numeric), ~round(., 3))))


################################################################################
# SIMULATION 3: Batch-Donor Confounding
#
# Realistic organoid hierarchy: donor -> batch -> organoid (replicate)
#
# With n_batches_per_donor = 1 (common in practice), batch and donor effects
# are completely confounded. Even with n_batches > 1, at small n_donors, the
# model cannot reliably separate sigma2_donor from sigma2_batch.
#
# We generate data with BOTH donor and batch variance, then show:
#   A one-way ANOVA (ignoring batch) inflates the apparent donor ICC.
#   A proper nested model cannot stably estimate either component.
################################################################################

cat("\n=== SIMULATION 3: Batch-Donor Confounding ===\n")

simulate_confounding <- function(n_donors_vec = c(2, 3, 5, 10, 20, 50),
                                 n_batches_per_donor = 2,
                                 n_org_per_batch = 3,
                                 icc_donor = 0.20,
                                 icc_batch = 0.15,
                                 n_sims = 3000,
                                 resid_sd = 1) {
  
  icc_resid <- 1 - icc_donor - icc_batch
  total_var <- resid_sd^2 / icc_resid
  donor_sd <- sqrt(icc_donor * total_var)
  batch_sd <- sqrt(icc_batch * total_var)
  
  cat(sprintf("True: ICC_donor=%.0f%%, ICC_batch=%.0f%%, ICC_resid=%.0f%%\n",
              icc_donor*100, icc_batch*100, icc_resid*100))
  cat(sprintf("SDs:  donor=%.3f, batch=%.3f, resid=%.3f\n",
              donor_sd, batch_sd, resid_sd))
  
  results <- data.frame()
  
  for (a in n_donors_vec) {
    b <- n_batches_per_donor
    m <- n_org_per_batch
    N <- a * b * m
    
    icc_oneway <- numeric(n_sims)  # naive one-way (donor only)
    icc_donor_nested <- numeric(n_sims)  # from nested model
    icc_batch_nested <- numeric(n_sims)  # from nested model
    p_oneway <- numeric(n_sims)
    
    for (sim in 1:n_sims) {
      # Generate hierarchical data
      donor_effects <- rnorm(a, 0, donor_sd)
      batch_effects <- rnorm(a * b, 0, batch_sd)
      resid_noise <- rnorm(N, 0, resid_sd)
      
      # Build y: donor_i + batch_ij + epsilon_ijk
      donor_idx <- rep(1:a, each = b * m)
      batch_idx <- rep(1:(a*b), each = m)
      y <- donor_effects[donor_idx] + batch_effects[batch_idx] + resid_noise
      
      # --- One-way ANOVA (ignoring batch, as commonly done) ---
      group_means_d <- tapply(y, donor_idx, mean)
      grand_mean <- mean(y)
      ss_between_d <- b * m * sum((group_means_d - grand_mean)^2)
      ss_within_d <- sum((y - group_means_d[donor_idx])^2)
      ms_between_d <- ss_between_d / (a - 1)
      ms_within_d <- ss_within_d / (N - a)
      
      f_d <- ms_between_d / ms_within_d
      p_oneway[sim] <- pf(f_d, a - 1, N - a, lower.tail = FALSE)
      
      icc_naive <- (ms_between_d - ms_within_d) /
        (ms_between_d + (b*m - 1) * ms_within_d)
      icc_oneway[sim] <- max(icc_naive, 0)
      
      # --- Nested ANOVA: batch nested within donor ---
      group_means_b <- tapply(y, batch_idx, mean)
      
      ss_donor <- b * m * sum((group_means_d - grand_mean)^2)
      ss_batch_in_donor <- m * sum((group_means_b -
                                      group_means_d[rep(1:a, each = b)])^2)
      ss_within <- sum((y - group_means_b[batch_idx])^2)
      
      df_donor <- a - 1
      df_batch <- a * (b - 1)
      df_within <- a * b * (m - 1)
      
      ms_donor_n <- ss_donor / df_donor
      ms_batch_n <- ss_batch_in_donor / df_batch
      ms_within_n <- ss_within / df_within
      
      # Variance component estimates (ANOVA method of moments)
      sigma2_within_est <- ms_within_n
      sigma2_batch_est <- (ms_batch_n - ms_within_n) / m
      sigma2_donor_est <- (ms_donor_n - ms_batch_n) / (b * m)
      
      # Truncate at zero
      sigma2_batch_est <- max(sigma2_batch_est, 0)
      sigma2_donor_est <- max(sigma2_donor_est, 0)
      
      total_est <- sigma2_donor_est + sigma2_batch_est + sigma2_within_est
      icc_donor_nested[sim] <- sigma2_donor_est / total_est
      icc_batch_nested[sim] <- sigma2_batch_est / total_est
    }
    
    results <- rbind(results, data.frame(
      n_donors = a,
      # One-way (naive) results
      mean_icc_oneway = mean(icc_oneway),
      median_icc_oneway = median(icc_oneway),
      q025_icc_oneway = quantile(icc_oneway, 0.025),
      q975_icc_oneway = quantile(icc_oneway, 0.975),
      power_oneway = mean(p_oneway < 0.05),
      # Nested model: donor component
      mean_icc_donor_nested = mean(icc_donor_nested),
      median_icc_donor_nested = median(icc_donor_nested),
      q025_icc_donor_nested = quantile(icc_donor_nested, 0.025),
      q975_icc_donor_nested = quantile(icc_donor_nested, 0.975),
      pct_donor_zero_nested = mean(icc_donor_nested == 0),
      # Nested model: batch component
      mean_icc_batch_nested = mean(icc_batch_nested),
      median_icc_batch_nested = median(icc_batch_nested),
      # Store distributions
      icc_oneway_vals = I(list(icc_oneway)),
      icc_donor_nested_vals = I(list(icc_donor_nested)),
      icc_batch_nested_vals = I(list(icc_batch_nested))
    ))
    
    cat(sprintf("n=%2d: naive ICC=%.1f%% [%.1f–%.1f] | nested donor=%.1f%% | nested batch=%.1f%%\n",
                a,
                mean(icc_oneway)*100,
                quantile(icc_oneway, 0.025)*100,
                quantile(icc_oneway, 0.975)*100,
                mean(icc_donor_nested)*100,
                mean(icc_batch_nested)*100))
  }
  
  return(results)
}

confound_results <- simulate_confounding()

cat("\n=== Generating Figures ===\n")

# ---------- FIGURE 1: Power curves across empirical ICC range ----------

fig1_data <- power_results %>%
  mutate(
    icc_label = factor(paste0(true_icc * 100, "%"),
                       levels = paste0(TRUE_ICCS * 100, "%"))
  )

p1 <- ggplot(fig1_data, aes(x = n_donors, y = power, color = icc_label)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "grey50",
             linewidth = 0.5) +
  annotate("text", x = 4.5, y = 0.83, label = "80% power threshold",
           color = "grey40", size = 3.5, hjust = 1) +
  geom_vline(xintercept = c(2, 3, 5, 10), linetype = "dotted",
             color = "grey70", linewidth = 0.3) +
  annotate("rect", xmin = 1.5, xmax = 5, ymin = 0, ymax = 1,
           fill = "red", alpha = 0.04) +
  annotate("text", x = 3, y = 0.97, label = "Typical iPSC study range",
           color = "#E63946", size = 3.5, fontface = "italic") +
  scale_x_continuous(breaks = c(2, 3, 5, 10, 20, 50),
                     trans = "log10") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_color_manual(values = palette_icc) +
  labs(
    title = "Power to detect donor-level genetic effects",
    x = "Number of donors (log scale)",
    y = "Statistical power",
    color = "True donor ICC\n(empirical range)"
  ) +
  theme_pub


tiff("../pics/fig1_power_donor_level.tiff", units = "in", res = 300, width = 6, height = 5)
print(p1)
dev.off()



# ---------- FIGURE 2: ICC estimate distributions ----------
# Show distributions at THREE representative ICCs: 0.05, 0.20, 0.50

fig2_iccs <- c(0.05, 0.20, 0.50)

icc_dist_data <- power_results %>%
  filter(true_icc %in% fig2_iccs) %>%
  select(true_icc, n_donors, icc_values) %>%
  unnest(icc_values) %>%
  rename(estimated_icc = icc_values) %>%
  mutate(
    n_label = factor(paste0("n=", n_donors),
                     levels = paste0("n=", N_DONORS)),
    icc_facet = factor(paste0("True ICC = ", true_icc * 100, "%"),
                       levels = paste0("True ICC = ", fig2_iccs * 100, "%"))
  )

p2 <- ggplot(icc_dist_data,
             aes(x = estimated_icc, fill = n_label)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, alpha = 0.85,
                 color = "white", linewidth = 0.2) +
  geom_vline(data = data.frame(
    icc_facet = factor(paste0("True ICC = ", fig2_iccs * 100, "%"),
                       levels = paste0("True ICC = ", fig2_iccs * 100, "%")),
    xint = fig2_iccs),
    aes(xintercept = xint), linetype = "dashed", color = "#E63946",
    linewidth = 0.8) +
  facet_grid(icc_facet ~ n_label, scales = "free_y") +
  scale_fill_manual(values = unname(palette_n)) +
  scale_x_continuous(labels = c(0,25,50,75,100), limits = c(-0.01, 1)) +
  scale_y_continuous(breaks = c(0, 5, 10),  limits = c(0,10))+
  labs(
    title = "Distribution of estimated donor ICC across simulations",
    x = "Estimated ICC",
    y = "Density"
  ) +
  theme_pub +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 9),
        strip.text.y = element_text(size = 9))

tiff("../pics/fig2_distribution_donor_level.tiff", units = "in", res = 300, width = 7.5, height = 5)
print(p2)
dev.off()



# ---------- FIGURE 3: Type M Error (Exaggeration) ----------

p3_data <- type_m_results %>%
  filter(true_icc %in% c(0.05, 0.10, 0.20, 0.30)) %>%
  mutate(
    icc_label = factor(paste0("True ICC = ", true_icc * 100, "%"),
                       levels = paste0("True ICC = ",
                                       c(5, 10, 20, 30), "%"))
  )

p3 <- ggplot(p3_data, aes(x = n_donors, y = type_m_median,
                          color = icc_label)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  annotate("text", x = 45, y = 1.15, label = "No exaggeration (ratio = 1)",
           color = "grey40", size = 3.5, hjust = 1) +
  annotate("rect", xmin = 1.5, xmax = 5, ymin = 0, ymax = Inf,
           fill = "red", alpha = 0.04) +
  scale_x_continuous(breaks = c(2, 3, 5, 10, 20, 50),
                     trans = "log10") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(values = c("#BFDBFE", "#93C5FD", "#3B82F6", "#1D4ED8")) +
  labs(
    title = "Type M error: significant results massively overestimate the truth",
    x = "Number of donors (log scale)",
    y = "Exaggeration ratio\n(estimated / true ICC)",
    color = NULL
  ) +
  theme_pub

ggsave("fig3_type_m_error.pdf", p3, width = 10, height = 6, dpi = 300)
ggsave("fig3_type_m_error.png", p3, width = 10, height = 6, dpi = 300)


# ---------- FIGURE 4: False Reassurance ----------
# Among non-significant results, the true ICC could still be substantial

reassurance_data <- power_results %>%
  filter(true_icc == 0.20) %>%
  select(n_donors, icc_values, p_values) %>%
  mutate(
    combined = mapply(function(icc, pv) {
      data.frame(estimated_icc = icc, p_value = pv)
    }, icc_values, p_values, SIMPLIFY = FALSE)
  ) %>%
  select(n_donors, combined) %>%
  unnest(combined) %>%
  mutate(
    significant = p_value < 0.05,
    n_label = factor(paste0("n=", n_donors),
                     levels = paste0("n=", N_DONORS))
  )

nonsig_data <- reassurance_data %>% filter(!significant)

nonsig_pct <- reassurance_data %>%
  group_by(n_donors, n_label) %>%
  summarize(pct_nonsig = mean(!significant), .groups = "drop")

p4 <- ggplot(nonsig_data %>% filter(n_donors <= 10),
             aes(x = estimated_icc, fill = n_label)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.85,
                 color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.20, linetype = "dashed", color = "#E63946",
             linewidth = 0.8) +
  geom_text(data = nonsig_pct %>% filter(n_donors <= 10),
            aes(x = 0.7, y = Inf,
                label = sprintf("%.0f%% of tests\nnon-significant",
                                pct_nonsig * 100)),
            vjust = 2, size = 3.5, fontface = "bold", color = "grey30",
            inherit.aes = FALSE) +
  facet_wrap(~n_label, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = unname(palette_n[1:4])) +
  scale_x_continuous(labels = c(0,25,50,75,100), limits = c(-0.01, 1)) +
  labs(
    title = "False reassurance",
    x = "Estimated ICC among non-significant results",
    y = "Density"
  ) +
  theme_pub +
  theme(legend.position = "none")

tiff("../pics/fig4_false_reassurance.tiff", units = "in", res = 300, width = 5, height = 5)
print(p4)
dev.off()


# ---------- FIGURE 5: Precision and bias across full ICC range ----------

precision_data <- power_results %>%
  mutate(
    n_label = factor(paste0("n=", n_donors),
                     levels = paste0("n=", N_DONORS)),
    icc_label = factor(paste0(true_icc * 100, "%"),
                       levels = paste0(TRUE_ICCS * 100, "%"))
  )

p5a <- ggplot(precision_data,
              aes(x = n_label, y = mean_est_icc, color = icc_label)) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = q025_icc, ymax = q975_icc),
                width = 0.3, linewidth = 0.7,
                position = position_dodge(width = 0.6)) +
  geom_hline(data = data.frame(true_icc = TRUE_ICCS,
                               icc_label = factor(paste0(TRUE_ICCS*100, "%"),
                                                  levels = paste0(TRUE_ICCS*100, "%"))),
             aes(yintercept = true_icc, color = icc_label),
             linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_color_manual(values = palette_icc) +
  labs(
    title = "Precision of donor ICC estimates by sample size",
    x = "Number of donors",
    y = "Estimated ICC",
    color = "True ICC"
  ) +
  theme_pub + 
  theme(plot.title = element_text(size = 14))

bias_data <- power_results %>%
  mutate(
    bias = mean_est_icc - true_icc,
    icc_label = factor(paste0(true_icc * 100, "%"),
                       levels = paste0(TRUE_ICCS * 100, "%"))
  )

p5b <- ggplot(bias_data, aes(x = n_donors, y = bias, color = icc_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_x_continuous(breaks = N_DONORS, trans = "log10") +
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(values = palette_icc) +
  labs(
    title = "Estimation bias: ANOVA-based ICC is downward-biased at small n",
    subtitle = "Due to truncation at zero, small samples systematically underestimate the donor contribution",
    x = "Number of donors (log scale)",
    y = "Bias (estimated \u2212 true ICC)",
    color = "True ICC"
  ) +
  theme_pub

p5 <- p5a / p5b + plot_annotation(
  title = "Small donor panels produce imprecise and biased ICC estimates",
  theme = theme(plot.title = element_text(face = "bold", size = 16))
)

tiff("../pics/fig5a_precision_donor_ICC.tiff", units = "in", res = 300, width = 6, height = 5)
print(p5a)
dev.off()



# ---------- FIGURE 6: Batch-Donor Confounding ----------

# Panel A: Naive vs nested ICC estimates
confound_long <- confound_results %>%
  select(n_donors, mean_icc_oneway, mean_icc_donor_nested,
         q025_icc_oneway, q975_icc_oneway,
         q025_icc_donor_nested, q975_icc_donor_nested) %>%
  pivot_longer(
    cols = c(mean_icc_oneway, mean_icc_donor_nested),
    names_to = "method",
    values_to = "mean_icc"
  ) %>%
  mutate(
    q025 = ifelse(method == "mean_icc_oneway",
                  q025_icc_oneway, q025_icc_donor_nested),
    q975 = ifelse(method == "mean_icc_oneway",
                  q975_icc_oneway, q975_icc_donor_nested),
    method_label = ifelse(method == "mean_icc_oneway",
                          "Na\u00EFve one-way ANOVA\n(batch ignored)",
                          "Nested ANOVA\n(batch modelled)"),
    n_label = factor(paste0("n=", n_donors),
                     levels = paste0("n=", N_DONORS))
  )

p6a <- ggplot(confound_long,
              aes(x = n_label, y = mean_icc, fill = method_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = q025, ymax = q975),
                position = position_dodge(width = 0.7),
                width = 0.2, linewidth = 0.5) +
  geom_hline(yintercept = 0.20, linetype = "dashed", color = "#E63946") +
  annotate("text", x = 6.5, y = 0.22, label = "True donor ICC = 20%",
           color = "#E63946", size = 3.5, fontface = "bold") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_fill_manual(values = c("#E63946", "#264653")) +
  labs(
    title = "Batch-donor confounding inflates apparent genetic effects",
    x = "Number of donors",
    y = "Estimated donor ICC",
    fill = NULL
  ) +
  theme_pub +
  theme(legend.position = "top")

# Panel B: Distribution of nested estimates at n=3 and n=10
confound_dist_data <- confound_results %>%
  filter(n_donors %in% c(3, 10)) %>%
  select(n_donors, icc_donor_nested_vals, icc_batch_nested_vals) %>%
  mutate(
    donor_vals = icc_donor_nested_vals,
    batch_vals = icc_batch_nested_vals
  ) %>%
  select(n_donors, donor_vals, batch_vals)

# Unnest both
dist_donor <- confound_dist_data %>%
  select(n_donors, donor_vals) %>%
  unnest(donor_vals) %>%
  rename(icc = donor_vals) %>%
  mutate(component = "Donor (genetic)")

dist_batch <- confound_dist_data %>%
  select(n_donors, batch_vals) %>%
  unnest(batch_vals) %>%
  rename(icc = batch_vals) %>%
  mutate(component = "Batch (technical)")

dist_combined <- bind_rows(dist_donor, dist_batch) %>%
  mutate(
    n_label = factor(paste0("n=", n_donors, " donors"),
                     levels = c("n=3 donors", "n=10 donors"))
  )

p6b <- ggplot(dist_combined, aes(x = icc, fill = component)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, alpha = 0.7,
                 position = "identity", color = "white", linewidth = 0.2) +
  geom_vline(data = data.frame(
    component = c("Donor (genetic)", "Batch (technical)"),
    true_val  = c(0.20, 0.15)),
    aes(xintercept = true_val, color = component),
    linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~n_label, ncol = 2) +
  scale_fill_manual(values = c("Batch (technical)" = "#F4A261",
                               "Donor (genetic)" = "#264653")) +
  scale_color_manual(values = c("Batch (technical)" = "#F4A261",
                                "Donor (genetic)" = "#264653"),
                     guide = "none") +
  scale_x_continuous(labels = scales::percent, limits = c(-0.01, 0.8)) +
  labs(
    title = "Variance component estimates from nested model",
    subtitle = paste0("At n=3, donor and batch components are effectively ",
                      "indistinguishable\n",
                      "Dashed lines = true values"),
    x = "Estimated ICC",
    y = "Density",
    fill = "Variance component"
  ) +
  theme_pub

p6 <- p6a / p6b + plot_annotation(
  tag_levels = "A",
  title = paste0("Ignoring batch structure inflates donor variance; ",
                 "accounting for it reveals estimation is futile at small n"),
  theme = theme(plot.title = element_text(face = "bold", size = 15))
)

ggsave("fig6_batch_confounding.pdf", p6, width = 12, height = 12, dpi = 300)



# ---------- FIGURE 7: Bayesian Update full ICC panel ----------

bayesian_data <- power_results %>%
  mutate(
    prior_h1 = 0.5,
    p_nonsig_h1 = 1 - power,
    p_nonsig_h0 = 0.95,
    posterior_h0 = (p_nonsig_h0 * (1 - prior_h1)) /
      (p_nonsig_h1 * prior_h1 + p_nonsig_h0 * (1 - prior_h1)),
    posterior_h1 = 1 - posterior_h0,
    n_label = factor(paste0("n=", n_donors),
                     levels = paste0("n=", N_DONORS)),
    icc_facet = factor(paste0("True ICC = ", true_icc * 100, "%"),
                       levels = paste0("True ICC = ", TRUE_ICCS * 100, "%"))
  )

p7_data <- bayesian_data %>%
  select(n_label, n_donors, icc_facet, prior_h1, posterior_h1) %>%
  pivot_longer(cols = c(prior_h1, posterior_h1),
               names_to = "belief", values_to = "probability") %>%
  mutate(belief = ifelse(
    belief == "prior_h1",
    "Prior: P(genetic effects exist)",
    "Posterior: P(genetic effects | non-sig. test)"
  ))

p7 <- ggplot(p7_data,
             aes(x = n_label, y = probability, fill = belief)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  facet_wrap(~icc_facet, ncol = 5) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_fill_manual(values = c("#E63946", "#264653")) +
  labs(
    title = "A non-significant donor test barely changes prior beliefs",
    x = "Number of donors",
    y = "Probability",
    fill = NULL
  ) +
  theme_pub +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.text.x = element_text(size=9))


tiff("../pics/fig7_bayesian_ICC.tiff", units = "in", res = 300, width = 7, height = 5)
print(p7)
dev.off()


# ---------- FIGURE 8: Summary heatmap ----------

heatmap_data <- power_results %>%
  mutate(
    n_label = factor(paste0("n=", n_donors),
                     levels = rev(paste0("n=", N_DONORS))),
    icc_label = factor(paste0(true_icc * 100, "%"),
                       levels = paste0(TRUE_ICCS * 100, "%")),
    power_label = sprintf("%.0f%%", power * 100)
  )

p8 <- ggplot(heatmap_data, aes(x = icc_label, y = n_label, fill = power)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = power_label),
            color = ifelse(heatmap_data$power > 0.5, "white", "grey20"),
            size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#FEE2E2", mid = "#F4A261", high = "#264653",
    midpoint = 0.5,
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  labs(
    title = "Statistical power to detect donor-level variance",
    x = "True genetic contribution (ICC)",
    y = "Number of donors",
    fill = "Power"
  ) +
  theme_pub +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  )



tiff("../pics/fig8_power_heatmap.tiff", units = "in", res = 300, width = 6, height = 5)
print(p8)
dev.off()


################################################################################
# SUMMARY STATISTICS FOR MANUSCRIPT TEXT
################################################################################

cat("\n\n")
cat("################################################################\n")
cat("# SUMMARY STATISTICS                    #\n")
cat("################################################################\n\n")

for (icc_val in TRUE_ICCS) {
  cat(sprintf("--- True ICC = %.0f%% ---\n", icc_val * 100))
  subset <- power_results %>% filter(true_icc == icc_val)
  for (i in 1:nrow(subset)) {
    row <- subset[i, ]
    cat(sprintf("n=%2d: Power=%5.1f%% | Mean ICC=%5.1f%% | 95%% range=[%.1f%%, %.1f%%] | ICC=0 in %.1f%%",
                row$n_donors, row$power * 100, row$mean_est_icc * 100,
                row$q025_icc * 100, row$q975_icc * 100, row$pct_icc_zero * 100))
    # Add Type M if available
    tm <- type_m_results %>%
      filter(true_icc == icc_val, n_donors == row$n_donors)
    if (nrow(tm) > 0) {
      cat(sprintf(" | Type M=%.1f\u00D7", tm$type_m_median[1]))
    }
    cat("\n")
  }
  cat("\n")
}

cat("\n--- Batch-Donor Confounding (true donor ICC=20%, batch ICC=15%) ---\n")
for (i in 1:nrow(confound_results)) {
  row <- confound_results[i, ]
  cat(sprintf("n=%2d: Na\u00EFve ICC=%.1f%% [%.1f\u2013%.1f%%] vs Nested donor ICC=%.1f%% [%.1f\u2013%.1f%%] (donor=0 in %.1f%%)\n",
              row$n_donors,
              row$mean_icc_oneway * 100,
              row$q025_icc_oneway * 100,
              row$q975_icc_oneway * 100,
              row$mean_icc_donor_nested * 100,
              row$q025_icc_donor_nested * 100,
              row$q975_icc_donor_nested * 100,
              row$pct_donor_zero_nested * 100))
}
