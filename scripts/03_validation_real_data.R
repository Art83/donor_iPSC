################################################################################
# Empirical Resampling: HipSci iPSC Gene expression data
#
# Uses processed RNA-seq data from Mirauta et al., eLife, 2020
# (Figure 1 source data 1 & 2): 202 iPSC lines from 151 donors.
#
# Strategy:
#   1. Map lines to donors using metadata
#   2. Compute gene-wise ICC from the 51 donors with two lines = reference
#   3. Check how much of that reference is shared batch rather than donor
#      (STEP 3b: negative control and within-batch ICC)
#   4. Repeatedly subsample to n = 2, 3, 5, 10, 20, 50 donors
#   5. Re-estimate gene-wise ICC from subsample
#   6. Compare distribution of subsampled estimates to the reference
#
# This is the "real data" validation
#
# Subsampling pool. A donor contributes to the ICC only if it has two lines.
# Only 51 of the 151 donors do. Drawing from all 151 means most small draws
# hold one replicated donor or none, so the estimate rests on fewer donors
# than the stated n and the result is specific to this dataset's missing
# replicates. SUBSAMPLE_POOL = "replicated" draws from the 51, so every sampled
# donor has two lines, as in a designed study. "all" makes it draw feom 151.
#
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

set.seed(42)

META_FILE <- "elife_data/elife-57390-fig1-data1-v3.xls"   # sample metadata
EXPR_FILE <- "elife_data/elife-57390-fig1-data2-v3.xls"   # expression matrix

SUBSAMPLE_POOL <- "replicated"   # "replicated" (51 two-line donors) or "all" (151)
N_CONTROL_PAIRINGS <- 20         # random pairings for the batch negative control

meta <- read.delim(META_FILE, row.names = 1, check.names = FALSE,
                   stringsAsFactors = FALSE)
expr <- read.delim(EXPR_FILE, row.names = 1, check.names = FALSE,
                   stringsAsFactors = FALSE)

cat(sprintf("Metadata: %d lines\n", nrow(meta)))
cat(sprintf("Expression matrix: %d genes x %d lines\n",
            nrow(expr), ncol(expr)))

line_names <- colnames(expr)

meta$line <- meta[["Public name"]]
meta$donor <- meta[["Donor"]]

# Match expression columns to donors (batch is the batch column of the source
# metadata, the only processing batch it records)
line_donor_map <- meta %>%
  select(line, donor, batch) %>%
  filter(line %in% line_names) %>%
  distinct()


missing_lines <- setdiff(line_names, line_donor_map$line)
if (length(missing_lines) > 0) {
  cat(sprintf("%d lines not in metadata, extracting donor from name\n",
              length(missing_lines)))
  extra_map <- data.frame(
    line = missing_lines,
    donor = sub(".*-([a-z]+)_.*", "\\1", missing_lines)
  )
  line_donor_map <- bind_rows(line_donor_map, extra_map)
}

# Final mapping
line_donor_map <- line_donor_map %>% filter(line %in% line_names)
n_lines <- nrow(line_donor_map)
n_donors <- length(unique(line_donor_map$donor))

cat(sprintf("Matched: %d lines from %d unique donors\n", n_lines, n_donors))


lines_per_donor <- line_donor_map %>%
  group_by(donor) %>%
  summarize(n_lines = n(), .groups = "drop")
cat(sprintf("Lines per donor: median = %d, range = %d-%d\n",
            median(lines_per_donor$n_lines),
            min(lines_per_donor$n_lines),
            max(lines_per_donor$n_lines)))


# ==============================================================================
# STEP 2: Average within donor (for donors with multiple lines)
# ==============================================================================

cat("\n=== Averaging within donors ===\n")

expr_mat <- as.matrix(expr[, line_donor_map$line])

# Create donor-level expression matrix
donors <- unique(line_donor_map$donor)
expr_by_donor <- matrix(NA, nrow = nrow(expr_mat), ncol = length(donors))
rownames(expr_by_donor) <- rownames(expr_mat)
colnames(expr_by_donor) <- donors

for (d in donors) {
  d_lines <- line_donor_map$line[line_donor_map$donor == d]
  if (length(d_lines) == 1) {
    expr_by_donor[, d] <- expr_mat[, d_lines]
  } else {
    expr_by_donor[, d] <- rowMeans(expr_mat[, d_lines], na.rm = TRUE)
  }
}

cat(sprintf("Donor-level matrix: %d genes x %d donors\n",
            nrow(expr_by_donor), ncol(expr_by_donor)))


# ==============================================================================
# STEP 3: Compute gene-wise ICC from full dataset
#
# For the full dataset, ICC is computed from the LINE-level data
# (not donor-averaged), treating donor as the grouping factor.
# This uses the same ANOVA-based ICC as the simulation.
#
# Only include donors with >= 2 lines for ICC computation on the full set.
# For the subsampling, we use the donor-averaged profiles and compute
# between-donor variance as a proportion of total variance.
# ==============================================================================

cat("\n=== Computing full-dataset gene-wise ICC ===\n")

# Method: for each gene, one-way ANOVA with donor as factor
# Using all lines (not averaged), restricted to donors with >= 2 lines
multi_line_donors <- lines_per_donor$donor[lines_per_donor$n_lines >= 2]
cat(sprintf("Donors with >= 2 lines: %d\n", length(multi_line_donors)))

multi_lines <- line_donor_map %>% filter(donor %in% multi_line_donors)
expr_multi <- expr_mat[, multi_lines$line]

# Compute ICC for each gene (vectorised where possible)
n_genes <- nrow(expr_multi)

compute_gene_icc <- function(expr_matrix, donor_ids) {
  # Vectorised one-way random-effects ICC, all genes at once.
  # expr_matrix: (n_genes x n_samples)
  # donor_ids:length-n_samples vector giving donor per column
  # Returns length-n_genes vector of ICC estimates, truncated at 0.
  # Returns all-NA if N == a (no within-donor df) or a < 2.
  # Assumes no NAs in expr_matrix.
  stopifnot(is.matrix(expr_matrix), ncol(expr_matrix) == length(donor_ids))
  donors  <- unique(donor_ids)
  a <- length(donors)
  N <- length(donor_ids)
  n_genes <- nrow(expr_matrix)
  
  if (a < 2 || N == a) return(rep(NA_real_, n_genes))
  
  grp_sizes <- as.integer(table(factor(donor_ids, levels = donors)))
  kbar <- mean(grp_sizes)
  
  # Indicator matrix (N x a): one 1 per row
  ind <- matrix(0, nrow = N, ncol = a)
  ind[cbind(seq_len(N), match(donor_ids, donors))] <- 1L
  
  # Donor means per gene
  grp_means <- sweep(expr_matrix %*% ind, 2, grp_sizes, "/")
  grand_means <- rowMeans(expr_matrix)
  
  dev_between <- sweep(grp_means, 1, grand_means, "-")
  ss_between <- rowSums(sweep(dev_between^2, 2, grp_sizes, "*"))
  
  idx <- match(donor_ids, donors)
  expanded <- grp_means[, idx, drop = FALSE]
  ss_within <- rowSums((expr_matrix - expanded)^2)
  
  ms_between <- ss_between / (a - 1)
  ms_within <- ss_within  / (N - a)
  
  icc <- (ms_between - ms_within) / (ms_between + (kbar - 1) * ms_within)
  pmax(icc, 0)
}

full_icc <- compute_gene_icc(expr_multi, multi_lines$donor)
names(full_icc) <- rownames(expr_multi)


full_icc <- full_icc[!is.na(full_icc)]
cat(sprintf("ICC computed for %d genes\n", length(full_icc)))
cat(sprintf("Median ICC = %.1f%%, Mean ICC = %.1f%%\n",
            median(full_icc) * 100, mean(full_icc) * 100))
cat(sprintf("ICC > 20%%: %d genes (%.1f%%)\n",
            sum(full_icc > 0.20), mean(full_icc > 0.20) * 100))
cat(sprintf("ICC > 50%%: %d genes (%.1f%%)\n",
            sum(full_icc > 0.50), mean(full_icc > 0.50) * 100))


# ==============================================================================
# STEP 3b: Batch structure of the replicate lines
#
# The reference ICC counts everything two lines of a donor share as donor-level
# variance. If both lines were processed in the same batch, shared batch is
# counted as donor. Three checks:
#   1   how many replicated donors have both lines in one batch
#   2  negative control: pair unrelated lines that share a batch and treat
#         each pair as a "donor". Its ICC, compared with random unrelated
#         pairs, is variance due to batch alone.
#   3 within-batch donor ICC: donors nested in batch, with the between-donor
#         mean square taken around the batch mean so that batch drops out.
# ==============================================================================

cat("\n=== Batch structure of replicate lines ===\n")

# do the lines of a donor share a batch
batches_per_donor <- multi_lines %>%
  group_by(donor) %>%
  summarize(n_batches = n_distinct(batch), .groups = "drop")
same_batch_donors <- batches_per_donor$donor[batches_per_donor$n_batches == 1]
cat(sprintf("Replicated donors with all lines in one batch: %d of %d\n",
            length(same_batch_donors), length(multi_line_donors)))

# negative control
pair_unrelated_within_batch <- function(map) {
  # Greedy random pairing of lines from different donors inside each batch
  pairs <- list()
  for (b in unique(map$batch)) {
    idx  <- which(map$batch == b)
    if (length(idx) > 1) idx <- sample(idx)
    used <- rep(FALSE, length(idx))
    pid  <- 0
    for (x in seq_along(idx)) {
      if (used[x]) next
      for (y in seq_along(idx)) {
        if (y <= x || used[y]) next
        if (map$donor[idx[x]] != map$donor[idx[y]]) {
          pid <- pid + 1
          pairs[[length(pairs) + 1]] <- data.frame(
            line = map$line[idx[c(x, y)]],
            pseudo_donor = paste(b, pid, sep = "_"))
          used[c(x, y)] <- TRUE
          break
        }
      }
    }
  }
  do.call(rbind, pairs)
}

pair_unrelated_at_random <- function(map, n_pairs) {
  # One line from each of 2 * n_pairs different donors, paired at random
  d <- sample(unique(map$donor), 2 * n_pairs)
  l <- sapply(d, function(x) {
    cand <- map$line[map$donor == x]
    if (length(cand) > 1) sample(cand, 1) else cand
  })
  data.frame(line = l, pseudo_donor = rep(seq_len(n_pairs), each = 2))
}

# Both controls use as many pairs as there are replicated donors, so their
# sampling noise matches that of the reference.
n_control_pairs <- length(multi_line_donors)
control_same_batch <- matrix(NA_real_, nrow(expr_mat), N_CONTROL_PAIRINGS)
control_random <- matrix(NA_real_, nrow(expr_mat), N_CONTROL_PAIRINGS)
for (i in seq_len(N_CONTROL_PAIRINGS)) {
  p1 <- pair_unrelated_within_batch(line_donor_map)
  keep <- sample(unique(p1$pseudo_donor), n_control_pairs)
  p1 <- p1[p1$pseudo_donor %in% keep, ]
  control_same_batch[, i] <- compute_gene_icc(expr_mat[, p1$line], p1$pseudo_donor)
  p2 <- pair_unrelated_at_random(line_donor_map, n_control_pairs)
  control_random[, i] <- compute_gene_icc(expr_mat[, p2$line], p2$pseudo_donor)
}
# median over genes and pairings, with the range of per-pairing medians
control_same_batch_median <- apply(control_same_batch, 2, median)
control_random_median <- apply(control_random, 2, median)

cat(sprintf("Negative control, unrelated lines sharing a batch (%d pairs, %d pairings): median ICC = %.1f%% [%.1f%%-%.1f%% across pairings]\n",
            n_control_pairs, N_CONTROL_PAIRINGS,
            median(control_same_batch) * 100,
            min(control_same_batch_median) * 100, max(control_same_batch_median) * 100))
cat(sprintf("Negative control, unrelated lines paired at random (%d pairs, %d pairings): median ICC = %.1f%% [%.1f%%-%.1f%% across pairings]\n",
            n_control_pairs, N_CONTROL_PAIRINGS,
            median(control_random) * 100,
            min(control_random_median) * 100, max(control_random_median) * 100))

# within-batch donor ICC
compute_gene_icc_within_batch <- function(expr_matrix, donor_ids, batch_ids) {
  # One-way random-effects ICC with donors nested in batch, all genes at once.
  # Every line of a donor must sit in one batch. Truncated at 0.
  d <- factor(donor_ids)
  fb <- factor(batch_ids)
  n_genes <- nrow(expr_matrix)
  n <- as.numeric(table(d))
  a <- nlevels(d)
  N <- length(d)
  stopifnot(all(tapply(as.character(batch_ids), d,
                       function(x) length(unique(x))) == 1))
  batch_of_donor <- factor(tapply(as.character(batch_ids), d, function(x) x[1]),
                           levels = levels(fb))
  nb <- nlevels(fb)
  if (a - nb < 1 || N == a) return(rep(NA_real_, n_genes))
  
  donor_means <- t(rowsum(t(expr_matrix), d)) / rep(n, each = n_genes)
  batch_means <- t(rowsum(t(expr_matrix), fb)) /
    rep(as.numeric(table(fb)), each = n_genes)
  
  ss_within <- rowSums((expr_matrix - donor_means[, as.integer(d), drop = FALSE])^2)
  dev <- donor_means - batch_means[, as.integer(batch_of_donor), drop = FALSE]
  ss_between <- rowSums(sweep(dev^2, 2, n, "*"))
  
  ms_within <- ss_within / (N - a)
  ms_between <- ss_between / (a - nb)
  
  # coefficient of the donor variance in E[MS_between] for the nested design
  k0 <- (N - sum(tapply(n^2, batch_of_donor, sum) /
                   tapply(n, batch_of_donor, sum))) / (a - nb)
  
  icc <- (ms_between - ms_within) / (ms_between + (k0 - 1) * ms_within)
  pmax(icc, 0)
}

# donors with all lines in one batch, in batches holding at least two of them
nested_lines <- multi_lines %>% filter(donor %in% same_batch_donors)
donors_per_batch <- nested_lines %>%
  group_by(batch) %>%
  summarize(n_d = n_distinct(donor), .groups = "drop")
nested_lines <- nested_lines %>%
  filter(batch %in% donors_per_batch$batch[donors_per_batch$n_d >= 2])

within_batch_icc <- compute_gene_icc_within_batch(
  expr_mat[, nested_lines$line], nested_lines$donor, nested_lines$batch)
names(within_batch_icc) <- rownames(expr_mat)
same_donors_icc <- compute_gene_icc(expr_mat[, nested_lines$line], nested_lines$donor)

cat(sprintf("Within-batch donor ICC (%d donors in %d batches): median = %.1f%%\n",
            n_distinct(nested_lines$donor), n_distinct(nested_lines$batch),
            median(within_batch_icc) * 100))
cat(sprintf("Same %d donors with batch ignored: median = %.1f%%\n",
            n_distinct(nested_lines$donor), median(same_donors_icc) * 100))
cat(sprintf("Share of the reference attributable to batch: %.1f percentage points\n",
            (median(same_donors_icc) - median(within_batch_icc)) * 100))

# ==============================================================================
#  STEP 4: Subsampling: repeatedly draw n donors, re-estimate ICC
# ==============================================================================

cat("\n=== Subsampling analysis ===\n")

subsample_pool <- if (SUBSAMPLE_POOL == "replicated") multi_line_donors else donors
cat(sprintf("Subsampling pool: %s (%d donors)\n", SUBSAMPLE_POOL, length(subsample_pool)))

n_subsample_vec <- c(2, 3, 5, 10, 20, 50)
n_iterations <- 1000  # per subsample size 

top_genes <- names(full_icc)
cat(sprintf("Using all %d genes for subsampling\n", length(top_genes)))

# Track subsamples that cannot be estimated (all sampled donors are single-line)
n_rejected <- integer(length(n_subsample_vec))
names(n_rejected) <- as.character(n_subsample_vec)

subsample_results <- list()

for (n_sub in n_subsample_vec) {
  if (n_sub > length(subsample_pool)) next
  
  cat(sprintf("n = %d donors: ", n_sub))
  flush.console()
  
  # Store median ICC per iteration (NA for rejected draws)
  iter_median_icc <- rep(NA_real_, n_iterations)
  # Store per-gene ICCs for each iteration
  all_iter_iccs <- matrix(NA_real_, nrow = length(top_genes), ncol = n_iterations)
  
  for (iter in 1:n_iterations) {
    sampled_donors <- sample(subsample_pool, n_sub)
    sampled_lines <- line_donor_map %>%
      filter(donor %in% sampled_donors)
    
    # Reject if no within-donor variance (all sampled donors are single-line)
    if (nrow(sampled_lines) <= n_sub) {
      n_rejected[as.character(n_sub)] <- n_rejected[as.character(n_sub)] + 1L
      next
    }
    
    sub_expr <- expr_mat[top_genes, sampled_lines$line, drop = FALSE]
    sub_donors <- sampled_lines$donor
    
    sub_icc <- compute_gene_icc(sub_expr, sub_donors)
    
    all_iter_iccs[, iter] <- sub_icc
    iter_median_icc[iter] <- median(sub_icc, na.rm = TRUE)
  }
  
  cat(sprintf("median of medians = %.1f%%\n",
              median(iter_median_icc, na.rm = TRUE) * 100))
  
  subsample_results[[as.character(n_sub)]] <- list(
    n_sub = n_sub,
    iter_medians = iter_median_icc,
    all_iccs = all_iter_iccs,
    gene_names = top_genes
  )
}


# ==============================================================================
# STEP 5: Figures
# ==============================================================================

cat("\n=== Generating figures ===\n")

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

# ---------- FIGURE D1: Full-dataset ICC distribution ----------

full_icc_df <- data.frame(icc = full_icc)

pD1 <- ggplot(full_icc_df, aes(x = icc)) +
  geom_histogram(bins = 80, fill = "#264653", color = "white",
                 linewidth = 0.2, alpha = 0.85) +
  geom_vline(xintercept = median(full_icc), linetype = "dashed",
             color = "#E63946", linewidth = 0.8) +
  annotate("text", x = median(full_icc) + 0.03, y = Inf,
           label = sprintf("Median ICC = %.0f%%", median(full_icc) * 100),
           color = "#E63946", vjust = 2, hjust = 0, size = 4,
           fontface = "bold") +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title = "Distribution of donor-level ICC across genes",
    x = "Donor ICC (fraction of variance explained by donor)",
    y = "Number of genes"
  ) +
  theme_pub + 
  theme(plot.title = element_text(size = 14))


tiff("../pics/figD_emp_median.tiff", units = "in", res = 300, width = 6, height = 5)
print(pD1)
dev.off()



# ---------- FIGURE D2: Subsampled ICC vs full ICC ----------
# For each gene, compare full-dataset ICC to subsampled ICC

# Pick genes with known high ICC and known low ICC
high_icc_genes <- names(sort(full_icc[top_genes], decreasing = TRUE))[1:500]
low_icc_genes <- names(sort(full_icc[top_genes]))[1:500]

scatter_data <- data.frame()

for (n_sub_char in names(subsample_results)) {
  res <- subsample_results[[n_sub_char]]
  # Average ICC across iterations for each gene
  mean_sub_icc <- rowMeans(res$all_iccs, na.rm = TRUE)
  names(mean_sub_icc) <- res$gene_names
  
  scatter_data <- rbind(scatter_data, data.frame(
    gene = res$gene_names,
    full_icc = full_icc[res$gene_names],
    sub_icc = mean_sub_icc,
    n_sub = res$n_sub,
    n_label = factor(paste0("n=", res$n_sub),
                     levels = paste0("n=", n_subsample_vec))
  ))
}

n3 <- scatter_data[scatter_data$n_sub == 3, ]
loess_n3 <- loess(sub_icc ~ full_icc, data = n3, span = 0.75) 
predict(loess_n3, newdata = data.frame(full_icc = c(0.05, 0.10, 0.50, 0.70, 0.95)))

# Genes with reference ICC near 10%
window_10 <- abs(n3$full_icc - 0.10) < 0.02
quantile(n3$sub_icc[window_10], c(0.025, 0.975), na.rm = TRUE)

window_70 <- abs(n3$full_icc - 0.70) < 0.02
quantile(n3$sub_icc[window_70], c(0.025, 0.975), na.rm = TRUE)


n3_all <- subsample_results[["3"]]$all_iccs  # genes x iterations
gene_names <- subsample_results[["3"]]$gene_names
ref <- full_icc[gene_names]

# For genes with reference ICC near 10%: distribution across all (gene, iteration) pairs
idx_10 <- which(abs(ref - 0.10) < 0.02)
vals_10 <- as.vector(n3_all[idx_10, ])
quantile(vals_10, c(0.025, 0.975), na.rm = TRUE)


# Flatten per-iteration ICCs with their reference values
ref_long <- rep(ref, ncol(n3_all))
sub_long <- as.vector(n3_all)
ok <- !is.na(sub_long)

# Genes for which some iteration produced a subsampled ICC near 60%
window_60 <- abs(sub_long - 0.60) < 0.02 & ok
quantile(ref_long[window_60], c(0.025, 0.975))


lo <- predict(loess_n3, newdata = data.frame(full_icc = c(0.1, 0.9)))
slope_n3 <- (lo[2] - lo[1]) / 0.8

pD2 <- ggplot(scatter_data %>% filter(n_sub %in% c(3, 5, 10)),
              aes(x = full_icc, y = sub_icc)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#E63946", linewidth = 0.5) +
  geom_point(alpha = 0.1, size = 0.5, color = "#264653") +
  geom_smooth(method = "loess", color = "#F4A261", linewidth = 1,
              se = FALSE) +
  facet_wrap(~n_label, ncol = 4) +
  scale_x_continuous(labels = c(0,25,50,75,100)) +
  scale_y_continuous(labels = c(0,25,50,75,100)) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Subsampled ICC vs reference ICC per gene",
    x = sprintf("Reference ICC (%d donors with two lines)", length(multi_line_donors)),
    y = "Mean ICC"
  ) +
  theme_pub +
  theme(legend.position = "none")


tiff("pics/figD_subsamples_new.tiff", units = "in", res = 300, width = 6, height = 5)
print(pD2)
dev.off()



# ---------- FIGURE D3: Distribution of subsampled median ICCs ----------

median_data <- data.frame()
for (n_sub_char in names(subsample_results)) {
  res <- subsample_results[[n_sub_char]]
  median_data <- rbind(median_data, data.frame(
    n_sub = res$n_sub,
    median_icc = res$iter_medians,
    n_label = factor(paste0("n=", res$n_sub),
                     levels = paste0("n=", n_subsample_vec))
  ))
}

full_median <- median(full_icc)

pD3 <- ggplot(median_data %>% filter(n_label != "n=50"), aes(x = median_icc, fill = n_label)) +
  geom_histogram(bins = 40, alpha = 0.85, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = full_median, linetype = "dashed",
             color = "#E63946", linewidth = 0.8) +
  facet_wrap(~n_label, ncol = 3, scales = "free_y") +
  scale_x_continuous(labels = c(0,25,50,75,100)) +
  scale_fill_manual(values = c("#E63946", "#F4A261", "#2A9D8F",
                               "#264653", "#6A4C93", "#1982C4")) +
  labs(
    title = "Distribution of genome-wide median ICC\n from subsampled donors",
    x = "Median ICC across genes (from subsampled donors)",
    y = "Count"
  ) +
  theme_pub +
  theme(legend.position = "none",
        plot.title = element_text(size = 14))


tiff("pics/figD_median_subsamples.tiff", units = "in", res = 300, width = 6, height = 5)
print(pD3)
dev.off()



# ---------- FIGURE D4: Misclassification rates ----------
# How often does subsampling lead to wrong conclusions about a gene?
# Define: "high ICC" = full ICC > 20%, "low ICC" = full ICC < 10%
# At each n, what fraction of high-ICC genes are estimated as <10%?
# And what fraction of low-ICC genes are estimated as >20%?

misclass_data <- data.frame()

for (n_sub_char in names(subsample_results)) {
  res <- subsample_results[[n_sub_char]]
  
  for (iter in 1:ncol(res$all_iccs)) {
    sub_iccs <- res$all_iccs[, iter]
    names(sub_iccs) <- res$gene_names
    true_iccs <- full_icc[res$gene_names]
    
    # High ICC genes misclassified as low
    truly_high <- true_iccs > 0.20
    called_low <- sub_iccs < 0.10
    miss_high <- mean(called_low[truly_high], na.rm = TRUE)
    
    # Low ICC genes misclassified as high
    truly_low <- true_iccs < 0.10
    called_high <- sub_iccs > 0.20
    miss_low <- mean(called_high[truly_low], na.rm = TRUE)
    
    misclass_data <- rbind(misclass_data, data.frame(
      n_sub = res$n_sub,
      miss_high_as_low = miss_high,
      miss_low_as_high = miss_low
    ))
  }
}

misclass_summary <- misclass_data %>%
  group_by(n_sub) %>%
  summarize(
    high_as_low_mean = mean(miss_high_as_low, na.rm = TRUE),
    low_as_high_mean = mean(miss_low_as_high, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(high_as_low_mean, low_as_high_mean),
               names_to = "error_type", values_to = "rate") %>%
  mutate(
    error_label = ifelse(
      error_type == "high_as_low_mean",
      "True ICC>20% called <10%\n(miss real genetic effect)",
      "True ICC<10% called >20%\n(false genetic effect)"
    ),
    n_label = factor(paste0("n=", n_sub),
                     levels = paste0("n=", n_subsample_vec))
  )

pD4 <- ggplot(misclass_summary, aes(x = n_label, y = rate, fill = error_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#E63946", "#264653")) +
  labs(
    title = "Gene-level misclassification rates from subsampled ICC",
#    subtitle = sprintf(
#      "At n=3, %.0f%% of genes with reference ICC above 20%% are estimated below 10%%",
#      100 * mean(misclass_data$miss_high_as_low[misclass_data$n_sub == 3], na.rm = TRUE)),
    x = "Number of donors sampled",
    y = "Misclassification rate",
    fill = NULL
  ) +
  theme_pub +
  theme(legend.position = "top")

batch_fig_levels <- c(
  "Reference, lines grouped by donor",
  "Donor within batch",
  "Unrelated lines sharing a batch",
  "Unrelated lines paired at random"
)



# The two controls pool all pairings, weighted so that counts are per pairing
batch_fig_data <- rbind(
  data.frame(icc = as.numeric(full_icc), grouping = batch_fig_levels[1], w = 1),
  data.frame(icc = as.numeric(within_batch_icc), grouping = batch_fig_levels[2], w = 1),
  data.frame(icc = as.vector(control_same_batch), grouping = batch_fig_levels[3],
             w = 1 / N_CONTROL_PAIRINGS),
  data.frame(icc = as.vector(control_random), grouping = batch_fig_levels[4],
             w = 1 / N_CONTROL_PAIRINGS)
) %>%
  mutate(grouping = factor(grouping, levels = batch_fig_levels))

batch_fig_medians <- batch_fig_data %>%
  group_by(grouping) %>%
  summarize(med = median(icc), .groups = "drop")

pD5 <- ggplot(batch_fig_data, aes(x = icc, weight = w)) +
  geom_histogram(binwidth = 1 / 60, boundary = 0, fill = "#264653", color = "white",
                 linewidth = 0.2, alpha = 0.85) +
  geom_vline(data = batch_fig_medians, aes(xintercept = med),
             linetype = "dashed", color = "#E63946", linewidth = 0.8) +
  geom_text(data = batch_fig_medians, inherit.aes = FALSE,
            aes(x = med + 0.03, y = Inf, label = sprintf("Median = %.0f%%", med * 100)),
            color = "#E63946", vjust = 2, hjust = 0, size = 4, fontface = "bold") +
  facet_wrap(~grouping, ncol = 1, scales = "free_y") +
  scale_x_continuous(labels = scales::percent) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "How much of the donor ICC is shared batch",
    x = "Gene-wise ICC",
    y = "Number of genes"
  ) +
  theme_pub +
  theme(plot.title = element_text(size = 14))

tiff("pics/figD_shared_batch.tiff", units = "in", res = 300, width = 6, height = 6)
print(pD5)
dev.off()

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n\n")
cat("################################################################\n")
cat("# EMPIRICAL RESAMPLING SUMMARY                                #\n")
cat("################################################################\n\n")

cat(sprintf("Dataset: %d iPSC lines, %d donors\n", n_lines, n_donors))
cat(sprintf("Genes with ICC computed: %d\n", length(full_icc)))
cat(sprintf("Median genome-wide ICC: %.1f%%\n", median(full_icc) * 100))
cat(sprintf("Mean genome-wide ICC: %.1f%%\n\n", mean(full_icc) * 100))

cat("Subsampling results (median of median ICCs across iterations):\n")
for (n_sub_char in names(subsample_results)) {
  res <- subsample_results[[n_sub_char]]
  rej <- n_rejected[n_sub_char]
  medians <- res$iter_medians[!is.na(res$iter_medians)]
  cat(sprintf("n=%2d: median ICC = %.1f%% [%.1f%%–%.1f%%] | rejected (all single-line): %d/%d (%.0f%%)\n",
              res$n_sub,
              median(medians) * 100,
              quantile(medians, 0.025) * 100,
              quantile(medians, 0.975) * 100,
              rej, n_iterations, 100 * rej / n_iterations))
}

cat("\nAgreement of subsampled with reference ICC, per gene:\n")
for (n_sub_char in names(subsample_results)) {
  res <- subsample_results[[n_sub_char]]
  d_n <- scatter_data[scatter_data$n_sub == res$n_sub, ]
  lin <- coef(lm(sub_icc ~ full_icc, data = d_n))[2]
  lo_n <- predict(loess(sub_icc ~ full_icc, data = d_n, span = 0.75),
                  newdata = data.frame(full_icc = c(0.1, 0.9)))
  cat(sprintf("n=%2d: slope %.2f (linear), %.2f (loess, reference 10%% to 90%%) | estimates at zero: %.0f%%\n",
              res$n_sub, lin, (lo_n[2] - lo_n[1]) / 0.8,
              100 * mean(res$all_iccs == 0, na.rm = TRUE)))
}
if (SUBSAMPLE_POOL == "replicated") {
  cat(sprintf("Note: the pool holds %d donors, so draws of 20 and 50 overlap heavily and\n",
              length(subsample_pool)))
  cat("their intervals are narrower than they would be for independent samples.\n")
}

cat("\nMisclassification rates:\n")
mc <- misclass_data %>%
  group_by(n_sub) %>%
  summarize(
    miss_high = mean(miss_high_as_low, na.rm = TRUE),
    miss_low = mean(miss_low_as_high, na.rm = TRUE),
    .groups = "drop")
for (i in 1:nrow(mc)) {
  cat(sprintf("n=%2d: %.0f%% of high-ICC genes called low | %.0f%% of low-ICC genes called high\n",
              mc$n_sub[i], mc$miss_high[i] * 100, mc$miss_low[i] * 100))
}
