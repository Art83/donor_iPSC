# ==============================================================================
#
# Effects against the variance each contrast faces, and the N it needs.
#
# Every contrast is fitted per protein as a linear model with the covariates a
# lab would use. The effect is standardised by the residual SD of that model:
#   imposed in one genome (mouse allele, diet, organoid time in culture)
#       residual SD = between replicates of that genome
#   donor attributes in iAstrocyte lines (APOE4, pathology)
#       residual SD = between donor lines
#   human brain, Diverse Cohorts frontal outside ROSMAP and Banner
#       residual SD = between donors
#
# Effects from small groups are inflated by noise, so each contrast is refitted
# with its label shuffled within strata and reported against that null. All
# relabellings are enumerated when there are MAX_EXACT or fewer, otherwise
# N_PERM random ones are drawn.
#
# N per group is the two-sample normal approximation for 80% power at two-sided
# alpha 0.05, n = 2 (z_0.975 + z_0.8)^2 / d^2, at the median and 90th percentile
# |d| with the same quantile of the null removed as a difference of squares,
# d^2 = observed^2 - null^2. No N is given when the observed quantile is not
# above the shuffled-label distribution at p <= 0.05.
#
# Each model contrast is also correlated with the human AD effect on shared
# proteins, against the same null.
# ==============================================================================

set.seed(20261026)

CACHE <- file.path("data", "prep_systems.rds")
TABLES <- "tables"

N_PERM <- 200
MAX_EXACT <- 2000
MIN_PRESENT <- 1.0      # proteins complete in the samples of a contrast
Z <- qnorm(0.975) + qnorm(0.80)
FC_CUTS <- c(1.1, 1.2, 1.5)

P <- readRDS(CACHE)
S <- P$samples


# ==============================================================================
# Helpers
# ==============================================================================

pinv <- function(A) {
  s <- svd(A)
  pos <- s$d > max(dim(A)) * max(s$d) * .Machine$double.eps
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}

design <- function(x, cov) {
  X <- cbind(1, x)
  for (v in names(cov)) {
    f <- factor(as.character(cov[[v]]))
    if (nlevels(f) > 1) X <- cbind(X, model.matrix(~f)[, -1, drop = FALSE])
  }
  unname(X)
}

fit <- function(Y, x, cov) {
  # OLS per protein. Y proteins x samples. Returns the coefficient of x, the
  # residual SD, t and the residual degrees of freedom.
  X <- design(x, cov)
  XtX_inv <- pinv(crossprod(X))
  B <- Y %*% X %*% XtX_inv
  df <- length(x) - qr(X)$rank
  sd <- sqrt(rowSums((Y - B %*% t(X))^2) / max(df, 1))
  list(beta = B[, 2], sd = sd, t = B[, 2] / (sd * sqrt(XtX_inv[2, 2])), df = df)
}

shuffle_within <- function(x, strata) {
  for (s in unique(strata)) {
    i <- which(strata == s)
    if (length(i) > 1) x[i] <- x[sample(i)]
  }
  x
}

relabellings <- function(x, strata) {
  # every within-strata relabelling of a binary x, or NULL above MAX_EXACT
  groups <- lapply(sort(unique(strata)), function(s) which(strata == s))
  ones <- sapply(groups, function(i) sum(x[i]))
  if (prod(choose(lengths(groups), ones)) > MAX_EXACT) return(NULL)
  sets <- mapply(function(i, k) {
    if (k == 0) return(list(integer(0)))
    if (length(i) == 1) return(list(i))
    combn(i, k, simplify = FALSE)
  }, groups, ones, SIMPLIFY = FALSE)
  grid <- expand.grid(lapply(sets, seq_along))
  lapply(seq_len(nrow(grid)), function(r) {
    xs <- numeric(length(x))
    for (g in seq_along(sets)) xs[sets[[g]][[grid[r, g]]]] <- 1
    xs
  })
}

q90 <- function(v) unname(quantile(v, 0.9))

contrast <- function(name, unit, mat, samp, xcol, covcols, strata_cols, binary = TRUE) {
  # samp: one row per sample, rownames are sample IDs
  samp <- samp[complete.cases(samp[, c(xcol, covcols), drop = FALSE]), ]
  samp <- samp[rownames(samp) %in% colnames(mat), ]
  Y <- mat[, rownames(samp), drop = FALSE]
  Y <- Y[rowMeans(!is.na(Y)) >= MIN_PRESENT, , drop = FALSE]
  x <- as.numeric(samp[[xcol]])
  if (!binary) x <- (x - mean(x)) / sd(x)            # effect per SD of the attribute
  cov <- samp[, covcols, drop = FALSE]
  strata <- if (length(strata_cols)) do.call(paste, c(samp[, strata_cols, drop = FALSE], sep = "|"))
            else rep("all", length(x))

  f <- fit(Y, x, cov)
  exact <- if (binary) relabellings(x, strata) else NULL
  perms <- if (!is.null(exact)) exact else replicate(N_PERM, shuffle_within(x, strata), simplify = FALSE)
  null_d  <- matrix(NA_real_, nrow(Y), length(perms), dimnames = list(rownames(Y), NULL))
  null_fc <- matrix(NA_real_, length(perms), 2)
  for (i in seq_along(perms)) {
    fs <- fit(Y, perms[[i]], cov)
    null_d[, i] <- fs$beta / fs$sd
    null_fc[i, ] <- c(median(abs(fs$beta)), q90(abs(fs$beta)))
  }
  n1 <- if (binary) sum(x > 0) else length(x)
  n0 <- if (binary) sum(x <= 0) else 0
  cat(sprintf("  %-40s %-11s n=%3d/%-3d proteins=%5d  df=%d  null: %s %d\n", name, unit, n1, n0,
              nrow(Y), f$df, if (!is.null(exact)) "exact" else "random", length(perms)))
  list(name = name, unit = unit, beta = f$beta, d = f$beta / f$sd, sd = f$sd,
       p = 2 * pt(-abs(f$t), f$df), null_d = null_d, null_fc = null_fc,
       n1 = n1, n0 = n0, exact = !is.null(exact))
}

by_sample <- function(d) {
  rownames(d) <- d$sample
  d
}


# ==============================================================================
# Contrasts
# ==============================================================================

C <- list()
cat("CONTRASTS\n")

# human brain, frontal cortex outside the ROSMAP and Banner cohorts
fm <- P$dc_fr_meta[!P$dc_fr_meta$discovery, ]
fm$ad <- c(AD = 1, Control = 0)[fm$ADoutcome]
fm$e4 <- ifelse(is.na(fm$apoeGenotype), NA, as.numeric(grepl("4", fm$apoeGenotype)))
C[[1]] <- contrast("human: AD vs control", "donor", P$DC_FR, fm, "ad", c("site", "sex"), "site")
C[[2]] <- contrast("human: APOE4 carrier vs not", "donor", P$DC_FR, fm, "e4", c("site", "sex"), "site")
HUMAN_AD <- C[[1]]$d

# iAstrocyte donor attributes
ph <- P$pheno
ph$e4    <- as.numeric(ph$apoe4)
ph$pmAD  <- suppressWarnings(as.numeric(ph$pmAD))
ph$gpath <- suppressWarnings(as.numeric(ph$gpath))
C[[3]] <- contrast("iA: APOE4 carrier vs not", "donor line", P$iA_bc, ph, "e4", c("batch", "sex"), "batch")
C[[4]] <- contrast("iA: pathological AD vs not", "donor line", P$iA_bc, ph, "pmAD", c("batch", "sex"), "batch")
C[[5]] <- contrast("iA: global pathology (per SD)", "donor line", P$iA_bc, ph, "gpath",
                   c("batch", "sex"), "batch", binary = FALSE)

# organoid time in culture, within line
os <- by_sample(S[S$system == "organoid", ])
os$t12 <- as.numeric(os$age == "12wk")
C[[6]] <- contrast("organoid: 12 vs 8 weeks, within line", "organoid", P$O1p, os, "t12", "genome", "genome")

# LOAD2
ls <- by_sample(S[S$system == "LOAD2", ])
ls$band <- as.character(cut(as.numeric(ls$age), c(0, 8, 16, 30),
                            labels = c("3-6 mo", "12-14 mo", "18-19 mo")))
ls <- ls[ls$genotype %in% c("LOAD1", "LOAD2"), ]
ls$hAb <- as.numeric(ls$genotype == "LOAD2")
ls$hfd <- as.numeric(ls$condition == "HFD")
for (b in c("3-6 mo", "12-14 mo", "18-19 mo")) {
  s <- ls[ls$band %in% b & ls$condition == "CD", ]
  C[[length(C) + 1]] <- contrast(sprintf("LOAD2: humanised Abeta, %s", b), "animal", P$M_L2, s,
                                 "hAb", "sex", "sex")
}
C[[length(C) + 1]] <- contrast("LOAD2: high-fat diet, 18 mo", "animal", P$M_L2,
                               ls[ls$band %in% "18-19 mo", ], "hfd", c("genotype", "sex"),
                               c("genotype", "sex"))

# Peng
ps <- by_sample(S[S$system == "Peng", ])
for (m in list(list("FAD", 6), list("FAD", 12), list("NLGF", 18))) {
  arm <- if (m[[1]] == "FAD") "5xFAD" else "APPKI"
  s <- ps[ps$condition == arm & as.numeric(ps$age) == m[[2]] & ps$genotype %in% c(m[[1]], "WT"), ]
  s$tg <- as.numeric(s$genotype == m[[1]])
  cov <- c("batch", if (length(unique(s$sex)) > 1) "sex")
  C[[length(C) + 1]] <- contrast(sprintf("Peng: %s vs WT, %d mo", if (m[[1]] == "FAD") "5xFAD" else "NLGF", m[[2]]),
                                 "animal", P$M_PG, s, "tg", cov, "batch")
}


# ==============================================================================
# Summary
# ==============================================================================

perm_p <- function(obs, null, exact) {
  k <- sum(null >= obs - 1e-12)
  if (exact) k / length(null) else (k + 1) / (length(null) + 1)
}

summarise <- function(c) {
  ad <- abs(c$d)
  null_meds <- apply(abs(c$null_d), 2, median)
  null_q90s <- apply(abs(c$null_d), 2, q90)
  null_med <- mean(null_meds)
  null_q90 <- mean(null_q90s)
  p_med <- perm_p(median(ad), null_meds, c$exact)
  p_q90 <- perm_p(q90(ad), null_q90s, c$exact)
  med_ex <- sqrt(max(median(ad)^2 - null_med^2, 0))
  q90_ex <- sqrt(max(q90(ad)^2 - null_q90^2, 0))
  fc <- abs(c$beta)

  rho <- rho_null95 <- p_conc <- p_min <- NA
  if (!startsWith(c$name, "human")) {
    j <- intersect(names(c$d), names(HUMAN_AD))
    if (length(j) >= 100) {
      rho <- cor(c$d[j], HUMAN_AD[j], method = "spearman")
      r_null <- as.vector(cor(c$null_d[j, , drop = FALSE], HUMAN_AD[j], method = "spearman"))
      k <- sum(abs(r_null) >= abs(rho) - 1e-12)
      n <- length(r_null)
      # under an exact null the observed labelling is one of those enumerated
      p_conc <- if (c$exact) k / n else (k + 1) / (n + 1)
      p_min  <- if (c$exact) 1 / n else 1 / (n + 1)
      rho_null95 <- unname(quantile(abs(r_null), 0.95))
    }
  }
  out <- data.frame(
    contrast = c$name, unit = c$unit, n_group1 = c$n1, n_group0 = c$n0, proteins = length(ad),
    median_abs_d = median(ad), null_median = null_med, q90_abs_d = q90(ad), null_q90 = null_q90,
    n_per_group_median = if (p_med <= 0.05 && med_ex > 0) 2 * Z^2 / med_ex^2 else NA,
    n_per_group_q90    = if (p_q90 <= 0.05 && q90_ex > 0) 2 * Z^2 / q90_ex^2 else NA,
    n_fdr05 = sum(p.adjust(c$p, "BH") < 0.05),
    median_resid_sd = median(c$sd),
    median_fold_change = 2^median(fc), null_median_fold_change = 2^mean(c$null_fc[, 1]),
    q90_fold_change = 2^q90(fc), null_q90_fold_change = 2^mean(c$null_fc[, 2]),
    rho_human_ad = rho, rho_null_95 = rho_null95, p_concordance = p_conc, p_min_attainable = p_min,
    stringsAsFactors = FALSE)
  for (cut in FC_CUTS) out[[sprintf("pct_fold_change_over_%s", cut)]] <- 100 * mean(fc > log2(cut))
  out$d_median_above_null <- med_ex
  out$d_q90_above_null <- q90_ex
  out$p_median_abs_d <- p_med
  out$p_q90_abs_d <- p_q90
  out
}
Sm <- do.call(rbind, lapply(C, summarise))

n_fmt <- function(v) if (is.na(v)) "-" else if (v > 10000) ">10000" else sprintf("%.0f", ceiling(v))

cat("\nSTANDARDISED EFFECTS |d| (effect / residual SD of the contrast), observed against shuffled labels\n")
cat(sprintf("  %-40s %5s %5s | %5s %5s | %7s %7s | %7s\n", "contrast", "med", "null", "q90", "null",
            "N@med", "N@q90", "FDR<.05"))
for (i in seq_len(nrow(Sm))) with(Sm[i, ], cat(sprintf(
  "  %-40s %5.2f %5.2f | %5.2f %5.2f | %7s %7s | %7d\n", contrast, median_abs_d, null_median,
  q90_abs_d, null_q90, n_fmt(n_per_group_median), n_fmt(n_per_group_q90), n_fdr05)))
cat("\n  N@med and N@q90: per group, 80% power at two-sided alpha 0.05 on one protein, using |d| above",
    "\n  the null. The unit is the unit of the contrast (animal, organoid of one line, donor line, donor).",
    "\n  A dash marks a quantile that is not above shuffled labels at p <= 0.05.\n")

cat("\nFOLD CHANGES (absolute, per protein), observed against shuffled labels\n")
cat(sprintf("  %-40s %6s %6s | %6s %6s | %7s %7s %7s\n", "contrast", "median", "null", "q90", "null",
            ">1.1", ">1.2", ">1.5"))
for (i in seq_len(nrow(Sm))) with(Sm[i, ], cat(sprintf(
  "  %-40s %6.3f %6.3f | %6.3f %6.3f | %6.1f%% %6.1f%% %6.1f%%\n", contrast, median_fold_change,
  null_median_fold_change, q90_fold_change, null_q90_fold_change, pct_fold_change_over_1.1,
  pct_fold_change_over_1.2, pct_fold_change_over_1.5)))

cat("\nCONCORDANCE WITH THE HUMAN AD EFFECT (Spearman of d on shared proteins, two-sided permutation p,",
    "\np_min is the smallest p the design can reach)\n")
for (i in which(is.finite(Sm$rho_human_ad))) with(Sm[i, ], cat(sprintf(
  "  %-40s rho %+.3f   null95 %.3f   p %.3f   p_min %.4f\n", contrast, rho_human_ad, rho_null_95,
  p_concordance, p_min_attainable)))

num <- sapply(Sm, is.numeric)
Sm[num] <- lapply(Sm[num], signif, 4)
write.csv(Sm, file.path(TABLES, "effects_summary.csv"), row.names = FALSE)

long <- do.call(rbind, lapply(C, function(c) data.frame(
  contrast = c$name, protein = names(c$d), log2fc = signif(c$beta, 5), d = signif(c$d, 5),
  p = signif(c$p, 5), resid_sd = signif(c$sd, 5), stringsAsFactors = FALSE)))
write.csv(long, gzfile(file.path(TABLES, "effects_per_protein.csv.gz")), row.names = FALSE)
cat(sprintf("\n%d contrasts written to %s\n", nrow(Sm), file.path(TABLES, "effects_summary.csv")))
