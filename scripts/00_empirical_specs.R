# ==============================================================================
# 00_empirical_specs.R
#
# Draws the specifications for the simulations from published iPSC datasets.
# So scripts 01 to 05 use empirical specs inputs (donor ICC, effect scale, interaction SD) rather than 
# made up numbers
# Terms
#   donor      the person/patient/organism
#   line       one iPSC clone derived from a donor.
#   replicate  a sample grown from one line (an organoid, a well)
#
#
# The public datasets have one sample per line, so the lines of a donor stand
# in as its replicates. With donor (D), clone (C) and residual (E) variance,
#   ICC from lines of a donor = D / (D + C + E)
#   ICC with one line per donor, organoids  = (D + C) / (D + C + E)
# so a line-based ICC is a lower bound for the ICC of a one-line-per-donor
# design, before any shared batch is counted.
#
# PART A  Donor ICC in iPSCs        -> scripts 01, 02, 03, 05
#         Mirauta et al. 2020        RNA-seq and proteomics of the same lines,
#                                    51 donors with two lines
#         Kilpinen et al. 2017       expression arrays, 165 donors with two or
#                                    three lines
#         Kilpinen et al. 2017       RNA-seq (Expression Atlas E-ENAD-35),
#                                    53 donors with two lines, no batch label
#         Carcamo-Orive et al. 2017  RNA-seq, 100 donors with two to seven lines
#
# PART B  SD of log2 abundance across lines, one line per donor  -> script 05
#         The six datasets above and below. Turns a fold change into SD units.
#
# PART C  Donor-by-treatment interaction in differentiated cells -> script 04
#         Alasoo et al. 2018   iPSC macrophages, 84 donors, naive / IFNg /
#                              Salmonella / both
#         Knowles et al. 2018  iPSC cardiomyocytes, 45 donors, intervention
#
# Interaction (PART C).
# Neither dataset has replicate samples of a donor within a condition, so the
# donor-specific response cannot be separated from noise in one treated-minus-
# control difference. Both have four samples per donor that form two response
# contrasts with no sample in common:
#     r1 = B - A      r2 = D - C
# Sample noise is independent between r1 and r2, so any covariance of r1 and r2
# across donors is response that belongs to the donor:
#     tau^2     = cov(r1, r2)              variance of the true donor response
#     sigma_e^2 = var(r2 - r1) / 4         residual variance of one sample
# In the notation of script 04, tau^2 = 2 * sigma_int^2, and the share of
# donors whose true response is opposite to the mean is pnorm(-|beta| / tau).
# Both estimates are conservative for the interaction. tau^2 only counts
# donor response common to the two contrasts, and sigma_e^2 also absorbs any
# donor response that differs between them.
# ==============================================================================

set.seed(20261124)

DATA_DIR <- "data"
TABLES_DIR <- "tables"

MIN_EXPR  <- 10      # RNA-seq gene kept if median CPM (or TPM) is at least this
N_PERM <- 20      # random pairings or permutations for the negative controls
FC_BREAKS <- c(1, 1.1, 1.2, 1.5, 2, Inf)   # bins of the mean fold change
FC_LABELS <- c("<1.1", "1.1-1.2", "1.2-1.5", "1.5-2", ">=2")
FOLD_CHANGES <- c(1.05, 1.1, 1.2, 1.5)     # as in 05_disease_contrast.R Part G
ARRAY_DETECTION <- 0.99    # array probe counted as detected above this value
ARRAY_MIN_SHARE <- 0.50    # and kept if detected in this share of arrays


# ==============================================================================
# Helpers
# ==============================================================================

read_tab <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path) else path
  read.delim(con, check.names = FALSE, stringsAsFactors = FALSE)
}

read_atlas <- function(path) {
  # Expression Atlas query result: comment lines, Gene ID, Gene Name, one
  # column per line. Blank cells are values below the Atlas cutoff.
  x <- read.delim(path, comment.char = "#", check.names = FALSE,
                  stringsAsFactors = FALSE)
  m <- as.matrix(x[, -(1:2)])
  rownames(m) <- x[[1]]
  colnames(m) <- sub(",.*$", "", colnames(m))   # drop ", skin" / ", blood"
  m
}

hipsci_donor <- function(line) sub("_[0-9]+$", "", sub("^HPSI[0-9]+i-", "", line))

centre_within <- function(x, block) {
  # subtract the block mean from each column of a genes x donors matrix
  b <- factor(block)
  bm <- t(rowsum(t(x), b)) / rep(as.numeric(table(b)), each = nrow(x))
  x - bm[, as.integer(b), drop = FALSE]
}

centred_contrasts <- function(A, B, C, D, block = NULL) {
  # A, B, C, D: genes x donors matrices of log2 expression, same donor order.
  # r1 = B - A and r2 = D - C share no sample.
  # block: optional batch label per donor. The contrasts are then centred on
  # the batch mean, and donors alone in a batch are dropped.
  r1 <- B - A
  r2 <- D - C
  n  <- ncol(r1)
  if (is.null(block)) block <- rep("all", n)
  tb <- table(block)
  ok <- block %in% names(tb)[tb >= 2]
  list(
    beta = (rowMeans(r1) + rowMeans(r2)) / 2,          # mean response, log2
    r1c = centre_within(r1[, ok, drop = FALSE], block[ok]),
    r2c = centre_within(r2[, ok, drop = FALSE], block[ok]),
    n = sum(ok),
    df = sum(ok) - length(unique(block[ok]))
  )
}

split_response <- function(A, B, C, D, block = NULL) {
  cc <- centred_contrasts(A, B, C, D, block)
  data.frame(
    gene = if (is.null(rownames(A))) seq_len(nrow(A)) else rownames(A),
    beta = cc$beta,
    tau2 = rowSums(cc$r1c * cc$r2c) / cc$df,           # donor response variance
    sigma_e2 = rowSums((cc$r2c - cc$r1c)^2) / cc$df / 4,   # residual variance
    n_donors = cc$n,
    df = cc$df,
    stringsAsFactors = FALSE
  )
}

add_derived <- function(e) {
  e$tau <- sqrt(pmax(e$tau2, 0))
  e$sigma_e <- sqrt(e$sigma_e2)
  e$sigma_int <- e$tau / sqrt(2)
  # in units of the residual SD, as in script 04 where sigma_res = 1
  e$beta_resid <- abs(e$beta) / e$sigma_e
  e$sigma_int_resid <- e$sigma_int / e$sigma_e
  # share of donors whose true response is opposite to the mean response
  e$frac_opposite <- ifelse(e$tau > 0, pnorm(-abs(e$beta) / e$tau), 0)
  e$fc_bin <- cut(2^abs(e$beta), FC_BREAKS, labels = FC_LABELS, right = FALSE)
  e
}

permutation_null <- function(A, B, C, D, block = NULL, n_perm = N_PERM) {
  # Negative control. The second contrast is reassigned to donors at random,
  # after centring on the batch mean. This breaks the pairing of the two
  # contrasts within a donor and keeps everything else, so the covariance
  # should centre on zero and be positive for about half of the genes.
  cc  <- centred_contrasts(A, B, C, D, block)
  out <- matrix(NA_real_, n_perm, 2,
                dimnames = list(NULL, c("median_tau2", "share_positive")))
  for (i in seq_len(n_perm)) {
    t2 <- rowSums(cc$r1c * cc$r2c[, sample(cc$n), drop = FALSE]) / cc$df
    out[i, ] <- c(median(t2), mean(t2 > 0))
  }
  out
}

summarise_by_effect <- function(e) {
  rows <- lapply(split(e, e$fc_bin), function(x) data.frame(
    fc_bin = x$fc_bin[1],
    n_genes = nrow(x),
    median_abs_log2fc = median(abs(x$beta)),
    median_sigma_e = median(x$sigma_e),
    median_tau = median(x$tau),
    pct_tau2_positive = 100 * mean(x$tau2 > 0),
    beta_resid_median = median(x$beta_resid),
    sigma_int_resid_q25 = unname(quantile(x$sigma_int_resid, 0.25)),
    sigma_int_resid_med = median(x$sigma_int_resid),
    sigma_int_resid_q75 = unname(quantile(x$sigma_int_resid, 0.75)),
    pct_opposite_median = 100 * median(x$frac_opposite),
    pct_genes_opposite_over_10 = 100 * mean(x$frac_opposite > 0.10)
  ))
  do.call(rbind, rows)
}

run_contrast <- function(dataset, contrast, batch_label, A, B, C, D, block = NULL) {
  e <- add_derived(split_response(A, B, C, D, block))
  null <- permutation_null(A, B, C, D, block)
  cat(sprintf("\n%s, %s, batch: %s\n", dataset, contrast, batch_label))
  cat(sprintf("  donors %d, genes %d\n", e$n_donors[1], nrow(e)))
  cat(sprintf("  donor response variance above zero: %.0f%% of genes (permuted donors: %.0f%%, range %.0f%%-%.0f%%)\n",
              100 * mean(e$tau2 > 0), 100 * median(null[, 2]),
              100 * min(null[, 2]), 100 * max(null[, 2])))
  cat(sprintf("  median tau2 %.4f (permuted donors: %.4f)\n",
              median(e$tau2), median(null[, 1])))
  cat(sprintf("  sigma_int in residual SD units: median %.2f, quartiles %.2f to %.2f\n",
              median(e$sigma_int_resid),
              quantile(e$sigma_int_resid, 0.25), quantile(e$sigma_int_resid, 0.75)))
  s <- summarise_by_effect(e)
  s <- cbind(dataset = dataset, contrast = contrast, batch = batch_label,
             n_donors = e$n_donors[1], s,
             pct_tau2_positive_permuted = 100 * median(null[, 2]))
  rownames(s) <- NULL
  shown <- s[, c("fc_bin", "n_genes", "beta_resid_median", "sigma_int_resid_q25",
                 "sigma_int_resid_med", "sigma_int_resid_q75", "pct_opposite_median",
                 "pct_genes_opposite_over_10")]
  shown[-(1:2)] <- round(shown[-(1:2)], 2)
  print(shown, row.names = FALSE)
  list(summary = s,
       genes = cbind(dataset = dataset, contrast = contrast, batch = batch_label, e))
}

line_sd_row <- function(dataset, cell_type, assay, unit, log2_mat) {
  # log2_mat: features x lines, one line per donor
  s <- apply(log2_mat, 1, sd)
  row <- data.frame(dataset = dataset, cell_type = cell_type, assay = assay, unit = unit,
                    n_lines = ncol(log2_mat), n_features = nrow(log2_mat),
                    sd_q25 = unname(quantile(s, 0.25)), sd_median = median(s),
                    sd_q75 = unname(quantile(s, 0.75)))
  for (fc in FOLD_CHANGES)
    row[[sprintf("fc_%s_in_sd", fc)]] <- log2(fc) / row$sd_median
  row
}

gene_icc <- function(mat, donor, batch = NULL) {
  # One-way random-effects ICC, all genes at once, truncated at 0.
  # batch = NULL  between-donor mean square around the grand mean, as
  #               compute_gene_icc() in script 03 (k0 in place of the mean
  #               group size, identical when donors have equal numbers of lines)
  # batch given   donors nested in batch, between-donor mean square around the
  #               batch mean, as compute_gene_icc_within_batch() in script 03.
  #               Every line of a donor must sit in one batch.
  d <- factor(donor)
  n <- as.numeric(table(d))
  a <- nlevels(d)
  N <- length(d)
  dm <- t(rowsum(t(mat), d)) / rep(n, each = nrow(mat))
  ms_within <- rowSums((mat - dm[, as.integer(d), drop = FALSE])^2) / (N - a)
  if (is.null(batch)) {
    centre <- matrix(rowMeans(mat), nrow(mat), a)
    nb <- 1
    k0 <- (N - sum(n^2) / N) / (a - 1)
  } else {
    stopifnot(all(tapply(as.character(batch), d, function(v) length(unique(v))) == 1))
    b  <- factor(batch)
    bd <- factor(tapply(as.character(batch), d, function(v) v[1]), levels = levels(b))
    bm <- t(rowsum(t(mat), b)) / rep(as.numeric(table(b)), each = nrow(mat))
    centre <- bm[, as.integer(bd), drop = FALSE]
    nb <- nlevels(b)
    k0 <- (N - sum(tapply(n^2, bd, sum) / tapply(n, bd, sum))) / (a - nb)
  }
  ms_between <- rowSums(sweep((dm - centre)^2, 2, n, "*")) / (a - nb)
  pmax((ms_between - ms_within) / (ms_between + (k0 - 1) * ms_within), 0)
}

icc_summary <- function(dataset, assay, grouping, icc, n_donors, n_lines,
                        n_batches = NA, same_donors_batch_ignored = NA) {
  row <- data.frame(dataset = dataset, assay = assay, grouping = grouping,
                    n_donors = n_donors, n_lines = n_lines, n_batches = n_batches,
                    n_features = length(icc),
                    icc_q25 = unname(quantile(icc, 0.25)), icc_median = median(icc),
                    icc_q75 = unname(quantile(icc, 0.75)),
                    pct_features_icc_over_20 = 100 * mean(icc > 0.20),
                    pct_features_icc_over_50 = 100 * mean(icc > 0.50),
                    icc_median_same_donors_batch_ignored = same_donors_batch_ignored,
                    control_min = NA, control_max = NA)
  cat(sprintf("  %-52s %3d donors  median ICC = %5.1f%% [quartiles %.1f%% to %.1f%%]",
              grouping, n_donors, 100 * row$icc_median, 100 * row$icc_q25,
              100 * row$icc_q75))
  if (!is.na(same_donors_batch_ignored))
    cat(sprintf("  (same donors, batch ignored: %.1f%%)", 100 * same_donors_batch_ignored))
  cat("\n")
  row
}

control_summary <- function(dataset, assay, grouping, medians, n_lines, n_features) {
  # negative control: median over pairings, with the range across pairings
  cat(sprintf("  %-52s             median ICC = %5.1f%% [%.1f%% to %.1f%% across %d pairings]\n",
              grouping, 100 * median(medians), 100 * min(medians), 100 * max(medians),
              length(medians)))
  data.frame(dataset = dataset, assay = assay, grouping = grouping,
             n_donors = NA, n_lines = n_lines, n_batches = NA, n_features = n_features,
             icc_q25 = NA, icc_median = median(medians), icc_q75 = NA,
             pct_features_icc_over_20 = NA, pct_features_icc_over_50 = NA,
             icc_median_same_donors_batch_ignored = NA,
             control_min = min(medians), control_max = max(medians))
}

pseudo_donors <- function(indiv, group, size) {
  # Groups of `size` lines from different individuals, formed within `group`.
  # Lines are shuffled within group and cut into consecutive sets. A set that
  # holds two lines of one individual is dropped. Returns NA for unused lines.
  out <- rep(NA_integer_, length(indiv))
  id <- 0L
  for (g in unique(group)) {
    idx <- which(group == g)
    if (length(idx) > 1) idx <- sample(idx)
    for (s in seq_len(length(idx) %/% size)) {
      set <- idx[(s - 1) * size + seq_len(size)]
      if (!anyDuplicated(indiv[set])) {
        id <- id + 1L
        out[set] <- id
      }
    }
  }
  out
}


donor_icc_block <- function(dataset, assay, mat, donor, batches = list(),
                            control_batch = NULL, control_size = 2) {
  # mat: features x lines (log2). donor: donor of each line.
  # batches: named list of batch labels per line, e.g. list("TMT batch" = b).
  # control_batch: name of the entry of `batches` used for the negative control.
  # The ICC rows use donors with two or more lines. The control rows use every
  # line. Returns one row per grouping.
  multi <- donor %in% names(which(table(donor) >= 2))
  m_mat <- mat[, multi, drop = FALSE]
  m_donor <- donor[multi]
  cat(sprintf("\n%s, %s: %d lines, %d donors with two or more lines, %d features\n",
              dataset, assay, length(donor), length(unique(m_donor)), nrow(mat)))

  rows <- list(icc_summary(dataset, assay, "lines grouped by donor, batch ignored",
                           gene_icc(m_mat, m_donor), length(unique(m_donor)), sum(multi)))

  for (nm in names(batches)) {
    b <- batches[[nm]][multi]
    # donors whose lines all sit in one recorded batch, in batches that hold at
    # least two such donors
    recorded <- tapply(b, m_donor, function(v) !anyNA(v))
    one <- tapply(b, m_donor, function(v) !anyNA(v) && length(unique(v)) == 1)
    nested <- m_donor %in% names(one)[one]
    per_batch <- tapply(m_donor[nested], b[nested], function(v) length(unique(v)))
    ok <- nested & b %in% names(per_batch)[per_batch >= 2]
    cat(sprintf("  %s recorded for %d of %d donors, all lines of the donor in the same one for %d\n",
                nm, sum(recorded), length(recorded), sum(one)))
    rows[[length(rows) + 1]] <- icc_summary(
      dataset, assay, paste("donor within", nm),
      gene_icc(m_mat[, ok], m_donor[ok], b[ok]),
      length(unique(m_donor[ok])), sum(ok), length(unique(b[ok])),
      same_donors_batch_ignored = median(gene_icc(m_mat[, ok], m_donor[ok])))
  }

  # negative controls: pseudo-donors made of unrelated lines
  control <- function(group, use) replicate(N_PERM, {
    pd <- pseudo_donors(donor[use], group[use], control_size)
    median(gene_icc(mat[, use, drop = FALSE][, !is.na(pd), drop = FALSE], pd[!is.na(pd)]))
  })
  if (!is.null(control_batch)) {
    b <- batches[[control_batch]]
    rows[[length(rows) + 1]] <- control_summary(
      dataset, assay, paste("control: unrelated lines in the same", control_batch),
      control(b, !is.na(b)), sum(!is.na(b)), nrow(mat))
  }
  all_lines <- rep(TRUE, length(donor))
  rows[[length(rows) + 1]] <- control_summary(
    dataset, assay, "control: unrelated lines grouped at random",
    control(rep("all", length(donor)), all_lines), length(donor), nrow(mat))
  do.call(rbind, rows)
}


# ==============================================================================
# Load: Mirauta et al. 2020, iPSC RNA-seq and proteomics
# ==============================================================================

cat("=== Mirauta et al. 2020 ===\n")
mi_meta <- read.delim(file.path(DATA_DIR, "Mirauta", "elife-57390-fig1-data1-v3.xls"),
                      row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
mi_expr <- as.matrix(read.delim(file.path(DATA_DIR, "Mirauta", "elife-57390-fig1-data2-v3.xls"),
                                row.names = 1, check.names = FALSE, stringsAsFactors = FALSE))
# Figure 1 source data 2 is RNA expression for the genes also detected as
# protein. Source data 3 is the protein matrix.
# line to donor and batch, as in 03_validation_real_data.R. batch is the only
# processing batch the source metadata records, the TMT batch of the protein
# sample. For the RNA matrix it is a proxy for lines handled together.
mi_map <- unique(data.frame(line = mi_meta[["Public name"]], donor = mi_meta[["Donor"]],
                            batch = as.character(mi_meta[["batch"]]),
                            stringsAsFactors = FALSE))
mi_map <- mi_map[match(colnames(mi_expr), mi_map$line), ]
mi_missing <- is.na(mi_map$line)
mi_map$line[mi_missing]  <- colnames(mi_expr)[mi_missing]
mi_map$donor[mi_missing] <- sub(".*-([a-z]+)_.*", "\\1", colnames(mi_expr)[mi_missing])
stopifnot(identical(mi_map$line, colnames(mi_expr)), !anyNA(mi_expr))
cat(sprintf("RNA-seq: %d lines, %d donors, %d genes\n",
            nrow(mi_map), length(unique(mi_map$donor)), nrow(mi_expr)))

# protein intensities after the authors' batch correction and quantile
# normalisation. Columns are line@TMT batch@channel. The 202 lines of the
# metadata are kept, which drops the repeats of the reference line.
mi_prot_all <- read.delim(file.path(DATA_DIR, "Mirauta", "elife-57390-fig1-data3.xls"),
                          row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(rownames(mi_meta) %in% colnames(mi_prot_all)))
mi_prot <- as.matrix(mi_prot_all[, rownames(mi_meta)])
colnames(mi_prot) <- mi_meta[["Public name"]]
mi_prot <- log2(mi_prot[rowSums(is.na(mi_prot) | mi_prot <= 0) == 0, ])
mi_prot_map <- mi_map[match(colnames(mi_prot), mi_map$line), ]
stopifnot(!anyNA(mi_prot_map$line))
cat(sprintf("Proteomics: %d lines, %d donors, %d of %d proteins quantified in every line\n",
            ncol(mi_prot), length(unique(mi_prot_map$donor)), nrow(mi_prot), nrow(mi_prot_all)))


# ==============================================================================
# Load: Kilpinen et al. 2017, iPSC expression arrays (E-MTAB-4057)
# ==============================================================================

cat("\n=== Kilpinen et al. 2017, expression arrays ===\n")
ka_dir   <- file.path(DATA_DIR, "Kilpinen1")
ka_cache <- file.path(ka_dir, "kilpinen_ipsc_arrays.rds")
if (file.exists(ka_cache)) {
  ka <- readRDS(ka_cache)
} else {
  # One GenomeStudio file per array. The table starts at the line that begins
  # with TargetID. Signal is the AVG_Signal column, the chip is the first part
  # of the array ID in its name. iPSC files only, empty files skipped.
  ka_files <- list.files(ka_dir, pattern = "^HPSI[0-9]{4}i-.*\\.txt$", full.names = TRUE)
  ka_files <- ka_files[file.size(ka_files) > 0]
  read_array <- function(f) {
    skip <- grep("^TargetID", readLines(f, n = 50))[1] - 1
    d <- read.delim(f, skip = skip, check.names = FALSE, stringsAsFactors = FALSE)
    sig <- grep("\\.AVG_Signal-", names(d))
    det <- grep("\\.Detection", names(d))
    stopifnot(length(sig) == 1, length(det) == 1)
    list(probe = d$PROBE_ID, signal = d[[sig]], detection = d[[det]],
         line  = sub("\\.AVG_Signal-.*$", "", names(d)[sig]),
         chip  = sub("_.*$", "", sub("^.*\\.AVG_Signal-", "", names(d)[sig])))
  }
  ka_raw <- lapply(ka_files, read_array)
  probes <- ka_raw[[1]]$probe
  signal <- sapply(ka_raw, function(a) a$signal[match(probes, a$probe)])
  detect <- sapply(ka_raw, function(a) a$detection[match(probes, a$probe)])
  rownames(signal) <- probes
  colnames(signal) <- sapply(ka_raw, `[[`, "line")
  keep <- rowMeans(detect > ARRAY_DETECTION, na.rm = TRUE) >= ARRAY_MIN_SHARE &
    rowSums(is.na(signal)) == 0
  ka <- list(mat = log2(signal[keep, ]), donor = hipsci_donor(colnames(signal)),
             line = colnames(signal), chip = sapply(ka_raw, `[[`, "chip"))
  saveRDS(ka, ka_cache)
}
cat(sprintf("%d arrays (one per line), %d donors, %d probes, %d chips\n",
            ncol(ka$mat), length(unique(ka$donor)), nrow(ka$mat), length(unique(ka$chip))))


# ==============================================================================
# Load: Expression Atlas tables (Kilpinen 2 (RNA-seq), Schwartzentruber)
# ==============================================================================

cat("\n=== Expression Atlas tables ===\n")
# TPM per line. Genes need a value in every line and a median TPM >= MIN_EXPR.
atlas_one_per_donor <- function(m) m[, !duplicated(hipsci_donor(colnames(m))), drop = FALSE]
atlas_log2 <- function(m) {
  ok <- rowSums(is.na(m)) == 0
  ok[ok] <- apply(m[ok, , drop = FALSE], 1, median) >= MIN_EXPR
  log2(m[ok, , drop = FALSE])
}
sn_tpm <- read_atlas(file.path(DATA_DIR, "Schwartzentruber", "E-ENAD-33-query-results.tpmss.tsv"))
ip_tpm <- read_atlas(file.path(DATA_DIR, "Kilpinen2", "E-ENAD-35-query-results.tpmss.tsv"))
ip_tpm <- ip_tpm[, grepl("^HPSI[0-9]+i-", colnames(ip_tpm))]   # iPSC lines only
sn_donor <- hipsci_donor(colnames(sn_tpm))
ip_donor <- hipsci_donor(colnames(ip_tpm))
cat(sprintf("Kilpinen RNA-seq: %d iPSC lines, %d donors, %d with two lines\n",
            ncol(ip_tpm), length(unique(ip_donor)), sum(table(ip_donor) == 2)))
cat(sprintf("Schwartzentruber sensory neurons: %d lines, %d donors, %d with two lines\n",
            ncol(sn_tpm), length(unique(sn_donor)), sum(table(sn_donor) == 2)))


# ==============================================================================
# Load: Carcamo-Orive et al. 2017, iPSCs
# ==============================================================================

cat("\n=== Carcamo-Orive et al. 2017 ===\n")
co_dir <- file.path(DATA_DIR, "Carcamo-Orive")

# sample annotation from the GEO series matrix
co_sm <- readLines(gzfile(file.path(co_dir, "GSE79636_series_matrix.txt.gz")))
geo_row <- function(key) {
  l <- co_sm[startsWith(co_sm, key)]
  lapply(strsplit(l, "\t"), function(v) gsub('"', "", v[-1]))
}
co_char <- do.call(rbind, geo_row("!Sample_characteristics_ch1"))
geo_field <- function(key) apply(co_char, 2, function(v) {
  # the row that holds a field differs between samples, so match on its name
  hit <- v[grepl(": ", v) & tolower(gsub("[. ]", " ", sub(": .*", "", v))) == key]
  if (length(hit)) sub("^[^:]*: ", "", hit[1]) else NA
})
co_meta <- data.frame(
  gsm = geo_row("!Sample_geo_accession")[[1]],
  indiv = geo_field("indiv"),
  seq_batch = geo_field("batch"),
  reprog_batch = geo_field("reprogramming batch"),
  harvest_date = geo_field("rna harvest date"),
  stringsAsFactors = FALSE
)
co_meta <- co_meta[!is.na(co_meta$indiv), ]   # the 317 lines of the main study
co_meta$harvest_date[co_meta$harvest_date == "NA"] <- NA

# counts: one featureCounts file per line inside GSE79636_RAW.tar
co_cache <- file.path(co_dir, "counts.rds")
if (file.exists(co_cache)) {
  co_count <- readRDS(co_cache)
} else {
  co_raw <- file.path(co_dir, "GSE79636_RAW")
  if (!dir.exists(co_raw)) untar(file.path(co_dir, "GSE79636_RAW.tar"), exdir = co_raw)
  co_files <- list.files(co_raw, pattern = "_counts\\.txt\\.gz$")
  names(co_files) <- sub("_.*$", "", co_files)
  stopifnot(all(co_meta$gsm %in% names(co_files)))
  co_read <- function(g) read.delim(gzfile(file.path(co_raw, co_files[[g]])),
                                    header = FALSE, stringsAsFactors = FALSE)
  co_count <- sapply(co_meta$gsm, function(g) co_read(g)[[2]])
  rownames(co_count) <- co_read(co_meta$gsm[1])[[1]]
  saveRDS(co_count, co_cache)
}
stopifnot(identical(colnames(co_count), co_meta$gsm))

co_cpm <- t(t(co_count) / colSums(co_count)) * 1e6
co_keep <- apply(co_cpm, 1, median) >= MIN_EXPR
co_log2 <- log2(t(t(co_count[co_keep, ] + 0.5) / (colSums(co_count) + 1)) * 1e6)

co_one <- function(x) tapply(x, co_meta$indiv, function(v) length(unique(v)))
cat(sprintf("%d lines, %d individuals, %d expressed genes\n",
            nrow(co_meta), length(unique(co_meta$indiv)), sum(co_keep)))
cat(sprintf("Individuals with all lines in one sequencing batch: %d of %d (%d batches)\n",
            sum(co_one(co_meta$seq_batch) == 1), length(unique(co_meta$indiv)),
            length(unique(co_meta$seq_batch))))
cat(sprintf("Individuals with all lines in one reprogramming batch: %d of %d (%d batches)\n",
            sum(co_one(co_meta$reprog_batch) == 1), length(unique(co_meta$indiv)),
            length(unique(co_meta$reprog_batch))))
co_dated <- !is.na(co_meta$harvest_date)
cat(sprintf("RNA harvest date recorded for %d lines. Individuals with all lines harvested on one date: %d of %d\n",
            sum(co_dated),
            sum(tapply(co_meta$harvest_date[co_dated], co_meta$indiv[co_dated],
                       function(v) length(unique(v))) == 1),
            length(unique(co_meta$indiv[co_dated]))))


# ==============================================================================
# Load: Alasoo et al. 2018, iPSC-derived macrophages
# ==============================================================================

cat("\n=== Alasoo et al. 2018 ===\n")
al_meta <- read_tab(file.path(DATA_DIR, "Alasoo", "RNA_sample_metadata.txt.gz"))
al_cqn <- as.matrix(read_tab(file.path(DATA_DIR, "Alasoo", "RNA_cqn_matrix.txt.gz")))
al_count <- as.matrix(read_tab(file.path(DATA_DIR, "Alasoo", "RNA_count_matrix.txt.gz")))
stopifnot(all(al_meta$sample_id %in% colnames(al_cqn)),
          identical(rownames(al_cqn), rownames(al_count)))

al_conditions <- c("naive", "IFNg", "SL1344", "IFNg_SL1344")
stopifnot(all(table(al_meta$donor, al_meta$condition_name) == 1))

# expressed genes: median CPM >= MIN_EXPR in at least one condition
al_cpm <- t(t(al_count) / colSums(al_count)) * 1e6
al_median_cpm <- sapply(al_conditions, function(cn)
  apply(al_cpm[, al_meta$sample_id[al_meta$condition_name == cn]], 1, median))
al_keep <- apply(al_median_cpm, 1, max) >= MIN_EXPR

al_donors <- sort(unique(al_meta$donor))
al_get <- function(cn) {
  s <- al_meta[al_meta$condition_name == cn, ]
  m <- al_cqn[al_keep, s$sample_id[match(al_donors, s$donor)]]
  colnames(m) <- al_donors
  m
}
al_naive <- al_get("naive")
al_ifng <- al_get("IFNg")
al_sl <- al_get("SL1344")
al_both <- al_get("IFNg_SL1344")
# batch: the date of stimulation, shared by the four samples of a donor
al_date <- al_meta$salmonella_date[match(al_donors, al_meta$donor)]

cat(sprintf("%d donors, %d conditions, %d expressed genes, %d stimulation dates\n",
            length(al_donors), length(al_conditions), sum(al_keep),
            length(unique(al_date))))


# ==============================================================================
# Load: Knowles et al. 2018, iPSC-derived cardiomyocytes
# ==============================================================================

cat("\n=== Knowles et al. 2018 ===\n")
kn_meta <- read_tab(file.path(DATA_DIR, "Knowles", "sample_table.txt"))
kn_expr <- as.matrix(read_tab(file.path(DATA_DIR, "Knowles", "log2_counts_per_million.txt.gz")))
stopifnot(identical(colnames(kn_expr), kn_meta$sample))

# one line per individual (one individual has two lines, the first is kept)
kn_first <- tapply(kn_meta$cell_line, kn_meta$individual, function(v) v[1])
kn_meta <- kn_meta[kn_meta$cell_line %in% kn_first, ]

kn_doses <- sort(unique(kn_meta$conc))
kn_median <- sapply(kn_doses, function(cc)
  apply(kn_expr[, kn_meta$sample[kn_meta$conc == cc]], 1, median))
kn_keep <- apply(kn_median, 1, max) >= log2(MIN_EXPR)

# lines with the four doses used for the two contrasts
kn_need <- c(0, 0.625, 2.5, 5)
kn_lines_with <- function(doses) {
  has <- tapply(kn_meta$conc, kn_meta$cell_line, function(v) all(doses %in% v))
  names(has)[has]
}
kn_lines <- kn_lines_with(kn_need)
kn_get <- function(cc, lines = kn_lines) {
  m <- kn_expr[kn_keep, kn_meta$sample[match(paste(lines, cc),
                                             paste(kn_meta$cell_line, kn_meta$conc))]]
  colnames(m) <- lines
  m
}
cat(sprintf("%d individuals, %d with doses %s, %d expressed genes\n",
            length(unique(kn_meta$individual)), length(kn_lines),
            paste(kn_need, collapse = ", "), sum(kn_keep)))


# ==============================================================================
# PART A: donor ICC in iPSCs
# ==============================================================================

cat("\n\n################ PART A: donor ICC in iPSCs ################\n")

ip_pairs <- atlas_log2(ip_tpm[, ip_donor %in% names(which(table(ip_donor) == 2))])

donor_icc <- rbind(
  donor_icc_block("Mirauta 2020", "RNA-seq", mi_expr, mi_map$donor,
                  batches = list("proteomics batch" = mi_map$batch),
                  control_batch = "proteomics batch", control_size = 2),
  donor_icc_block("Kilpinen 2017", "expression array", ka$mat, ka$donor,
                  batches = list("array chip" = ka$chip),
                  control_batch = "array chip", control_size = 2),
  donor_icc_block("Kilpinen 2017", "RNA-seq", ip_pairs, hipsci_donor(colnames(ip_pairs)),
                  control_size = 2),
  donor_icc_block("Carcamo-Orive 2017", "RNA-seq", co_log2, co_meta$indiv,
                  batches = list("sequencing batch" = co_meta$seq_batch,
                                 "reprogramming batch" = co_meta$reprog_batch,
                                 "RNA harvest date" = co_meta$harvest_date),
                  control_batch = "sequencing batch", control_size = 3)
)
# The random number state is restored after the protein block, so the random
# pairings and permutations of the other datasets do not depend on it.
rng_state <- .Random.seed
donor_icc <- rbind(
  donor_icc,
  donor_icc_block("Mirauta 2020", "proteomics", mi_prot, mi_prot_map$donor,
                  batches = list("TMT batch" = mi_prot_map$batch),
                  control_batch = "TMT batch", control_size = 2)
)
.Random.seed <- rng_state
rownames(donor_icc) <- NULL
write.csv(donor_icc, file.path(TABLES_DIR, "ext_donor_icc.csv"), row.names = FALSE)


# ==============================================================================
# PART B: SD of log2 abundance across lines
# ==============================================================================

cat("\n\n################ PART B: SD across lines, one line per donor ################\n")

sn_log2 <- atlas_log2(atlas_one_per_donor(sn_tpm))
ip_log2 <- atlas_log2(atlas_one_per_donor(ip_tpm))
shared <- intersect(rownames(sn_log2), rownames(ip_log2))

line_sd <- rbind(
  # The log base of the Mirauta RNA matrix is not stated. Against Kilpinen
  # log2 TPM on shared genes and lines its slope is 0.88, so it is taken as log2.
  line_sd_row("Mirauta 2020", "iPSCs", "RNA-seq", "log2 expression",
              mi_expr[, !duplicated(mi_map$donor)]),
  line_sd_row("Mirauta 2020", "iPSCs", "proteomics", "log2 intensity",
              mi_prot[, !duplicated(mi_prot_map$donor)]),
  line_sd_row("Kilpinen 2017", "iPSCs", "expression array", "log2 intensity",
              ka$mat[, !duplicated(ka$donor)]),
  line_sd_row("Kilpinen 2017", "iPSCs", "RNA-seq", "log2 TPM", ip_log2),
  line_sd_row("Carcamo-Orive 2017", "iPSCs", "RNA-seq", "log2 CPM",
              co_log2[, !duplicated(co_meta$indiv)]),
  line_sd_row("Alasoo 2018", "macrophages, naive", "RNA-seq", "log2 cqn", al_naive),
  line_sd_row("Knowles 2018", "cardiomyocytes, untreated", "RNA-seq", "log2 CPM",
              kn_expr[kn_keep, kn_meta$sample[kn_meta$conc == 0]]),
  line_sd_row("Schwartzentruber 2018", "sensory neurons", "RNA-seq", "log2 TPM", sn_log2),
  line_sd_row("Schwartzentruber 2018", "sensory neurons, genes shared with iPSCs",
              "RNA-seq", "log2 TPM", sn_log2[shared, ]),
  line_sd_row("Kilpinen 2017", "iPSCs, genes shared with sensory neurons",
              "RNA-seq", "log2 TPM", ip_log2[shared, ])
)
rownames(line_sd) <- NULL
print(line_sd[, 1:9], digits = 3, row.names = FALSE)
write.csv(line_sd, file.path(TABLES_DIR, "ext_line_sd.csv"), row.names = FALSE)


# ==============================================================================
# PART C: donor-by-treatment interaction
# ==============================================================================

cat("\n\n################ PART C: donor-by-treatment interaction ################\n")

partA <- list(
  # IFNg response: without Salmonella (r1) and with Salmonella (r2)
  run_contrast("Alasoo 2018, macrophages", "IFNg", "within stimulation date",
               A = al_naive, B = al_ifng, C = al_sl, D = al_both, block = al_date),
  run_contrast("Alasoo 2018, macrophages", "IFNg", "ignored",
               A = al_naive, B = al_ifng, C = al_sl, D = al_both),
  # Salmonella response: without IFNg (r1) and with IFNg (r2)
  run_contrast("Alasoo 2018, macrophages", "Salmonella", "within stimulation date",
               A = al_naive, B = al_sl, C = al_ifng, D = al_both, block = al_date),
  run_contrast("Alasoo 2018, macrophages", "Salmonella", "ignored",
               A = al_naive, B = al_sl, C = al_ifng, D = al_both),
  # Doxorubicin response: 2.5 against 0 uM (r1) and 5 against 0.625 uM (r2).
  # There is one untreated sample per line, so the second contrast cannot use
  # it. No batch label is available, so any batch difference in response is
  # counted as donor response here.
  run_contrast("Knowles 2018, cardiomyocytes", "Doxorubicin", "not available",
               A = kn_get(0), B = kn_get(2.5), C = kn_get(0.625), D = kn_get(5))
)

# Sensitivity of the Knowles estimate to the choice of doses. The four doses
# can be paired into two contrasts in several ways, and the result is reported
# for each so that the dependence on this choice is visible.
kn_pairings <- list(c(0, 2.5, 0.625, 5), c(0, 5, 0.625, 2.5), c(0, 1.25, 0.625, 2.5),
                    c(0, 0.625, 1.25, 2.5), c(0, 1.25, 2.5, 5))
knowles_dose_pairing <- do.call(rbind, lapply(kn_pairings, function(d) {
  L <- kn_lines_with(d)
  A <- kn_get(d[1], L); B <- kn_get(d[2], L); C <- kn_get(d[3], L); D <- kn_get(d[4], L)
  e <- add_derived(split_response(A, B, C, D))
  null <- permutation_null(A, B, C, D)
  data.frame(
    r1 = sprintf("%s - %s", d[2], d[1]), r2 = sprintf("%s - %s", d[4], d[3]),
    n_lines = length(L),
    median_tau2 = median(e$tau2),
    pct_tau2_positive = 100 * mean(e$tau2 > 0),
    pct_tau2_positive_permuted_min = 100 * min(null[, 2]),
    pct_tau2_positive_permuted_max = 100 * max(null[, 2]),
    sigma_int_resid_med = median(e$sigma_int_resid))
}))
cat("\nKnowles, sensitivity to the pairing of doses\n")
print(cbind(knowles_dose_pairing[1:3], round(knowles_dose_pairing[-(1:3)], 3)), row.names = FALSE)
write.csv(knowles_dose_pairing,
          file.path(TABLES_DIR, "ext_knowles_dose_pairing.csv"), row.names = FALSE)

interaction_by_effect <- do.call(rbind, lapply(partA, `[[`, "summary"))
interaction_per_gene  <- do.call(rbind, lapply(partA, `[[`, "genes"))

write.csv(interaction_by_effect,
          file.path(TABLES_DIR, "ext_interaction_by_effect.csv"), row.names = FALSE)
num <- sapply(interaction_per_gene, is.numeric)
interaction_per_gene[num] <- lapply(interaction_per_gene[num], signif, 5)
write.csv(interaction_per_gene,
          gzfile(file.path(TABLES_DIR, "ext_interaction_per_gene.csv.gz")), row.names = FALSE)


# ==============================================================================
# Specifications table for scripts 01 to 05
# ==============================================================================

spec <- function(parameter, feeds, dataset, cell_type, assay, condition,
                 value, lower = NA, upper = NA, range = NA, n_donors = NA) {
  data.frame(parameter = parameter, feeds = feeds, dataset = dataset,
             cell_type = cell_type, assay = assay, condition = condition,
             value = value, lower = lower, upper = upper, range = range,
             n_donors = n_donors, stringsAsFactors = FALSE)
}

is_control <- grepl("^control", donor_icc$grouping)
ie <- interaction_by_effect
ie_cell <- sub("^.*, ", "", ie$dataset)
ie_set  <- sub(",.*$", "", ie$dataset)
ie_cond <- sprintf("%s, batch %s, mean fold change %s", ie$contrast, ie$batch, ie$fc_bin)

empirical_specs <- rbind(
  with(donor_icc[!is_control, ],
       spec("donor_icc", "01 02 03 05", dataset, "iPSCs", assay, grouping,
            icc_median, icc_q25, icc_q75, "quartiles across features", n_donors)),
  with(donor_icc[is_control, ],
       spec("control_icc", "01", dataset, "iPSCs", assay, grouping,
            icc_median, control_min, control_max,
            "lowest and highest median across pairings")),
  with(line_sd,
       spec("line_sd_log2", "05", dataset, cell_type, assay, unit,
            sd_median, sd_q25, sd_q75, "quartiles across features", n_lines)),
  spec("beta_resid", "04", ie_set, ie_cell, "RNA-seq", ie_cond,
       ie$beta_resid_median, n_donors = ie$n_donors),
  spec("sigma_int_resid", "04", ie_set, ie_cell, "RNA-seq", ie_cond,
       ie$sigma_int_resid_med, ie$sigma_int_resid_q25, ie$sigma_int_resid_q75,
       "quartiles across genes", ie$n_donors),
  spec("pct_donors_opposite", "04", ie_set, ie_cell, "RNA-seq", ie_cond,
       ie$pct_opposite_median, n_donors = ie$n_donors)
)
num <- sapply(empirical_specs, is.numeric)
empirical_specs[num] <- lapply(empirical_specs[num], signif, 4)
write.csv(empirical_specs, file.path(TABLES_DIR, "empirical_specs.csv"), row.names = FALSE)

cat(sprintf("\n%d specifications written to %s\n", nrow(empirical_specs),
            file.path(TABLES_DIR, "empirical_specs.csv")))
print(table(empirical_specs$parameter))

