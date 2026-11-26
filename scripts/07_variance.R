# ==============================================================================
#
# Where the variance sits in each system: technical, within one genome and
# between genomes.
#
# Components, per protein
#   technical       the same material measured twice
#                     ROSMAP DLPFC, DC   GIS pairs within a TMT batch
#                     LOAD2              12 animals measured in both TMT
#                                        experiments, experiment offset removed
#   within genome   between replicates of one genome in one condition
#                     organoid   between organoids, within line x timepoint
#                     LOAD2      between animals, within genotype x diet x age
#                                band x sex
#                     Peng       between animals, within batch x genotype x age
#   between genome  between donors or lines
#                     organoid   line x timepoint means after timepoint
#                                centring, minus the within-line share
#                     iA, iN     across donor lines, one sample per line, so
#                                differentiation and technical variance are in it
#                     brain      across donors minus technical. DC uses donors
#                                outside the ROSMAP and Banner cohorts.
#
# ==============================================================================

library(data.table)

CACHE <- file.path("data", "prep_systems.rds")
DATA_DC <- file.path("data", "DiverseCohorts")
TABLES  <- "tables"
MIN_CELL <- 2      # a group needs this many samples to count as replicated

P <- readRDS(CACHE)
S <- P$samples


# ==============================================================================
# Helpers
# ==============================================================================

row_var <- function(m) {
  n <- rowSums(!is.na(m))
  v <- rowSums((m - rowMeans(m, na.rm = TRUE))^2, na.rm = TRUE) / (n - 1)
  v[n < 2] <- NA
  v
}

centre_within <- function(m, group) {
  for (g in unique(group)) {
    cols <- which(group == g)
    m[, cols] <- m[, cols, drop = FALSE] - rowMeans(m[, cols, drop = FALSE], na.rm = TRUE)
  }
  m
}

group_means <- function(m, group) {
  sapply(unique(group), function(g) rowMeans(m[, group == g, drop = FALSE], na.rm = TRUE))
}

pooled_within <- function(X, group) {
  # pooled within-group variance per protein, from groups of at least MIN_CELL
  size <- table(group)
  keep <- group %in% names(size)[size >= MIN_CELL]
  X <- X[, keep, drop = FALSE]
  group <- group[keep]
  ss <- df <- numeric(nrow(X))
  for (g in unique(group)) {
    x <- X[, group == g, drop = FALSE]
    ss <- ss + rowSums((x - rowMeans(x, na.rm = TRUE))^2, na.rm = TRUE)
    df <- df + pmax(rowSums(!is.na(x)) - 1, 0)
  }
  list(var = setNames(ifelse(df > 0, ss / df, NA), rownames(X)),
       n_groups = length(unique(group)), n_samples = sum(keep))
}

components <- list()
add <- function(system, component, v, note = "") {
  v <- v[!is.na(v)]
  components[[length(components) + 1]] <<- data.frame(
    system = system, component = component, protein = names(v), var = unname(v),
    note = note, stringsAsFactors = FALSE)
}

after_bar <- function(x) ifelse(grepl("|", x, fixed = TRUE), sub("^[^|]*\\|([^|]*).*$", "\\1", x), x)


# ==============================================================================
# Human brain
# ==============================================================================

BR <- P$BR_bc
tech_br <- P$tech_sd^2
tech_br <- tech_br[!duplicated(names(tech_br))]
add("ROSMAP DLPFC", "technical", tech_br, "within-batch GIS pairs")
add("ROSMAP DLPFC", "between genome", pmax(row_var(BR) - tech_br[rownames(BR)], 0),
    "across donors minus technical")
add("ROSMAP DLPFC", "total", row_var(BR), "across donors")

dc_technical <- function(file) {
  # technical variance from batches that carry two reference channels
  asy <- fread(file.path(DATA_DC, "AMP-AD_DiverseCohorts_assay_TMTproteomics_metadata.csv"),
               colClasses = "character", data.table = FALSE)
  control <- asy$specimenID[toupper(asy$isAssayControl) == "TRUE"]
  raw <- fread(file.path(DATA_DC, file), check.names = FALSE, data.table = FALSE,
               integer64 = "double")
  ids <- after_bar(as.character(raw[[1]]))
  X <- as.matrix(raw[, -1])
  storage.mode(X) <- "double"
  colnames(X) <- trimws(colnames(X))
  X <- X[, colnames(X) %in% control, drop = FALSE]
  X[!is.na(X) & X <= 0] <- NA
  X <- sweep(log2(X), 2, apply(log2(X), 2, median, na.rm = TRUE))
  if (anyDuplicated(ids)) {
    present <- !is.na(X)
    X <- rowsum(replace(X, !present, 0), ids, reorder = FALSE) /
      rowsum(present * 1, ids, reorder = FALSE)
  } else {
    rownames(X) <- ids
  }
  batch <- sub("\\..*$", "", colnames(X))
  pairs <- lapply(sort(unique(batch)), function(b) {
    cols <- which(batch == b)
    if (length(cols) >= 2) (X[, cols[1]] - X[, cols[2]])^2 / 2 else NULL
  })
  pairs <- do.call(cbind, pairs[!sapply(pairs, is.null)])
  list(var = rowMeans(pairs, na.rm = TRUE), n_pairs = ncol(pairs))
}

for (dc in list(list("DC frontal", P$DC_FR, P$dc_fr_meta, "normAbundances_post-FP_frontal.csv"),
                list("DC temporal", P$DC_TE, P$dc_te_meta, "normAbundances_post-FP_temporal.csv"))) {
  X <- dc[[2]][, rownames(dc[[3]])[!dc[[3]]$discovery], drop = FALSE]
  tech <- dc_technical(dc[[4]])
  cat(sprintf("%s: %d donors, technical from %d GIS pairs\n", dc[[1]], ncol(X), tech$n_pairs))
  add(dc[[1]], "technical", tech$var, sprintf("%d within-batch GIS pairs", tech$n_pairs))
  add(dc[[1]], "between genome", pmax(row_var(X) - tech$var[rownames(X)], 0),
      "across donors minus technical")
  add(dc[[1]], "total", row_var(X), "across donors")
}


# ==============================================================================
# iPSC lines and organoids
# ==============================================================================

add("iA", "between genome", row_var(P$iA_bc), "across 44 lines, includes differentiation and technical")
add("iN", "between genome", row_var(P$iN), "across 40 lines, includes differentiation and technical")

O  <- P$O1p
os <- S[S$system == "organoid", ]
os <- os[match(colnames(O), os$sample), ]
cell <- paste(os$genome, os$age, sep = "|")
w_org <- pooled_within(O, cell)
add("organoid", "within genome", w_org$var, "between organoids within line x timepoint")

between_lines <- function(strata) {
  # variance of line x timepoint means after removing the strata means, minus
  # what the within-line variance contributes to a mean of n organoids
  means <- group_means(centre_within(O, strata), cell)
  n0 <- mean(table(cell))
  var_means <- rowSums(means^2, na.rm = TRUE) / (rowSums(!is.na(means)) - length(unique(strata)))
  pmax(var_means - w_org$var / n0, 0)
}
b_org    <- between_lines(os$age)
b_org_dx <- between_lines(paste(os$condition, os$age, sep = "|"))
add("organoid", "between genome", b_org,
    "line x timepoint means after timepoint centring, minus within/n")
add("organoid", "between genome, within diagnosis", b_org_dx,
    "as above, after removing diagnosis x timepoint means")
add("organoid", "total", b_org + w_org$var, "one organoid per line, between plus within")
cat(sprintf("organoid: %d line x timepoint cells, %d organoids, about %.1f per cell\n",
            w_org$n_groups, w_org$n_samples, mean(table(cell))))


# ==============================================================================
# Mouse
# ==============================================================================

ls <- S[S$system == "LOAD2", ]
ls <- ls[match(colnames(P$M_L2), ls$sample), ]
band <- cut(as.numeric(ls$age), c(0, 8, 16, 30), labels = c("young", "mid", "old"))
w_l2 <- pooled_within(P$M_L2, paste(ls$genotype, ls$condition, band, ls$sex, sep = "|"))
add("LOAD2", "within genome", w_l2$var, "between animals within genotype x diet x age band x sex")
cat(sprintf("LOAD2: %d cells, %d animals\n", w_l2$n_groups, w_l2$n_samples))

l24 <- S[S$system == "LOAD2_24", ]
bridge <- sort(intersect(l24$animal, ls$animal))
A <- P$M_L2[, ls$sample[match(bridge, ls$animal)]]
B <- P$M_L2_24[, l24$sample[match(bridge, l24$animal)]]
common <- intersect(rownames(A), rownames(B))
D <- A[common, ] - B[common, ]
t_l2 <- row_var(D - rowMeans(D, na.rm = TRUE)) / 2
add("LOAD2", "technical", t_l2, sprintf("%d animals measured in both TMT experiments", length(bridge)))
add("LOAD2", "within genome, biological", pmax(w_l2$var[common] - t_l2, 0),
    "within genome minus technical")
cat(sprintf("LOAD2 bridging animals: %d\n", length(bridge)))

ps <- S[S$system == "Peng", ]
ps <- ps[match(colnames(P$M_PG), ps$sample), ]
w_pg <- pooled_within(P$M_PG, paste(ps$batch, ps$genotype, ps$age, sep = "|"))
add("Peng", "within genome", w_pg$var, "between animals within batch x genotype x age")
cat(sprintf("Peng: %d cells, %d animals\n", w_pg$n_groups, w_pg$n_samples))

V <- do.call(rbind, components)
V$sd <- sqrt(V$var)


# ==============================================================================
# Summary
# ==============================================================================

sd_of <- function(system, component) {
  d <- V[V$system == system & V$component == component, ]
  setNames(d$sd, d$protein)
}
q <- function(x, p) unname(quantile(x, p))

UNITS <- c("ROSMAP DLPFC" = ncol(BR), "DC frontal" = sum(!P$dc_fr_meta$discovery),
           "DC temporal" = sum(!P$dc_te_meta$discovery), "iA" = ncol(P$iA_bc), "iN" = ncol(P$iN),
           "organoid" = length(unique(os$genome)), "LOAD2" = w_l2$n_samples, "Peng" = w_pg$n_samples)
TISSUE <- c("ROSMAP DLPFC" = "human brain", "DC frontal" = "human brain",
            "DC temporal" = "human brain", "iA" = "iAstrocytes", "iN" = "iNeurons",
            "organoid" = "cerebral organoids", "LOAD2" = "mouse brain", "Peng" = "mouse brain")
ASSAY <- c("ROSMAP DLPFC" = "TMT ratio", "DC frontal" = "TMT ratio", "DC temporal" = "TMT ratio",
           "iA" = "TMT ratio", "iN" = "label-free", "organoid" = "DIA", "LOAD2" = "TMT ratio",
           "Peng" = "TMT intensity")

spec <- function(parameter, system, condition, value, lower = NA, upper = NA, range = NA,
                 feeds = "05") {
  data.frame(parameter = parameter, feeds = feeds, dataset = system,
             cell_type = unname(TISSUE[system]), assay = unname(ASSAY[system]),
             condition = condition, value = value, lower = lower, upper = upper,
             range = range,
             n_donors = if (grepl("technical", condition)) NA else unname(UNITS[system]),
             stringsAsFactors = FALSE)
}
specs <- list()

cat("\nMEDIAN SD (log2) PER COMPONENT, all proteins each system quantifies\n")
cat(sprintf("  %-14s %-34s %7s %6s\n", "system", "component", "median", "n"))
keys <- unique(V[, c("system", "component")])
for (i in seq_len(nrow(keys))) {
  d <- sd_of(keys$system[i], keys$component[i])
  cat(sprintf("  %-14s %-34s %7.3f %6d\n", keys$system[i], keys$component[i], median(d), length(d)))
  specs[[length(specs) + 1]] <- spec("sd_log2", keys$system[i], keys$component[i], median(d),
                                     q(d, .25), q(d, .75), "quartiles across proteins")
}

cat("\nPROPORTIONS WITHIN ONE SYSTEM, median over proteins\n")
t <- sd_of("LOAD2", "technical"); w <- sd_of("LOAD2", "within genome")
j <- intersect(names(t), names(w))
share <- pmin(t[j]^2 / w[j]^2, 1)
cat(sprintf("  LOAD2        technical share of within-genome variance   %.2f   (n=%d)\n",
            median(share), length(j)))
specs[[length(specs) + 1]] <- spec("technical_share", "LOAD2", "of within-genome variance",
                                   median(share), q(share, .25), q(share, .75),
                                   "quartiles across proteins")
for (tag in c("ROSMAP DLPFC", "DC frontal", "DC temporal")) {
  t <- sd_of(tag, "technical"); b <- sd_of(tag, "between genome")
  j <- intersect(names(t), names(b))
  share <- t[j]^2 / (t[j]^2 + b[j]^2)
  cat(sprintf("  %-12s technical share of total                    %.2f   (n=%d)\n",
              tag, median(share), length(j)))
  specs[[length(specs) + 1]] <- spec("technical_share", tag, "of total variance", median(share),
                                     q(share, .25), q(share, .75), "quartiles across proteins")
}
w <- sd_of("organoid", "within genome")
for (comp in c("between genome", "between genome, within diagnosis")) {
  b <- sd_of("organoid", comp)
  j <- intersect(names(w), names(b))
  icc <- b[j]^2 / (w[j]^2 + b[j]^2)
  label <- if (comp == "between genome") "all lines" else "within diagnosis"
  cat(sprintf("  organoid     line ICC, %-18s                  %.2f   [quartiles %.2f, %.2f]  (n=%d)\n",
              label, median(icc), q(icc, .25), q(icc, .75), length(j)))
  cat(sprintf("  organoid     effective sample size of 5 organoids of one line (%s): %.2f\n",
              label, 5 / (1 + 4 * median(icc))))
  specs[[length(specs) + 1]] <- spec("line_icc", "organoid", label, median(icc), q(icc, .25),
                                     q(icc, .75), "quartiles across proteins", feeds = "01 05")
}

cat("\nTMT RATIO SYSTEMS ON SHARED PROTEINS, SD relative to the DC frontal between-donor SD\n")
ref <- sd_of("DC frontal", "between genome")
pairs <- list(c("DC temporal", "between genome"), c("ROSMAP DLPFC", "between genome"),
              c("iA", "between genome"), c("LOAD2", "within genome"),
              c("LOAD2", "within genome, biological"), c("LOAD2", "technical"),
              c("DC frontal", "technical"))
shared <- Reduce(intersect, c(list(names(ref)), lapply(pairs, function(p) names(sd_of(p[1], p[2])))))
cat(sprintf("  shared proteins: %d\n", length(shared)))
for (p in pairs) {
  r <- sd_of(p[1], p[2])[shared] / ref[shared]
  r <- r[is.finite(r)]
  cat(sprintf("  %-14s %-28s %5.2f  [quartiles %.2f, %.2f]\n", p[1], p[2], median(r), q(r, .25), q(r, .75)))
  specs[[length(specs) + 1]] <- spec("sd_relative_to_brain", p[1], p[2], median(r), q(r, .25),
                                     q(r, .75), "quartiles across shared proteins")
}
cat("\n  organoid (DIA) and Peng (TMT intensities) are not in this block because their scale differs.\n")

# SD across donors or lines as a study with one sample per donor would see it
for (s in list(c("DC frontal", "total"), c("ROSMAP DLPFC", "total"), c("iA", "between genome"),
               c("iN", "between genome"), c("organoid", "total"))) {
  d <- sd_of(s[1], s[2])
  specs[[length(specs) + 1]] <- spec("line_sd_log2", s[1], "one sample per donor or line", median(d),
                                     q(d, .25), q(d, .75), "quartiles across proteins")
}

specs <- do.call(rbind, specs)
num <- sapply(specs, is.numeric)
specs[num] <- lapply(specs[num], signif, 4)
write.csv(specs, file.path(TABLES, "variance_specs.csv"), row.names = FALSE)

V$var <- signif(V$var, 5)
V$sd  <- signif(V$sd, 5)
write.csv(V, gzfile(file.path(TABLES, "variance_per_protein.csv.gz")), row.names = FALSE)
cat(sprintf("\n%d rows written to %s\n", nrow(specs), file.path(TABLES, "variance_specs.csv")))
