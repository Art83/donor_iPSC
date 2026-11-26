# ==============================================================================
# 06_prep.R
#
# Builds the protein matrices for the human, iPSC, organoid and mouse systems
# and saves them in one file, data/prep_systems.rds. Follwing scripts will read it.
#
#
#   iA          ROSMAP iAstrocytes, TMT log2 ratio, 44 donor lines in 3 batches
#   iA_bc       the same with the batch mean of each protein removed
#   iN          ROSMAP iNeurons, label-free, log2 and median-centred
#   BR_bc       ROSMAP DLPFC, batch-corrected log2 ratio (C2), 400 donors
#   tech_sd     technical SD per protein from the GIS pairs of the raw DLPFC file
#   DC_FR, DC_TE  AMP-AD Diverse Cohorts frontal and temporal cortex
#   O1p         organoids, 5 lines x 2 timepoints x 5 organoids, DIA proteomics
#   M_L2        MODEL-AD LOAD2 mice, re-keyed on human UniProt
#   M_L2_24     the 24-month LOAD2 experiment, kept for its bridging animals
#   M_PG        Peng mice (5xFAD and APP knock-in), re-keyed on human UniProt
#   samples     one row per sample of every system
# ==============================================================================

library(data.table)

DATA_R <- file.path("data", "rosmap")
DATA_O <- file.path("data", "organoids")
DATA_M <- file.path("data", "mouse")
DATA_DC <- file.path("data", "DiverseCohorts")
HCOP <- file.path("data", "human_mouse_hcop_fifteen_column.txt")
CACHE <- file.path("data", "prep_systems.rds")

MIN_FRAC_BRAIN <- 0.50
MIN_FRAC_IN <- 0.80
MIN_FRAC_ORG <- 0.50
MIN_HCOP_SUPPORT <- 3
ROSMAP_BANNER <- c("ROS", "MAP", "Banner")


# ==============================================================================
# Helpers
# ==============================================================================

read_table <- function(path, ...) {
  as.data.frame(fread(path, check.names = FALSE, data.table = FALSE, integer64 = "double", ...))
}

as_matrix <- function(df, ids) {
  m <- as.matrix(df)
  storage.mode(m) <- "double"
  rownames(m) <- ids
  m
}

centre_columns <- function(m) sweep(m, 2, apply(m, 2, median, na.rm = TRUE))

keep_present <- function(m, min_frac) m[rowMeans(!is.na(m)) >= min_frac, , drop = FALSE]

mean_by_id <- function(m, ids) {
  if (!anyDuplicated(ids)) {
    rownames(m) <- ids
    return(m)
  }
  present <- !is.na(m)
  total <- rowsum(replace(m, !present, 0), ids, reorder = FALSE)
  count <- rowsum(present * 1, ids, reorder = FALSE)
  out <- total / count
  out[count == 0] <- NA
  out
}

centre_within <- function(m, group) {
  for (g in unique(group)) {
    cols <- which(group == g)
    m[, cols] <- m[, cols, drop = FALSE] - rowMeans(m[, cols, drop = FALSE], na.rm = TRUE)
  }
  m
}

after_bar  <- function(x) ifelse(grepl("|", x, fixed = TRUE), sub("^[^|]*\\|([^|]*).*$", "\\1", x), x)
before_bar <- function(x) sub("\\|.*$", "", x)


# ==============================================================================
# ROSMAP
# ==============================================================================

ia_meta <- read_table(file.path(DATA_R, "assay_meta_ROSMAP_iA_TMT-MS.csv"))
ia_meta$batch <- gsub("[^A-Za-z0-9_]", "", ia_meta$batch, useBytes = TRUE)   # trailing non-breaking space
fq <- read_table(file.path(DATA_R, "fastq.map.csv"))
br2sub <- unique(rbind(ia_meta[, c("BRID", "SubjectID")], fq[, c("BRID", "SubjectID")]))
br2sub <- setNames(br2sub$SubjectID, br2sub$BRID)

pheno <- ia_meta[, c("BRID", "SubjectID", "sex", "batch", "dx", "apoe_genotype",
                     "cogn_global_random_slope", "gpath", "amyloid", "tangles",
                     "nft", "pmAD", "pmi")]
pheno$apoe4 <- as.integer(grepl("4", as.character(pheno$apoe_genotype)))
rownames(pheno) <- pheno$BRID

ia <- read_table(file.path(DATA_R, "ROSMAP-iA-TMT-MS.csv"), encoding = "UTF-8")
ia_sym <- setNames(ia$symbol, ia$uniprot)
iA <- as_matrix(ia[, setdiff(names(ia), c("uniprot", "symbol"))], ia$uniprot)
colnames(iA) <- sapply(strsplit(colnames(iA), "_"), `[`, 2)
iA <- iA[, pheno$BRID]
iA_bc <- centre_within(iA, pheno$batch) + rowMeans(iA, na.rm = TRUE)

inn <- read_table(file.path(DATA_R, "data_proteomics.csv"), encoding = "UTF-8")
iN <- as_matrix(inn[, -1], inn[[1]])
colnames(iN) <- sub("d21iNprot_", "", colnames(iN), fixed = TRUE)
iN[iN == 0] <- NA
n_in <- nrow(iN)
iN <- centre_columns(log2(keep_present(iN, MIN_FRAC_IN)))
cat(sprintf("iN: %d of %d proteins kept (>=%.0f%% quantified), residual missing %.3f\n",
            nrow(iN), n_in, 100 * MIN_FRAC_IN, mean(is.na(iN))))

tq <- read_table(file.path(DATA_R, "ROSMAP_assay_proteomics_TMTquantitation_metadata.csv"))
tq <- tq[tq$TMTdataSubmissionBatch == 1, ]
tq$individualID <- ifelse(grepl("\\.R[0-9]+$", tq$specimenID),
                          sub("^.*\\.(R[0-9]+)$", "\\1", tq$specimenID), NA)
rownames(tq) <- tq$batchChannel

c2 <- read_table(file.path(DATA_R, paste0("C2.median_polish_corrected_log2_abundanceRatioCenteredOnMedian",
                                          "OfBatchMediansPerProtein_-8817x400.csv")))
BR_sym <- setNames(before_bar(c2[[1]]), after_bar(c2[[1]]))
BR <- as_matrix(c2[, -1], after_bar(c2[[1]]))
colnames(BR) <- tq[colnames(BR), "individualID"]
BR <- BR[, !duplicated(colnames(BR))]
BR_bc <- keep_present(BR, MIN_FRAC_BRAIN)

a1 <- read_table(file.path(DATA_R, paste0("A1.PD-RAW_normalized_abundance-11672x500TMTchannels-",
                                          "no_batch_correction_1_.csv")))
a1_log <- log2(as_matrix(a1[, -1], after_bar(a1[[1]])))
a1_tq <- tq[colnames(a1_log), ]
gis <- sapply(sort(unique(a1_tq$batch)), function(b) {
  ch <- a1_tq$batchChannel[a1_tq$batch == b & as.logical(a1_tq$isAssayControl)]
  (a1_log[, ch[1]] - a1_log[, ch[2]]) / sqrt(2)
})
tech_sd <- apply(gis, 1, sd, na.rm = TRUE)
tech_sd[!is.finite(tech_sd)] <- NA
cat(sprintf("brain: %d proteins x %d donors, median technical SD %.3f log2\n",
            nrow(BR_bc), ncol(BR_bc), median(tech_sd, na.rm = TRUE)))

ia_brain <- pheno$BRID[pheno$SubjectID %in% colnames(BR_bc)]
in_brain <- colnames(iN)[br2sub[colnames(iN)] %in% colnames(BR_bc)]
in_ia <- intersect(colnames(iN), colnames(iA))
cat(sprintf("donor overlaps: iA-brain %d, iN-brain %d, iN-iA %d\n",
            length(ia_brain), length(in_brain), length(in_ia)))


# ==============================================================================
# Organoids
# ==============================================================================

load_dia <- function(path, keep_cols = NULL, min_frac = MIN_FRAC_ORG) {
  df <- read_table(path, sep = "\t")
  val <- grep("\\.raw$", names(df), value = TRUE)
  up  <- sub("-.*$", "", sub(";.*$", "", as.character(df[["Protein.Group"]])))
  sym <- trimws(sub(";.*$", "", as.character(df[["Genes"]])))
  X <- as_matrix(df[, val], NULL)
  colnames(X) <- sub("\\.raw$", "", sub("^.*\\\\", "", val))
  if (!is.null(keep_cols)) X <- X[, intersect(keep_cols, colnames(X)), drop = FALSE]
  X[X == 0] <- NA
  ok <- rowMeans(!is.na(X)) >= min_frac & !is.na(up) & up != ""
  X <- centre_columns(log2(X[ok, , drop = FALSE]))
  smap <- setNames(sym, up)
  list(mat = mean_by_id(X, up[ok]), sym = smap[!duplicated(names(smap))])
}

wells <- read_table(file.path(DATA_O, "meta_wells.csv"))
org_meta <- read_table(file.path(DATA_O, "meta.csv"))
org_meta <- org_meta[!duplicated(org_meta$donor), ]
wp <- wells[wells$layer == "prot", ]
o1 <- load_dia(file.path(DATA_O, "proteomics", "report.pg_matrix.tsv"), keep_cols = wp$ms_file)
O1p <- o1$mat
O1p_sym <- o1$sym
s1p <- wp[match(colnames(O1p), wp$ms_file), ]
s1p <- cbind(s1p, org_meta[match(s1p$donor, org_meta$donor), c("condition", "sex", "genotype")])
rownames(s1p) <- s1p$ms_file
cat(sprintf("organoids: %d proteins x %d organoids, %d lines\n",
            nrow(O1p), ncol(O1p), length(unique(s1p$donor))))


# ==============================================================================
# Mouse
# ==============================================================================

# one-to-one orthologs, mouse symbol -> human symbol
hc <- read_table(HCOP, sep = "\t", colClasses = "character", quote = "")
hc <- hc[hc$human_symbol != "-" & hc$mouse_symbol != "-", ]
hc <- hc[lengths(strsplit(hc$support, ",")) >= MIN_HCOP_SUPPORT, ]
hc <- hc[!duplicated(hc[, c("human_symbol", "mouse_symbol")]), ]
one <- table(hc$human_symbol)[hc$human_symbol] == 1 & table(hc$mouse_symbol)[hc$mouse_symbol] == 1
o2o <- setNames(hc$human_symbol[one], hc$mouse_symbol[one])
cat(sprintf("\nHCOP: %d pairs with support >= %d, %d one-to-one\n",
            nrow(hc), MIN_HCOP_SUPPORT, length(o2o)))

# human symbol -> the UniProt ID the human matrices use
sm <- c(BR_sym, ia_sym, O1p_sym)
sm <- sm[!duplicated(names(sm)) & !is.na(sm)]
sm <- data.frame(up = names(sm), sym = unname(sm), stringsAsFactors = FALSE)
sm$iso <- grepl("-", sm$up, fixed = TRUE)
sm$not_brain <- !(sm$up %in% rownames(BR_bc))
sm <- sm[order(sm$sym, sm$not_brain, sm$iso, method = "radix"), ]
sm <- sm[!duplicated(sm$sym), ]
sym2up <- setNames(sm$up, sm$sym)

to_human <- function(X, msym, acc, tag) {
  # one entry per mouse gene (fewest missing, then Swiss-Prot accession), mapped
  # through a one-to-one ortholog onto the human UniProt ID
  d <- data.frame(row = seq_len(nrow(X)), msym = msym, acc = acc,
                  nmiss = rowSums(is.na(X)),
                  noncanon = !grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$", acc),
                  stringsAsFactors = FALSE)
  d <- d[!is.na(d$msym) & !(d$msym %in% c("nan", "NA")), ]
  d <- d[order(d$msym, d$nmiss, d$noncanon, method = "radix"), ]
  d <- d[!duplicated(d$msym), ]
  d$hsym <- unname(o2o[d$msym])
  d$up <- unname(sym2up[d$hsym])
  cat(sprintf("  %-8s %5d mouse genes, %5d with a one-to-one ortholog, %5d on a human UniProt ID\n",
              tag, nrow(d), sum(!is.na(d$hsym)), sum(!is.na(d$up))))
  d <- d[!is.na(d$up), ]
  d <- d[!duplicated(d$up), ]
  out <- X[d$row, , drop = FALSE]
  rownames(out) <- d$up
  out
}

cat("mouse matrices re-keyed on human UniProt\n")
L2DIR <- file.path(DATA_M, "load2")
ind <- read_table(file.path(L2DIR, "individual_animal_metadata_Jax.IU.Pitt_LOAD2.csv"),
                  colClasses = list(character = "individualID"))
ind$individualID <- trimws(ind$individualID)
ind$geno <- ifelse(grepl("hAPP-3/3_homozygous", ind$genotype), "LOAD2",
            ifelse(grepl("hAPP-3/3_WT", ind$genotype), "LOAD1",
            ifelse(grepl("hAPP-3/3_heterozygous", ind$genotype), "het_hAPP", "other")))
ind$diet <- ifelse(grepl("high", tolower(ind$rodentDiet)), "HFD", "CD")
rownames(ind) <- ind$individualID

load_load2 <- function(file, tag) {
  # rows are samples in the deposited files
  w <- read_table(file.path(L2DIR, file), colClasses = list(character = c("specimenID", "individualID")))
  X <- t(as_matrix(w[, -(1:2)], NULL))
  colnames(X) <- w$specimenID
  list(mat = to_human(X, before_bar(rownames(X)), after_bar(rownames(X)), tag),
       samples = data.frame(specimenID = w$specimenID, individualID = trimws(w$individualID),
                            stringsAsFactors = FALSE))
}
l2    <- load_load2("TMT_normAbundances_Batch_regressed_Jax.IU.Pitt_LOAD2.csv", "LOAD2")
l2_24 <- load_load2("TMT_normAbundances.Batch_regressed_24Months_Jax.IU.Pitt_LOAD2.csv", "LOAD2_24")
M_L2 <- l2$mat
M_L2_24 <- l2_24$mat

PGDIR <- file.path(DATA_M, "peng")
pg <- read_table(file.path(PGDIR, "peng_whole_proteome.csv"), encoding = "UTF-8")
parts <- strsplit(as.character(pg[[1]]), "|", fixed = TRUE)
part <- function(i) vapply(parts, function(p) {
  j <- if (i == "last") length(p) else i
  if (length(p) >= max(j, 1)) p[j] else NA_character_
}, character(1))
pg <- as_matrix(pg[, -1], NULL)
keep_pg <- !(part(1) %in% "cu")                   # custom human and mutant entries
pg[!is.na(pg) & pg <= 0] <- NA
pg <- centre_columns(log2(pg[keep_pg, ]))
M_PG <- to_human(pg, part("last")[keep_pg], part(2)[keep_pg], "Peng")
ps <- read_table(file.path(PGDIR, "peng_samples.csv"))
ps$arm <- ifelse(as.integer(substring(ps$batch, 2)) >= 6, "APPKI", "5xFAD")   # B6 holds WT only


# ==============================================================================
# AMP-AD Diverse Cohorts
# ==============================================================================

dca <- read_table(file.path(DATA_DC, "AMP-AD_DiverseCohorts_assay_TMTproteomics_metadata.csv"),
                  colClasses = "character", na.strings = c("", "NA"))
dcb <- read_table(file.path(DATA_DC, "AMP-AD_DiverseCohorts_biospecimen_metadata.csv"),
                  colClasses = "character", na.strings = c("", "NA"))
dci <- read_table(file.path(DATA_DC, "AMP-AD_DiverseCohorts_individual_metadata.csv"),
                  colClasses = "character", na.strings = c("", "NA"))
dcb <- dcb[!duplicated(dcb$specimenID), c("specimenID", "individualID", "tissue")]
dca <- cbind(dca, dcb[match(dca$specimenID, dcb$specimenID), c("individualID", "tissue")])
dca <- cbind(dca, dci[match(dca$individualID, dci$individualID), setdiff(names(dci), "individualID")])
dca <- dca[!duplicated(dca$specimenID), ]
rownames(dca) <- dca$specimenID
dca$site <- sub("_?b[0-9]+$", "", dca$batch)
dca$site[dca$site == ""] <- "stg"
# cohorts that also supply the ROSMAP DLPFC matrix and the iPSC donors
dca$discovery <- !is.na(dca$cohort) & dca$cohort %in% ROSMAP_BANNER

load_dc <- function(file, tag) {
  raw <- read_table(file.path(DATA_DC, file))
  ids <- as.character(raw[[1]])
  X <- as_matrix(raw[, -1], NULL)
  colnames(X) <- trimws(colnames(X))
  control <- toupper(dca[colnames(X), "isAssayControl"]) %in% "TRUE"
  X <- X[, colnames(X) %in% dca$specimenID & !control, drop = FALSE]
  X[!is.na(X) & X <= 0] <- NA
  X <- centre_columns(log2(X))
  ok <- rowMeans(!is.na(X)) >= MIN_FRAC_BRAIN
  X <- centre_within(X[ok, ], sub("\\..*$", "", colnames(X)))
  X <- mean_by_id(X, after_bar(ids[ok]))
  m <- dca[colnames(X), ]
  cat(sprintf("  %-9s %5d proteins x %4d samples, missing %.3f, ROSMAP and Banner cohorts %d, independent %d\n",
              tag, nrow(X), ncol(X), mean(is.na(X)), sum(m$discovery), sum(!m$discovery)))
  list(mat = X, meta = m)
}

cat("\nAMP-AD Diverse Cohorts\n")
dc_fr <- load_dc("normAbundances_post-FP_frontal.csv", "frontal")
dc_te <- load_dc("normAbundances_post-FP_temporal.csv", "temporal")


# ==============================================================================
# One table of every sample
# ==============================================================================

sample_cols <- c("system", "sample", "genome", "unit", "genotype", "condition", "sex",
                 "batch", "age", "animal", "cohort", "site", "discovery")
sample_rows <- function(...) {
  x <- data.frame(..., stringsAsFactors = FALSE)
  x[setdiff(sample_cols, names(x))] <- NA
  x[, sample_cols]
}
load2_rows <- function(tag, s) {
  a <- ind[s$individualID, ]
  sample_rows(system = tag, sample = s$specimenID, genome = "B6J.APOE4.Trem2R47H",
              unit = "animal", genotype = a$geno, condition = a$diet, sex = a$sex,
              batch = sub("\\..*$", "", s$specimenID),
              age = as.character(suppressWarnings(as.numeric(a$ageDeath))),
              animal = s$individualID)
}
dc_rows <- function(tag, m) {
  sample_rows(system = tag, sample = rownames(m), genome = m$individualID, unit = "donor",
              genotype = m$apoeGenotype, condition = m$ADoutcome, sex = m$sex,
              batch = sub("\\..*$", "", rownames(m)), age = m$ageDeath,
              cohort = m$cohort, site = m$site, discovery = m$discovery)
}

samples <- rbind(
  sample_rows(system = "iA", sample = pheno$BRID, genome = pheno$SubjectID, unit = "donor line",
              genotype = as.character(pheno$apoe_genotype), condition = pheno$dx,
              sex = pheno$sex, batch = pheno$batch),
  sample_rows(system = "iN", sample = colnames(iN), unit = "donor line",
              genome = ifelse(is.na(br2sub[colnames(iN)]), colnames(iN), br2sub[colnames(iN)])),
  sample_rows(system = "DLPFC", sample = colnames(BR_bc), genome = colnames(BR_bc), unit = "donor"),
  sample_rows(system = "organoid", sample = rownames(s1p), genome = s1p$donor, unit = "organoid",
              genotype = s1p$genotype, condition = s1p$condition, sex = s1p$sex,
              age = s1p$timepoint),
  load2_rows("LOAD2", l2$samples),
  load2_rows("LOAD2_24", l2_24$samples),
  sample_rows(system = "Peng", sample = ps$sample,
              genome = ifelse(ps$arm == "5xFAD", "C57BL/6J", "C57BL/6J (APP-KI)"),
              unit = "animal", genotype = ps$geno, condition = ps$arm,
              sex = ifelse(ps$sex == "M", "male", ifelse(ps$sex == "F", "female", ps$sex)),
              batch = ps$batch, age = as.character(ps$age), animal = ps$sample),
  dc_rows("DC_frontal", dc_fr$meta),
  dc_rows("DC_temporal", dc_te$meta)
)
rownames(samples) <- NULL
cat(sprintf("\nsample table: %d samples\n", nrow(samples)))
print(table(samples$system))

saveRDS(list(
  iA = iA, iA_bc = iA_bc, iA_sym = ia_sym, iN = iN, BR_bc = BR_bc, BR_sym = BR_sym,
  tech_sd = tech_sd, pheno = pheno, br2sub = br2sub,
  ia_brain = ia_brain, in_brain = in_brain, in_ia = in_ia,
  O1p = O1p, O1p_sym = O1p_sym, s1p = s1p, org_meta = org_meta,
  M_L2 = M_L2, M_L2_24 = M_L2_24, M_PG = M_PG, o2o = o2o, sym2up = sym2up,
  DC_FR = dc_fr$mat, DC_TE = dc_te$mat, dc_fr_meta = dc_fr$meta, dc_te_meta = dc_te$meta,
  samples = samples), CACHE)
cat(sprintf("\nwrote %s\n", CACHE))
