################################################################################
#
# Human brain, iPSC-derived cells, organoids and mice on the same axes
#
#
# Part A: where the variance sits in each system
# Part B: units per group needed to detect a fold change, at the SD one unit
#         of each system carries
# Part C: what each contrast shows against shuffled labels, and its agreement
#         with the human AD effect
#
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

ALPHA <- 0.05
TARGET <- 0.80
K_ORGANOIDS <- 5
FOLD_CHANGES <- c(1.05, 1.1, 1.2, 1.5)
N_SHOWN <- c(3, 5, 10)

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

variance_tbl <- read.csv("tables/variance_specs.csv", stringsAsFactors = FALSE)
variance_protein <- read.csv(gzfile("tables/variance_per_protein.csv.gz"), stringsAsFactors = FALSE)
effects_tbl <- read.csv("tables/effects_summary.csv", stringsAsFactors = FALSE)

################################################################################
# PART A: WHERE THE VARIANCE SITS
################################################################################

cat("=== PART A: SD of one protein, by source and system (log2) ===\n\n")

SYSTEM_NAMES <- c(
  "DC frontal" = "Human brain, frontal",
  "ROSMAP DLPFC" = "Human brain, DLPFC",
  "DC temporal" = "Human brain, temporal",
  "iA" = "iAstrocytes",
  "iN" = "iNeurons",
  "organoid" = "Cerebral organoids",
  "LOAD2" = "Mouse brain, LOAD2",
  "Peng" = "Mouse brain, Peng"
)
COMPONENTS <- c(
  "technical" = "Repeat measurement of one sample",
  "within genome, biological" = "Same genome, measurement noise removed",
  "within genome" = "Same genome (animals of a strain, organoids of a line)",
  "between genome" = "Different genomes (donors, donor lines)"
)
palette_component <- setNames(c("grey45", "#56B4E9", "#0072B2", "#D55E00"), COMPONENTS)
shape_component <- setNames(c(4, 1, 16, 17), COMPONENTS)

sd_tbl <- variance_tbl %>%
  filter(parameter == "sd_log2", condition %in% names(COMPONENTS)) %>%
  mutate(
    system = sprintf("%s (%s)", SYSTEM_NAMES[dataset], assay),
    component = factor(COMPONENTS[condition], levels = COMPONENTS)
  )
sd_tbl$system <- factor(sd_tbl$system,
                        levels = rev(unique(sd_tbl$system[order(match(sd_tbl$dataset, names(SYSTEM_NAMES)))])))

for (s in rev(levels(sd_tbl$system))) {
  r <- sd_tbl[sd_tbl$system == s, ]
  cat(sprintf("  %-38s %s\n", s,
              paste(sprintf("%s %.3f", r$condition, r$value), collapse = " | ")))
}

share <- variance_tbl[variance_tbl$parameter == "technical_share", ]
cat("\nShare of variance that is measurement noise\n")
for (i in seq_len(nrow(share))) {
  cat(sprintf("  %-14s %-26s %4.0f%%  (quartiles %2.0f%% to %2.0f%%)\n", share$dataset[i],
              share$condition[i], share$value[i] * 100, share$lower[i] * 100, share$upper[i] * 100))
}

pM1 <- ggplot(sd_tbl, aes(x = system, y = value, ymin = lower, ymax = upper,
                          color = component, shape = component)) +
  geom_linerange(position = position_dodge(width = 0.6), linewidth = 0.8) +
  geom_point(position = position_dodge(width = 0.6), size = 3, stroke = 1.2) +
  coord_flip() +
  scale_color_manual(values = palette_component) +
  scale_shape_manual(values = shape_component) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.03))) +
  guides(color = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
  labs(
    title = "Where the variance sits in each system",
    subtitle = "Median protein, bars span the quartiles across proteins",
    x = NULL,
    y = "SD of one protein (log2)",
    color = NULL,
    shape = NULL
  ) +
  theme_pub

tiff("pics/figM_variance.tiff", units = "in", res = 300, width = 10, height = 7)
print(pM1)
dev.off()

################################################################################
# PART B: UNITS PER GROUP NEEDED TO DETECT A FOLD CHANGE
################################################################################

cat("\n=== PART B: Units per group needed to detect a fold change ===\n\n")

power_two_group <- function(n, d, alpha = ALPHA) {
  df <- 2 * n - 2
  tcrit <- qt(1 - alpha / 2, df)
  ncp <- d * sqrt(n / 2)
  pt(tcrit, df, ncp = ncp, lower.tail = FALSE) + pt(-tcrit, df, ncp = ncp)
}

units_required <- function(d, target = TARGET, n_max = 20000) {
  n <- 2:n_max
  hit <- which(power_two_group(n, d) >= target)
  if (length(hit)) n[hit[1]] else NA_integer_
}

median_sd <- function(dataset, condition) {
  variance_tbl$value[variance_tbl$parameter == "sd_log2" & variance_tbl$dataset == dataset &
                       variance_tbl$condition == condition]
}

organoid_var <- variance_protein %>%
  filter(system == "organoid", component %in% c("within genome", "between genome")) %>%
  select(protein, component, var) %>%
  pivot_wider(names_from = component, values_from = var) %>%
  filter(is.finite(`within genome`), is.finite(`between genome`))
organoid_line_sd <- function(k) {
  median(sqrt(pmax(organoid_var$`between genome`, 0) + organoid_var$`within genome` / k))
}

UNITS <- data.frame(
  system = c("Human brain, one sample per donor (TMT)",
             "iAstrocytes, one sample per line (TMT)",
             "iNeurons, one sample per line (label-free)",
             "Organoids, one organoid per line (DIA)",
             sprintf("Organoids, %d organoids per line (DIA)", K_ORGANOIDS),
             "Mouse LOAD2, one animal (TMT)",
             "Mouse Peng, one animal (TMT)"),
  group = c("Human brain", "iPSC monoculture", "iPSC monoculture", "Organoids", "Organoids",
            "Mouse", "Mouse"),
  unit = c("donor", "donor line", "donor line", "donor line", "donor line", "animal", "animal"),
  sd = c(median_sd("DC frontal", "total"),
         median_sd("iA", "between genome"),
         median_sd("iN", "between genome"),
         organoid_line_sd(1),
         organoid_line_sd(K_ORGANOIDS),
         median_sd("LOAD2", "within genome"),
         median_sd("Peng", "within genome")),
  stringsAsFactors = FALSE
)
UNITS$system <- factor(UNITS$system, levels = UNITS$system)

ad <- effects_tbl[effects_tbl$contrast == "human: AD vs control", ]
AD_FC <- c("AD, median protein" = ad$median_fold_change, "AD, top tenth" = ad$q90_fold_change)

fc_table <- merge(UNITS, data.frame(fold_change = sort(c(FOLD_CHANGES, AD_FC)))) %>%
  mutate(
    d = log2(fold_change) / sd,
    units_per_group = sapply(d, units_required)
  )
for (n in N_SHOWN) fc_table[[sprintf("power_n%d", n)]] <- power_two_group(n, fc_table$d)
fc_table <- fc_table %>% arrange(system, fold_change)

fc_cols <- sort(unique(fc_table$fold_change))
cat(sprintf("  %-44s %-10s %6s  %s\n", "system", "unit", "SD",
            paste(sprintf("FC %-5.3f", fc_cols), collapse = " ")))
for (s in levels(UNITS$system)) {
  r <- fc_table[fc_table$system == s, ]
  cat(sprintf("  %-44s %-10s %6.3f  %s\n", s, r$unit[1], r$sd[1],
              paste(sprintf("%-8s", ifelse(is.na(r$units_per_group), ">20000", r$units_per_group)),
                    collapse = " ")))
}
cat(sprintf("\n  FC %.3f and %.3f are the median and 90th percentile protein in human AD cortex\n",
            AD_FC[1], AD_FC[2]))

cat("\nPower with five units per group\n")
for (s in levels(UNITS$system)) {
  r <- fc_table[fc_table$system == s, ]
  cat(sprintf("  %-44s %-10s %6s  %s\n", s, "", "",
              paste(sprintf("%-8s", sprintf("%.1f%%", r$power_n5 * 100)), collapse = " ")))
}

curve_grid <- merge(UNITS, data.frame(fold_change = 2^seq(log2(1.02), log2(2), length.out = 80))) %>%
  mutate(units_per_group = sapply(log2(fold_change) / sd, units_required))

palette_group <- c("Human brain" = "black", "iPSC monoculture" = "#0072B2",
                   "Organoids" = "#D55E00", "Mouse" = "#009E73")
color_system <- setNames(palette_group[UNITS$group], UNITS$system)
linetype_system <- setNames(c("solid", "solid", "dashed", "solid", "dashed", "solid", "dashed"),
                            UNITS$system)

pM2 <- ggplot(curve_grid, aes(x = fold_change, y = units_per_group, color = system, linetype = system)) +
  annotate("rect", xmin = 1.02, xmax = 2, ymin = 3, ymax = 5, fill = "grey50", alpha = 0.15) +
  geom_vline(xintercept = AD_FC, linetype = "dotted", color = "grey40") +
  geom_line(linewidth = 1) +
  scale_x_continuous(trans = "log10", breaks = c(1.05, 1.1, 1.2, 1.5, 2)) +
  scale_y_continuous(trans = "log10", breaks = c(2, 3, 5, 10, 20, 50, 100, 500, 2000, 10000),
                     labels = scales::comma) +
  scale_color_manual(values = color_system) +
  scale_linetype_manual(values = linetype_system) +
  guides(color = guide_legend(ncol = 2), linetype = guide_legend(ncol = 2)) +
  labs(
    title = "Units per group needed to detect a fold change",
    subtitle = paste0("80% power on one protein at the SD one unit carries. Grey band marks 3 to 5 per group.\n",
                      sprintf("Dotted lines mark the median and top-tenth protein in human AD cortex (%.2f and %.2f).",
                              AD_FC[1], AD_FC[2])),
    x = "Fold change (log scale)",
    y = "Units per group (log scale)",
    color = NULL,
    linetype = NULL
  ) +
  theme_pub +
  theme(legend.key.width = unit(1.6, "cm"))

tiff("pics/figM_units_by_fold_change.tiff", units = "in", res = 300, width = 10, height = 7.5)
print(pM2)
dev.off()

write.csv(fc_table %>% select(-group), "tables/systems_units_required.csv", row.names = FALSE)

################################################################################
# PART C: WHAT EACH CONTRAST SHOWS
################################################################################

cat("\n=== PART C: Observed effects against shuffled labels ===\n\n")

CONTRAST_SYSTEMS <- c("human" = "Human brain", "iA" = "iAstrocytes", "organoid" = "Organoids",
                      "LOAD2" = "Mouse, LOAD2", "Peng" = "Mouse, Peng")

contrasts <- effects_tbl %>%
  mutate(
    system = factor(CONTRAST_SYSTEMS[sub(":.*$", "", contrast)], levels = CONTRAST_SYSTEMS),
    label = sprintf("%s (%s)", sub("^[^:]*: ", "", contrast),
                    ifelse(n_group0 > 0, paste(n_group1, "v", n_group0), n_group1))
  )
contrasts$label <- factor(contrasts$label, levels = rev(unique(contrasts$label)))

cat(sprintf("  %-46s %-15s %-15s %-15s %-15s %8s %7s %7s\n", "contrast", "FC median", "FC top tenth",
            "SD median", "SD top tenth", "FDR 5%", "rho", "null95"))
for (i in seq_len(nrow(contrasts))) {
  r <- contrasts[i, ]
  cat(sprintf("  %-46s %-15s %-15s %-15s %-15s %8d %7s %7s\n",
              paste0(sub(":.*$", "", r$contrast), ": ", r$label),
              sprintf("%.3f / %.3f", r$median_fold_change, r$null_median_fold_change),
              sprintf("%.3f / %.3f", r$q90_fold_change, r$null_q90_fold_change),
              sprintf("%.2f / %.2f", r$median_abs_d, r$null_median),
              sprintf("%.2f / %.2f", r$q90_abs_d, r$null_q90),
              r$n_fdr05,
              ifelse(is.na(r$rho_human_ad), "", sprintf("%+.2f", r$rho_human_ad)),
              ifelse(is.na(r$rho_null_95), "", sprintf("%.2f", r$rho_null_95))))
}
cat("\n  Each pair is observed / shuffled labels\n")

STATS <- c("Median protein", "Top tenth of proteins")
KINDS <- c("Observed", "Shuffled labels")

pair_long <- function(obs_median, null_median, obs_q90, null_q90) {
  bind_rows(
    contrasts %>% transmute(system, label, stat = STATS[1], Observed = .data[[obs_median]],
                            `Shuffled labels` = .data[[null_median]]),
    contrasts %>% transmute(system, label, stat = STATS[2], Observed = .data[[obs_q90]],
                            `Shuffled labels` = .data[[null_q90]])
  ) %>%
    mutate(stat = factor(stat, levels = STATS))
}

pair_panel <- function(d, title, x_label, show_labels = FALSE) {
  points <- d %>%
    pivot_longer(all_of(KINDS), names_to = "kind", values_to = "value") %>%
    mutate(kind = factor(kind, levels = KINDS))
  layers <- lapply(seq_along(STATS), function(i) {
    nudge <- position_nudge(y = c(0.2, -0.2)[i])
    list(
      geom_linerange(data = d[d$stat == STATS[i], ],
                     aes(y = label, xmin = `Shuffled labels`, xmax = Observed),
                     color = "grey55", linewidth = 0.7, position = nudge),
      geom_point(data = points[points$stat == STATS[i], ],
                 aes(x = value, y = label, shape = stat, fill = kind),
                 color = "#B2182B", size = 2.6, stroke = 0.9, position = nudge)
    )
  })
  ggplot() +
    layers +
    facet_grid(system ~ ., scales = "free_y", space = "free_y") +
    scale_shape_manual(values = setNames(c(21, 24), STATS), drop = FALSE) +
    scale_fill_manual(values = setNames(c("#B2182B", "white"), KINDS), drop = FALSE) +
    guides(fill = guide_legend(override.aes = list(shape = 21)),
           shape = guide_legend(override.aes = list(fill = "#B2182B"))) +
    labs(title = title, x = x_label, y = NULL, shape = NULL, fill = NULL) +
    theme_pub +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      strip.text.y = element_blank(),
      axis.text.y = if (show_labels) element_text() else element_blank(),
      plot.margin = margin(8, 6, 8, 6)
    )
}

pC1 <- pair_panel(pair_long("median_fold_change", "null_median_fold_change",
                            "q90_fold_change", "null_q90_fold_change"),
                  "Fold change", "Absolute fold change", show_labels = TRUE)

pC2 <- pair_panel(pair_long("median_abs_d", "null_median", "q90_abs_d", "null_q90"),
                  "Effect in SD of the contrast's own unit", "|d|")

pC3 <- ggplot(contrasts, aes(y = label)) +
  geom_vline(xintercept = 0, color = "grey60") +
  geom_linerange(aes(xmin = -rho_null_95, xmax = rho_null_95), color = "grey80", linewidth = 4,
                 na.rm = TRUE) +
  geom_point(aes(x = rho_human_ad), shape = 21, fill = "#B2182B", color = "#B2182B", size = 2.8,
             na.rm = TRUE) +
  facet_grid(system ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(-0.5, 0.5), breaks = c(-0.4, 0, 0.4)) +
  labs(title = "Agreement with human AD", x = "Spearman rho", y = NULL) +
  theme_pub +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    strip.text.y = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(8, 6, 8, 6)
  )

pC4 <- ggplot(contrasts, aes(x = 0, y = label, label = scales::comma(n_fdr05))) +
  geom_text(size = 4) +
  facet_grid(system ~ ., scales = "free_y", space = "free_y") +
  labs(title = "FDR 5%", x = "proteins", y = NULL) +
  theme_pub +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    strip.text.y = element_text(angle = 0, hjust = 0),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(8, 6, 8, 6)
  )

pM3 <- pC1 + pC2 + pC3 + pC4 +
  plot_layout(widths = c(3, 3, 2.2, 0.8), guides = "collect") +
  plot_annotation(
    title = "What each contrast shows against shuffled labels",
    caption = "Agreement is the correlation of protein effects with the human AD effect. The grey band spans 95% of shuffled-label correlations.",
    theme = theme(plot.title = element_text(face = "bold", size = 16),
                  plot.caption = element_text(size = 10, color = "grey40", hjust = 0))
  ) &
  theme(legend.position = "bottom")

tiff("pics/figM_effects.tiff", units = "in", res = 300, width = 15, height = 8)
print(pM3)
dev.off()

write.csv(contrasts %>%
            transmute(system, contrast, unit, n_group1, n_group0, proteins,
                      median_fold_change, null_median_fold_change, q90_fold_change, null_q90_fold_change,
                      median_abs_d, null_median, q90_abs_d, null_q90, n_fdr05,
                      rho_human_ad, rho_null_95, p_concordance, p_min_attainable),
          "tables/systems_effects.csv", row.names = FALSE)
