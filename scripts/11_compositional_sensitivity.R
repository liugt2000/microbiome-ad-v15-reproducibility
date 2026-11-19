# Public-cohort differential-abundance sensitivity analysis
# Inputs: previously rebuilt ASV-to-genus raw count matrices and harmonized group labels.
# Methods: ANCOM-BC2 and LinDA, independently rerun without changing source files.

suppressPackageStartupMessages({
  library(ANCOMBC)
  library(MicrobiomeStat)
})

REPO <- getwd()  # run from repository root
PREP <- file.path(REPO, "data", "compositional_prepared")
OUT <- file.path(REPO, "results", "sensitivity", "compositional")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
cohorts <- c("beijing", "hangzhou", "lishui", "shanghai", "turkey")
all_results <- list()
manifest <- list()
set.seed(20261119)

for (ch in cohorts) {
  message("\n=== ", ch, " ===")
  d <- readRDS(file.path(PREP, paste0(ch, "_prepared.rds")))
  genus <- d$genus
  meta <- d$meta
  rownames(meta) <- meta$sample_id
  common <- intersect(rownames(genus), rownames(meta))
  genus <- genus[common, , drop = FALSE]
  meta <- meta[common, , drop = FALSE]
  meta$group <- factor(meta$case_control, levels = c("HC", "AD"))
  feat <- t(round(genus))
  keep <- rowMeans(feat > 0) >= 0.10
  feat_keep <- feat[keep, , drop = FALSE]
  manifest[[ch]] <- data.frame(
    cohort = ch,
    n_samples = ncol(feat),
    n_AD = sum(meta$group == "AD"),
    n_HC = sum(meta$group == "HC"),
    n_genera_input = nrow(feat),
    n_genera_prevalence_10pct = sum(keep),
    median_library_size = median(colSums(feat)),
    stringsAsFactors = FALSE
  )

  message("ANCOM-BC2: ", nrow(feat), " input genera")
  ancom <- ancombc2(
    data = feat,
    taxa_are_rows = TRUE,
    meta_data = meta,
    fix_formula = "group",
    p_adj_method = "BH",
    prv_cut = 0.10,
    lib_cut = 0,
    s0_perc = 0.05,
    group = "group",
    struc_zero = FALSE,
    neg_lb = FALSE,
    alpha = 0.05,
    n_cl = 1,
    verbose = FALSE,
    global = FALSE,
    pairwise = FALSE
  )
  saveRDS(ancom, file.path(OUT, paste0(ch, "_ancombc2_full.rds")))
  ar <- ancom$res
  aout <- data.frame(
    cohort = ch,
    genus = ar$taxon,
    method = "ANCOMBC2",
    n = ncol(feat),
    effect = ar$lfc_groupAD,
    se = ar$se_groupAD,
    W = ar$W_groupAD,
    p = ar$p_groupAD,
    fdr = ar$q_groupAD,
    diff_abundant = ar$diff_groupAD,
    passed_ss = ar$passed_ss_groupAD,
    stringsAsFactors = FALSE
  )
  write.csv(aout, file.path(OUT, paste0(ch, "_ancombc2.csv")), row.names = FALSE, na = "")
  all_results[[paste0(ch, "_ancom")]] <- aout

  message("LinDA: ", nrow(feat_keep), " prevalence-filtered genera")
  linda_fit <- linda(
    feature.dat = feat_keep,
    meta.dat = meta,
    formula = "~ group",
    feature.dat.type = "count",
    prev.filter = 0,
    mean.abund.filter = 0,
    max.abund.filter = 0,
    is.winsor = TRUE,
    outlier.pct = 0.03,
    adaptive = TRUE,
    zero.handling = "pseudo-count",
    pseudo.cnt = 0.5,
    corr.cut = 0.1,
    p.adj.method = "BH",
    alpha = 0.05,
    n.cores = 1,
    verbose = FALSE
  )
  saveRDS(linda_fit, file.path(OUT, paste0(ch, "_linda_full.rds")))
  lo <- linda_fit$output
  if (is.list(lo) && !is.data.frame(lo)) {
    candidate_names <- names(lo)
    target <- candidate_names[grepl("group", candidate_names, ignore.case = TRUE)]
    if (length(target) == 0) target <- candidate_names[1]
    lr <- lo[[target[1]]]
  } else {
    lr <- lo
  }
  lr$genus <- rownames(lr)
  stat_col <- intersect(c("stat", "t", "z"), colnames(lr))[1]
  lout <- data.frame(
    cohort = ch,
    genus = lr$genus,
    method = "LinDA",
    n = ncol(feat_keep),
    effect = lr$log2FoldChange,
    se = if ("lfcSE" %in% colnames(lr)) lr$lfcSE else NA_real_,
    statistic = if (!is.na(stat_col)) lr[[stat_col]] else NA_real_,
    p = lr$pvalue,
    fdr = lr$padj,
    stringsAsFactors = FALSE
  )
  write.csv(lout, file.path(OUT, paste0(ch, "_linda.csv")), row.names = FALSE, na = "")
  all_results[[paste0(ch, "_linda")]] <- lout
  rm(ancom, linda_fit, d, genus, feat, feat_keep)
  gc(verbose = FALSE)
}

res <- do.call(rbind, all_results)
man <- do.call(rbind, manifest)
write.csv(res, file.path(OUT, "all_cohort_compositional_results.csv"), row.names = FALSE, na = "")
write.csv(man, file.path(OUT, "analysis_manifest.csv"), row.names = FALSE)

# Core results and robust random-effects synthesis with REML + Hartung-Knapp.
core_names <- c(
  "Blautia", "Romboutsia", "Butyricicoccus", "Coprococcus", "Roseburia",
  "Faecalibacterium", "Dialister", "Agathobacter", "Lachnospiraceae NK4A136 group",
  "Lachnospiraceae UCG-004", "Eubacterium", "Dorea", "Anaerostipes",
  "Pseudobutyrivibrio", "Anaerococcus", "Bifidobacterium", "Bacteroides",
  "Parabacteroides", "Alistipes", "Corynebacterium", "Clostridium", "Ruminococcus",
  "Lachnospiraceae", "Lactobacillus", "Collinsella", "Erysipelatoclostridium",
  "Enterococcus", "Limosilactobacillus", "UBA1819", "Incertae Sedis",
  "Family XIII AD3011 group", "Lachnospiraceae FCS020 group", "Lachnospira",
  "Haemophilus", "Parasutterella"
)
suppressPackageStartupMessages(library(metafor))
core <- res[res$genus %in% core_names & is.finite(res$effect) & is.finite(res$se) & res$se > 0, ]
meta_rows <- list()
for (g in core_names) {
  for (m in c("ANCOMBC2", "LinDA")) {
    x <- core[core$genus == g & core$method == m, ]
    if (nrow(x) < 2) next
    fit <- rma.uni(yi = x$effect, sei = x$se, method = "REML", test = "knha")
    meta_rows[[paste(g, m)]] <- data.frame(
      genus = g,
      method = m,
      k = nrow(x),
      pooled_effect = as.numeric(fit$b),
      se = fit$se,
      ci_low = fit$ci.lb,
      ci_high = fit$ci.ub,
      p = fit$pval,
      tau2 = fit$tau2,
      I2 = fit$I2,
      Q = fit$QE,
      Q_p = fit$QEp,
      stringsAsFactors = FALSE
    )
  }
}
meta <- do.call(rbind, meta_rows)
meta$fdr <- ave(meta$p, meta$method, FUN = function(z) p.adjust(z, "BH"))
write.csv(meta, file.path(OUT, "core_genera_reml_hk.csv"), row.names = FALSE, na = "")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("\nDONE: ", nrow(res), " cohort-method-genus results; ", nrow(meta), " core meta rows\n", sep = "")
