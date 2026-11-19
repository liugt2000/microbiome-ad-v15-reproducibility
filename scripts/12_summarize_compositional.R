# Summarize completed ANCOM-BC2 and LinDA per-cohort outputs after heterogeneous schemas.
suppressPackageStartupMessages(library(metafor))
REPO <- getwd()  # run from repository root
OUT <- file.path(REPO, "results", "sensitivity", "compositional")
cohorts <- c("beijing", "hangzhou", "lishui", "shanghai", "turkey")
files <- as.vector(outer(cohorts, c("ancombc2", "linda"), function(a,b) file.path(OUT, paste0(a,"_",b,".csv"))))
if (!all(file.exists(files))) stop("Missing per-cohort result files: ", paste(files[!file.exists(files)], collapse="; "))
all_results <- lapply(files, read.csv, stringsAsFactors=FALSE, check.names=FALSE)
cols <- unique(unlist(lapply(all_results, names)))
all_results <- lapply(all_results, function(x) {
  for (nm in setdiff(cols, names(x))) x[[nm]] <- NA
  x[cols]
})
res <- do.call(rbind, all_results)
write.csv(res, file.path(OUT,"all_cohort_compositional_results.csv"), row.names=FALSE, na="")
man <- do.call(rbind, lapply(split(res,res$cohort), function(x) data.frame(
  cohort=x$cohort[1], n_samples=x$n[1], n_genera_ANCOMBC2=sum(x$method=="ANCOMBC2"),
  n_genera_LinDA=sum(x$method=="LinDA"), stringsAsFactors=FALSE)))
write.csv(man,file.path(OUT,"analysis_manifest.csv"),row.names=FALSE)
core_names <- c("Blautia", "Romboutsia", "Butyricicoccus", "Coprococcus", "Roseburia",
  "Faecalibacterium", "Dialister", "Agathobacter", "Lachnospiraceae NK4A136 group",
  "Lachnospiraceae UCG-004", "Eubacterium", "Dorea", "Anaerostipes",
  "Pseudobutyrivibrio", "Anaerococcus", "Bifidobacterium", "Bacteroides",
  "Parabacteroides", "Alistipes", "Corynebacterium", "Clostridium", "Ruminococcus",
  "Lachnospiraceae", "Lactobacillus", "Collinsella", "Erysipelatoclostridium",
  "Enterococcus", "Limosilactobacillus", "UBA1819", "Incertae Sedis",
  "Family XIII AD3011 group", "Lachnospiraceae FCS020 group", "Lachnospira",
  "Haemophilus", "Parasutterella")
core <- res[res$genus %in% core_names & is.finite(res$effect) & is.finite(res$se) & res$se > 0, ]
meta_rows <- list(); failures <- list()
for (g in core_names) for (m in c("ANCOMBC2","LinDA")) {
  x <- core[core$genus==g & core$method==m,]
  if (nrow(x)<2) next
  fit <- tryCatch(rma.uni(yi=x$effect,sei=x$se,method="REML",test="knha"), error=function(e) e)
  if (inherits(fit,"error")) { failures[[paste(g,m)]] <- data.frame(genus=g,method=m,error=conditionMessage(fit)); next }
  meta_rows[[paste(g,m)]] <- data.frame(genus=g,method=m,k=nrow(x),pooled_effect=as.numeric(fit$b),
    se=fit$se,ci_low=fit$ci.lb,ci_high=fit$ci.ub,p=fit$pval,tau2=fit$tau2,I2=fit$I2,Q=fit$QE,Q_p=fit$QEp)
}
if (length(meta_rows)) {
  meta <- do.call(rbind,meta_rows); meta$fdr <- ave(meta$p,meta$method,FUN=function(z)p.adjust(z,"BH"))
  write.csv(meta,file.path(OUT,"core_genera_reml_hk.csv"),row.names=FALSE,na="")
}
if (length(failures)) write.csv(do.call(rbind,failures),file.path(OUT,"meta_fit_failures.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(OUT,"sessionInfo_summary.txt"))
cat("DONE:",nrow(res),"cohort-method-genus results;",length(meta_rows),"core meta rows\n")
