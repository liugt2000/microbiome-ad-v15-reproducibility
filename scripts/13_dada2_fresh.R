# Fresh raw-read DADA2 rerun for the 84-sample discovery cohort.
# Input: vendor raw FASTQ only. Output: independent validation directory.
suppressPackageStartupMessages({library(dada2); library(Biostrings)})
set.seed(42)

RAW <- Sys.getenv("FASTQ_DIR", unset="")
REF <- Sys.getenv("SILVA_TRAIN", unset="")
OUT <- Sys.getenv("DADA2_OUT", unset=file.path(getwd(), "results", "raw_reanalysis", "dada2_fresh"))
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

stopifnot(dir.exists(RAW), file.exists(REF))
fnFs <- sort(list.files(RAW, pattern="R1.raw.fastq.gz$", full.names=TRUE))
fnRs <- sort(list.files(RAW, pattern="R2.raw.fastq.gz$", full.names=TRUE))
stopifnot(length(fnFs)==84, length(fnRs)==84)
sample.names <- sub(".338F_806R.R1.raw.fastq.gz$", "", basename(fnFs))
stopifnot(identical(sample.names, sub(".338F_806R.R2.raw.fastq.gz$", "", basename(fnRs))))

manifest <- data.frame(
  sample=rep(sample.names, 2),
  direction=rep(c("R1","R2"), each=length(sample.names)),
  file=c(basename(fnFs), basename(fnRs)),
  bytes=c(file.info(fnFs)$size, file.info(fnRs)$size),
  md5=c(unname(tools::md5sum(fnFs)), unname(tools::md5sum(fnRs))),
  stringsAsFactors=FALSE
)
write.csv(manifest, file.path(OUT, "input_fastq_manifest.csv"), row.names=FALSE)
ref_manifest <- data.frame(file=basename(REF), bytes=file.info(REF)$size,
  md5=unname(tools::md5sum(REF)), stringsAsFactors=FALSE)
write.csv(ref_manifest, file.path(OUT, "reference_manifest.csv"), row.names=FALSE)

filtFs <- file.path(OUT, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(OUT, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names; names(filtRs) <- sample.names

cat(format(Sys.time()), "FILTER_START\n")
filter.out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
  truncLen=c(280,220), maxN=0, maxEE=c(2,2), truncQ=2,
  rm.phix=TRUE, compress=TRUE, multithread=12)
filter.track <- data.frame(sample=sample.names, reads.in=filter.out[,"reads.in"],
  reads.filtered=filter.out[,"reads.out"], stringsAsFactors=FALSE)
write.csv(filter.track, file.path(OUT, "filter_track.csv"), row.names=FALSE)
saveRDS(filter.out, file.path(OUT, "filter_out.rds"))

cat(format(Sys.time()), "LEARN_ERRORS_START\n")
errF <- learnErrors(filtFs, multithread=12, randomize=TRUE)
saveRDS(errF, file.path(OUT, "errF.rds"))
errR <- learnErrors(filtRs, multithread=12, randomize=TRUE)
saveRDS(errR, file.path(OUT, "errR.rds"))

cat(format(Sys.time()), "DADA_START\n")
dadaFs <- dada(filtFs, err=errF, multithread=12, pool=FALSE)
saveRDS(dadaFs, file.path(OUT, "dadaF.rds"))
dadaRs <- dada(filtRs, err=errR, multithread=12, pool=FALSE)
saveRDS(dadaRs, file.path(OUT, "dadaR.rds"))

cat(format(Sys.time()), "MERGE_START\n")
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
saveRDS(mergers, file.path(OUT, "mergers.rds"))
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, file.path(OUT, "seqtab_all.rds"))
cat("ASVs before length/abundance filters:", ncol(seqtab), "\n")
seqtab <- seqtab[, nchar(getSequences(seqtab)) >= 250 & nchar(getSequences(seqtab)) <= 500, drop=FALSE]
seqtab <- seqtab[, colSums(seqtab) >= 5, drop=FALSE]
saveRDS(seqtab, file.path(OUT, "seqtab_prefiltered.rds"))
cat("ASVs after filters:", ncol(seqtab), "\n")

cat(format(Sys.time()), "CHIMERA_START\n")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=12, verbose=TRUE)
saveRDS(seqtab.nochim, file.path(OUT, "seqtab_nochim.rds"))
cat("Final dimensions:", dim(seqtab.nochim), "reads retained:", sum(seqtab.nochim)/sum(seqtab), "\n")

cat(format(Sys.time()), "TAXONOMY_START\n")
taxa <- assignTaxonomy(seqtab.nochim, REF, multithread=12, minBoot=80, tryRC=TRUE)
saveRDS(taxa, file.path(OUT, "taxa.rds"))

asv_sequences <- colnames(seqtab.nochim)
writeXStringSet(DNAStringSet(asv_sequences), file.path(OUT, "asv_sequences.fasta"))
write.csv(data.frame(sample=rownames(seqtab.nochim), nonchim_reads=rowSums(seqtab.nochim)),
  file.path(OUT, "nonchim_reads.csv"), row.names=FALSE)

getN <- function(x) sum(getUniques(x))
merged_reads <- sapply(mergers, getN)
track <- cbind(filter.out, merged=merged_reads,
  nonchim=rowSums(seqtab.nochim)[sample.names])
track <- as.data.frame(track)
track$sample <- sample.names
track$pct_retained_vs_input <- round(100*track$nonchim/track$reads.in, 2)
write.csv(track, file.path(OUT, "read_tracking_final.csv"), row.names=FALSE)

genus <- ifelse(is.na(taxa[,"Genus"]), "Unclassified", taxa[,"Genus"])
genus.count <- rowsum(t(seqtab.nochim), group=genus)
genus.count.t <- t(genus.count)
genus.relabund <- genus.count.t / rowSums(genus.count.t)
write.csv(data.frame(sample=rownames(genus.relabund), genus.relabund, check.names=FALSE),
  file.path(OUT, "genus_table_dada2.csv"), row.names=FALSE)
write.csv(data.frame(sample=rownames(genus.count.t), genus.count.t, check.names=FALSE),
  file.path(OUT, "genus_table_dada2_count.csv"), row.names=FALSE)

blautia <- data.frame(sample=rownames(seqtab.nochim), Group=sub("_.*$", "", rownames(seqtab.nochim)),
  total_reads=rowSums(seqtab.nochim), Blautia_count=seqtab.nochim[, genus=="Blautia", drop=FALSE] |> rowSums(),
  Blautia_relabund=genus.relabund[,"Blautia"], stringsAsFactors=FALSE)
write.csv(blautia, file.path(OUT, "blautia_per_sample_dada2.csv"), row.names=FALSE)

summary <- data.frame(
  metric=c("samples","asv_before_filters","asv_after_filters","asv_nonchim",
    "input_reads","filtered_reads","nonchim_reads","mean_pct_retained","genera"),
  value=c(nrow(seqtab.nochim), ncol(readRDS(file.path(OUT,"seqtab_all.rds"))),
    ncol(seqtab), ncol(seqtab.nochim), sum(filter.out[,"reads.in"]),
    sum(filter.out[,"reads.out"]), sum(seqtab.nochim),
    mean(track$pct_retained_vs_input), ncol(genus.relabund))
)
write.csv(summary, file.path(OUT, "run_summary.csv"), row.names=FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat(format(Sys.time()), "DONE\n")
