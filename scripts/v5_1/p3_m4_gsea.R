# SAP-consistent GSEA: rank genes by per-contrast dPDUI (and omnibus-directed stat),
# test 4 annotation axes with fgsea. Shortening = negative dPDUI -> negative NES = axis genes
# shift toward shortening. Output: results/p3_m4_gsea.tsv
# Run: /home/taylor/miniconda3/envs/rlimma/bin/Rscript p3_m4_gsea.R
suppressMessages(library(fgsea))
MATRIX <- "/mnt/d/stroke_apa/results/p2_full_pdui_matrix.tsv"
OUTD <- "/mnt/d/stroke_apa/results"
mat <- read.delim(MATRIX, check.names = FALSE)
samples <- c("sham1","sham2","day1_rep1","day1_rep2","day3_rep1","day3_rep2",
             "day7_rep1","day7_rep2","day21_rep1","day21_rep2","day60_rep1","day60_rep2")
P <- mat[, samples]
sym <- sub(".*\\|", "", mat$Gene)
cat("events:", nrow(mat), "genes:", length(unique(sym)), "\n")

# gene-level mean dPDUI per contrast (mean of event dPDUIs mapped to gene)
sham <- rowMeans(P[, c("sham1","sham2")], na.rm = TRUE)
pathsets <- list(
  PAP_localized = sub("\n", "", readLines("/mnt/d/stroke_apa/results/p3_m4_set1_pap_enriched.txt")),
  endfoot_stroke_responsive = readLines("/mnt/d/stroke_apa/results/p3_m4_set2_endfoot_de.txt"),
  zone_cortex = readLines("/mnt/d/stroke_apa/results/p3_m4_set3_stroke_responsive_cortex.txt"),
  zone_whitematter = readLines("/mnt/d/stroke_apa/results/p3_m4_set3_stroke_responsive_whitematter.txt"))

rows <- data.frame()
run_gsea <- function(stats, tag) {
  stats <- stats[!is.na(stats)]
  stats <- stats[order(stats, decreasing = TRUE)]
  # collapse duplicate gene names by mean
  stats <- tapply(stats, names(stats), mean)
  for (nm in names(pathsets)) {
    set <- pathsets[[nm]]
    if (!any(set %in% names(stats))) next
    fg <- fgseaSimple(pathways = list(axis = set), stats = stats, minSize = 10, maxSize = Inf, nperm = 10000)
    fg$contrast <- tag; fg$axis <- nm
    rows <<- rbind(rows, fg[, c("contrast","axis","ES","NES","pval","padj","size")])
  }
}
for (tp in c("day1","day3","day7","day21","day60")) {
  tpv <- rowMeans(P[, paste0(tp, "_rep1"), drop = FALSE])
  tpv <- rowMeans(P[, paste0(tp, c("_rep1","_rep2"))], na.rm = TRUE)
  d <- tpv - sham
  gs <- tapply(d, sym, mean, na.rm = TRUE)
  run_gsea(gs, tp)
}
# omnibus-directed: gene-level max |dPDUI| with sign of strongest contrast
best <- rep(NA_real_, nrow(mat))
for (tp in c("day1","day3","day7","day21","day60")) {
  tpv <- rowMeans(P[, paste0(tp, c("_rep1","_rep2")), drop = FALSE], na.rm = TRUE)
  d <- tpv - sham
  take <- is.na(best) & !is.na(d)
  best[take] <- d[take]
  upd <- !is.na(d) & !is.na(best) & abs(d) > abs(best)
  best[upd] <- d[upd]
}
run_gsea(tapply(best, sym, mean, na.rm = TRUE), "omnibus_directed")

rows <- rows[order(rows$pval), ]
write.table(rows, file.path(OUTD, "p3_m4_gsea.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(rows, digits = 3)
cat("[DONE] fgsea\n")
