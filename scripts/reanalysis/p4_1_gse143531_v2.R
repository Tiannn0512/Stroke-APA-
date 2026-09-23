#!/usr/bin/env Rscript
# Task A (03_执行方Agent_下一步任务书 §4): GSE143531 PAP vs Full, edgeR TMM.
# Unit of analysis = biological LIBRARY (unpaired 3 PAP vs 3 Full; PAP libraries
# pool 2 mice each per Mazaré 2020 methods, so C1/C2/C3 labels alone do not
# justify a paired model). Technical files summed per gene in the input table.
.libPaths(c("~/Rlibs", .libPaths()))
suppressPackageStartupMessages({library(edgeR)})

BASE <- "/mnt/d/stroke_apa_reanalysis"
OUT <- file.path(BASE, "results/04_public_pap")
cnt <- read.delim(file.path(OUT, "counts_by_biological_library.tsv"),
                  row.names = 1, check.names = FALSE)
qcmeta <- read.delim(file.path(OUT, "library_qc.tsv"), stringsAsFactors = FALSE)
qcmeta$library_id <- as.character(qcmeta$library_id)
stopifnot(ncol(cnt) == nrow(qcmeta))
stopifnot(identical(colnames(cnt), qcmeta$library_id))

meta <- qcmeta
meta$group <- factor(ifelse(meta$cell_part == "PAP", "PAP", "Full"), levels = c("Full", "PAP"))
cat("== design ==\n")
print(data.frame(library = colnames(cnt), group = meta$group, lib_size = meta$total_counts))

y <- DGEList(counts = cnt, group = meta$group)
keep <- filterByExpr(y, group = meta$group)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
cat(sprintf("kept genes after filterByExpr: %d / %d\n", nrow(y), nrow(cnt)))
print(y$samples[, c("group", "lib.size", "norm.factors")])
write.table(y$samples[, c("group", "lib.size", "norm.factors")],
            file.path(OUT, "tmm_samples.tsv"), sep = "\t")

design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- levels(meta$group)
write.table(design, file.path(OUT, "design_matrix.tsv"), sep = "\t", quote = FALSE)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit, contrast = makeContrasts(PAP - Full, levels = design))
tab <- topTags(qlf, n = Inf, sort.by = "none")$table
write.csv(data.frame(gene_id = rownames(tab), tab, check.names = FALSE),
          file.path(OUT, "edgeR_TMM_PAP_vs_Full.tsv"), row.names = FALSE)

# ---- candidate summary v2 (normalized CPM + edgeR test) ----
sym <- read.delim(file.path(BASE, "input_links/dapars2_event_coordinates.tsv"),
                  stringsAsFactors = FALSE)
sym$ensembl <- sub(".*\\|(.*)\\|.*", "\\1", sym$event_id)
sym2gene <- unique(sym[, c("gene_symbol", "ensembl")])
sym2gene <- sym2gene[!duplicated(sym2gene$gene_symbol), ]
cpmm <- cpm(y)
libs <- colnames(cpmm)
cands <- c("Atp2a2","Agpat3","Aplp1","Sirt2","Ndrg2","Cnp","Pea15a","Kazn","Plec","Fam107a")
rows <- list()
for (g in cands) {
  gid <- sym2gene$ensembl[sym2gene$gene_symbol == g]
  if (length(gid) != 1 || !gid %in% rownames(tab)) {
    rows[[length(rows)+1]] <- data.frame(gene_symbol = g, ensembl = gid,
                                         status = "filtered_out_or_missing")
    next
  }
  t <- tab[gid, ]
  cv <- as.numeric(cpmm[gid, ]); names(cv) <- libs
  getc <- function(lib) round(cv[paste0("X", lib) %in% names(cv) |
                                  lib %in% names(cv)], 2)
  # column names keep original library ids (check.names=FALSE on input)
  raw <- as.numeric(cpmm[gid, ])
  names(raw) <- libs
  rows[[length(rows)+1]] <- data.frame(
    gene_symbol = g, ensembl = gid, status = "ok",
    edgeR_log2FC_PAP_vs_Full = round(t$logFC, 3), edgeR_FDR = signif(t$FDR, 3),
    edgeR_logCPM = round(t$logCPM, 2),
    C1_Full_CPM = round(raw["2017545"], 2), C1_PAP_CPM = round(raw["2017547"], 2),
    C2_Full_CPM = round(raw["2017553"], 2), C2_PAP_CPM = round(raw["2017555"], 2),
    C3_Full_CPM = round(raw["2017561"], 2), C3_PAP_CPM = round(raw["2017563"], 2))
}
cand <- do.call(rbind, rows)
write.table(cand, file.path(OUT, "candidate_gene_summary.v2.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- descriptive genome percentiles from NORMALIZED edgeR logFC ----
ok <- is.finite(tab$logFC)
bg <- sort(tab$logFC[ok])
out <- list()
for (g in cands) {
  gid <- sym2gene$ensembl[sym2gene$gene_symbol == g]
  if (length(gid) == 1 && gid %in% rownames(tab)) {
    out[[length(out)+1]] <- data.frame(gene_symbol = g,
      edgeR_log2FC = round(tab[gid, "logFC"], 3),
      descriptive_genome_percentile = round(100 * mean(bg <= tab[gid, "logFC"]), 1))
  }
}
res <- do.call(rbind, out)
write.table(res, file.path(OUT, "candidate_percentiles_normalized.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("== descriptive percentiles (edgeR logFC; NOT a significance test) ==\n")
print(res)

# ---- MDS ----
if (requireNamespace("png", quietly = TRUE) || TRUE) {
  ok_dev <- tryCatch({ png(file.path(OUT, "library_mds.png"), width = 1200,
                           height = 900, res = 150); TRUE }, error = function(e) FALSE)
  if (ok_dev) {
    plotMDS(y, labels = colnames(y),
            col = ifelse(meta$group == "PAP", "#e45756", "#4c78a8"),
            main = "GSE143531 biological libraries (edgeR TMM) - MDS")
    dev.off()
    cat("library_mds.png written\n")
  }
}

cat("== sessionInfo ==\n"); print(sessionInfo())
