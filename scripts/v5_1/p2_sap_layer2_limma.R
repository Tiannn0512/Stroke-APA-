# SAP layer 2 (v1.1): event-level moderated statistics via limma (user-approved R exception, 2026-09-19).
# Design: ~0 + group (6 levels). eBayes(robust=TRUE, trend=TRUE) per audit 3.3.
# Missingness: contrast layer requires 4/4 valid PDUI (sham1,sham2,tp1,tp2) else not_testable; no imputation.
# Omnibus: complete-case events (12/12 valid), joint F-test across all group coefficients.
# Effect threshold |dPDUI| >= 0.1 evaluated on the LINEAR PDUI scale (logit only for the test).
# Multiple testing: PRIMARY = BH across the union of 5 contrast p-vectors (joint); per-contrast BH recorded as auxiliary.
# Replicate concordance: both replicates must move in the same direction as the mean.
# Run: /home/taylor/miniconda3/envs/rlimma/bin/Rscript p2_sap_layer2_limma.R
options(stringsAsFactors = FALSE)
EPS <- 0.01  # boundary-corrected logit, pre-frozen
MATRIX <- "/mnt/d/stroke_apa/results/p2_full_pdui_matrix.tsv"
OUTD <- "/mnt/d/stroke_apa/results"
library(limma)

mat <- read.delim(MATRIX, check.names = FALSE)
samples <- c("sham1","sham2","day1_rep1","day1_rep2","day3_rep1","day3_rep2",
             "day7_rep1","day7_rep2","day21_rep1","day21_rep2","day60_rep1","day60_rep2")
stopifnot(all(samples %in% colnames(mat)))
P <- mat[, samples]
P[] <- lapply(P, function(x) as.numeric(as.character(x)))
cat("events loaded:", nrow(mat), "\n")

logit_eps <- function(x) log((x * (1 - 2 * EPS) + EPS) / (1 - (x * (1 - 2 * EPS) + EPS)))
L <- as.data.frame(lapply(P, function(x) { xi <- suppressWarnings(as.numeric(x)); ifelse(is.na(xi), NA, logit_eps(xi)) }))
group_of <- function(s) sub("_rep[12]$", "", sub("^(sham)[12]$", "\\1", s))
grp <- sapply(samples, group_of)

TPS <- c("day1","day3","day7","day21","day60")
contrast_rows <- list()
for (tp in TPS) {
  cols <- c("sham1","sham2", paste0(tp,"_rep1"), paste0(tp,"_rep2"))
  sub <- L[, cols]
  ok <- rowSums(is.na(sub)) == 0
  cat(tp, ": testable events (4/4 valid) =", sum(ok), "\n")
  if (sum(ok) < 10) { contrast_rows[[tp]] <- NULL; next }
  y <- as.matrix(sub[ok, ])
  g <- factor(c("sham","sham", tp, tp), levels = c("sham", tp))
  design <- model.matrix(~0 + g); colnames(design) <- c("sham", tp)
  fit <- lmFit(y, design)
  cfit <- contrasts.fit(fit, makeContrasts(paste0(tp, "-sham"), levels = design))
  efit <- eBayes(cfit, robust = TRUE, trend = TRUE)
  tt <- topTable(efit, number = Inf, sort.by = "none")
  out <- data.frame(event = rownames(tt), contrast = tp,
                    logitFC = tt$logFC, t = tt$t, p_within = tt$P.Value, padj_within = tt$adj.P.Val,
                    dPDUI = (rowMeans(P[ok, c(paste0(tp,"_rep1"), paste0(tp,"_rep2")), drop = FALSE]) -
                             rowMeans(P[ok, c("sham1","sham2"), drop = FALSE])),
                    d_rep1 = as.numeric(P[ok, paste0(tp,"_rep1")]) - as.numeric(P[ok, "sham1"]),
                    d_rep2 = as.numeric(P[ok, paste0(tp,"_rep2")]) - as.numeric(P[ok, "sham2"]))
  contrast_rows[[tp]] <- out
}
allc <- do.call(rbind, contrast_rows)
allc$padj_joint <- p.adjust(allc$p_within, method = "BH")  # PRIMARY: joint BH across events x contrasts
allc$sig_joint <- allc$padj_joint < 0.1 & abs(allc$dPDUI) >= 0.1 &
  sign(allc$d_rep1) == sign(allc$dPDUI) & sign(allc$d_rep2) == sign(allc$dPDUI) & allc$dPDUI != 0
write.table(allc, file.path(OUTD, "p2_sap_layer2_events.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n[contrast summary] joint-BH FDR<0.1 & |dPDUI|>=0.1 & reps concordant:\n")
print(table(allc$contrast[allc$sig_joint], allc$dPDUI[allc$sig_joint] < 0))

# omnibus on complete cases
ok12 <- rowSums(is.na(L)) == 0
cat("omnibus complete-case events:", sum(ok12), "\n")
y12 <- as.matrix(L[ok12, ])
g12 <- factor(grp, levels = c("sham", TPS))
design12 <- model.matrix(~0 + g12); colnames(design12) <- c("sham", TPS)
fit12 <- lmFit(y12, design12)
e12 <- eBayes(fit12, robust = TRUE, trend = TRUE)
F <- topTable(e12, number = Inf, sort.by = "none", coef = 1:6)
om <- data.frame(event = rownames(F), F = F$F, p_omnibus = F$P.Value, padj_omnibus = F$adj.P.Val)
write.table(om, file.path(OUTD, "p2_sap_layer2_omnibus.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("omnibus padj<0.1:", sum(om$padj_omnibus < 0.1), "\n[DONE] layer2 limma\n")
