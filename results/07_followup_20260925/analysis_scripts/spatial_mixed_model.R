# spatial_mixed_model.R — Gate C 定位分析框架（预定义，无数据不填结果）
#
# 主检验（预先固定，见 spatial_analysis_plan.md）：
#   卒中状态(day7 vs sham) × 区室(PAP vs soma) 对 L 版本占比(logit) 的交互作用
# 动物 = 生物学单位；细胞/RNA 点嵌套于动物内。
#
# 输入：spatial_counts.tsv，列 =
#   animal_id, group(sham|day7), compartment(PAP|soma), batch,
#   n_L, n_S, n_M, total_expression, pap_area_um2, spots_total, probe_detect_pct
# （由影像定量脚本产出；本框架不生成数据）

suppressPackageStartupMessages({
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("需要 lme4；未安装时先 install.packages('lme4')")
})

run_main_model <- function(path) {
  d <- read.delim(path, stringsAsFactors = FALSE)
  d$L_prop <- (d$n_L + 0.5) / (d$n_L + d$n_S + d$n_M + 1)   # 加平滑防 0/1
  d$group  <- factor(d$group, levels = c("sham", "day7"))
  d$compartment <- factor(d$compartment, levels = c("soma", "PAP"))
  m <- lme4::glmer(cbind(n_L, n_S + n_M) ~ group * compartment +
                     total_expression + (1 | animal) + (1 | batch),
                   data = d, family = binomial)
  print(summary(m))
  it <- m@beta["groupday7:compartmentPAP"]
  message("交互作用估计与 CI 请以 confint(m) 报告；效应量以 OR 及其 CI 报告。")
  invisible(m)
}

# 敏感性1：动物级汇总后 Wilcoxon/置换
# 敏感性2：排除血管周围终足行后的复算（行已在影像脚本中剔除）
# 敏感性3：S 版本占比同模型复算（S 用 J_S 推断，需检出效率校正先行）

if (sys.nframe() == 0) {
  message("用法：run_main_model('spatial_counts.tsv')；数据到位前不运行。")
}
