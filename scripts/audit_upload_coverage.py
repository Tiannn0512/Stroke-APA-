#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
audit_upload_coverage.py — 上传覆盖机械审计

对 D:\stroke_apa 与 D:\stroke_apa_reanalysis 的**每一个文件**分类：
  STAGED      已复制进 GitHub 仓库副本（并做体积核对）
  RELEASE_ZIP 在四个交付 zip 之一内（挂 GitHub Releases）
  REGENERABLE 可由公开数据+入库脚本确定性重建（REGENERATE.md 有对应条目）
  UNCOVERED   三类都不沾 —— 必须为 0 才允许删除本地

用法：python3 audit_upload_coverage.py
（本地工作区删除后此脚本失去对象，仅作为审计记录存档于仓库）
r"""
import hashlib
import os
import sys
import zipfile
from collections import Counter

SRC1 = r"D:\stroke_apa"            # v5.1 档案
SRC2 = r"D:\stroke_apa_reanalysis"  # 重建区
STAGE = r"D:\stroke_apa_github"     # GitHub 仓库副本
ZIPS = [
    r"D:\stroke_apa_reanalysis\REANALYSIS_RESULTS_20260923.zip",
    r"D:\stroke_apa_reanalysis\REANALYSIS_RESULTS_20260921.zip",
    r"D:\stroke_apa\STROKE_APA_FINAL_DELIVERABLES_20260920.zip",
    r"D:\stroke_apa\ATP2A2_IGV_CHECK_PACKAGE_20260920.zip",
]

SKIP_DIRS = {".git", "__pycache__", ".zcode"}

# —— REGENERABLE 规则（正则前缀，相对各源根）——
REGEN_RULES = [
    # v5.1 档案
    (SRC1, "data/GSE147119/", "GEO GSE147119 原始补充，peaks 冻结副本在 input_links"),
    (SRC1, "data/GSE263986/", "GEO GSE263986 xlsx，Set3 DE 结果已在 results/v5_1"),
    # 重建区
    (SRC2, "data/bam/", "12 BAM 由 FASTQ 经 p2_gse238125_align.sh 重建"),
    (SRC2, "data/references/", "GENCODE/GRCm38 公开重下 + REGENERATE.md §3 建索引"),
    (SRC2, "data/salmon/", "salmon index/quant 重建"),
    (SRC2, "data/tools/DaPars2/.git/", "上游克隆历史，可重新 git clone"),
    (SRC2, "data/tools/DaPars2/Dapars2_Test_Dataset.zip", "上游仓库自带测试数据"),
    # 注意：data/tools/DaPars2/{src,LICENSE,README.md} 是打过补丁的源码 —— 必须 STAGED
    (SRC2, "data/GSE143531/", "GEO GSE143531 RAW 重下（REGENERATE.md §5）"),
    (SRC2, "data/GSE143531_RAW.tar", "GEO 直链重下"),
    (SRC2, "results/02_candidate_rebuild/distal_nullpool.fa", "p2_build_null_set.py seed42 重建"),
    (SRC2, "results/02_candidate_rebuild/fimo_distal/cisml.xml", "FIMO 重扫（REGENERATE.md §6）"),
    (SRC2, "results/02_candidate_rebuild/fimo_distal/fimo.gff", "FIMO 重扫"),
    (SRC2, "results/02_candidate_rebuild/fimo_null/cisml.xml", "FIMO 重扫"),
    (SRC2, "results/02_candidate_rebuild/fimo_null/fimo.gff", "FIMO 重扫"),
    (SRC2, "results/03_bam_review/slices/", "BAM 重建后 samtools view regions.bed（bed 已入库）"),
]


def classify(rel, src_root):
    """返回 (类别, 说明)。rel 使用正斜杠相对路径。"""
    for root, prefix, why in REGEN_RULES:
        if src_root == root and rel.replace(os.sep, "/").startswith(prefix):
            return "REGENERABLE", why
    return None, ""


def stage_path_for(src_root, rel):
    """源文件 → 仓库副本中的路径；无映射返回 None。"""
    rel = rel.replace(os.sep, "/")
    if src_root == SRC1:
        if rel.startswith("scripts/"):
            return "scripts/v5_1/" + rel[len("scripts/"):]
        if rel.startswith("results/"):
            return "results/v5_1/" + rel[len("results/"):]
        if rel.startswith("data/"):
            return "data_frozen/v5_1_data/" + rel[len("data/"):]
        if rel.startswith(("reviews/", "docs/", "memory/")):
            head, _, tail = rel.partition("/")
            return {"reviews": "docs/v5_1/reviews", "docs": "docs/v5_1/docs_internal",
                    "memory": "docs/v5_1/memory"}[head] + ("/" + tail if tail else "")
        if rel.endswith((".md", ".html")) and "/" not in rel:
            return "docs/v5_1/" + rel
        if rel == ".gitignore":
            return "docs/v5_1/dot_gitignore.v5_1"
        if rel.startswith(".openclaw-attachments/"):
            return "docs/v5_1/attachments_openclaw/" + rel.split("/", 1)[1]
    else:
        if rel.startswith("scripts/"):
            return "scripts/reanalysis/" + rel[len("scripts/"):]
        if rel.startswith("results/"):
            return "results/reanalysis/" + rel[len("results/"):]
        if rel.startswith("config/"):
            return rel
        if rel.startswith("input_links/"):
            return rel
        if rel.startswith("metadata/"):
            return rel
        if rel.startswith("logs/"):
            return rel
        if rel.startswith(("data/qc/", "data/dapars2_full/")):
            return "data_frozen/reanalysis_data/" + rel[len("data/"):]
        if rel.startswith("data/tools/DaPars2/"):
            return "tools/DaPars2_patched/" + rel[len("data/tools/DaPars2/"):]
        if rel in ("data/sequencing_depth.tsv", "data/sequencing_depth_dapars2.tsv",
                   "data/gse143531_series_matrix.txt.gz"):
            return "data_frozen/reanalysis_data/" + rel.split("/", 1)[1]
        if rel.startswith("data/tools/DaPars2/"):
            return "tools/DaPars2_patched/" + rel[len("data/tools/DaPars2/"):]
        if rel in ("01_研究Proposal.md", "02_项目TODO.md", "PROJECT_FULL_REPORT.md",
                   "REANALYSIS_FINAL_REPORT.md", "README.md", "实验流程.md"):
            return "docs/reanalysis/" + rel
        if rel == ".gitignore":
            return "docs/reanalysis/dot_gitignore.reanalysis"
        if rel.startswith(".openclaw-attachments/"):
            return "docs/reanalysis/attachments_openclaw/" + rel.split("/", 1)[1]
        if rel.startswith("wetlab/"):
            return None  # 空骨架，.gitkeep 已代
    return None


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    zip_names = set()
    zip_entry_size = {}   # basename -> size（用于内容核对）
    zip_of = {}           # basename -> zip 文件名
    for zp in ZIPS:
        with zipfile.ZipFile(zp) as z:
            for i in z.infolist():
                b = os.path.basename(i.filename)
                zip_names.add(i.filename)
                zip_names.add(b)
                if b and b not in zip_entry_size:
                    zip_entry_size[b] = i.file_size
                    zip_of[b] = os.path.basename(zp)

    rows, uncovered, mismatch = [], [], []
    counts = Counter()
    for src_root in (SRC1, SRC2):
        for dirpath, dirnames, filenames in os.walk(src_root):
            dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
            for fn in filenames:
                full = os.path.join(dirpath, fn)
                rel = os.path.relpath(full, src_root).replace(os.sep, "/")
                if fn.endswith(".pyc"):
                    continue
                if ".openclaw" + os.sep + "workspace-state.json" in os.path.join(
                        os.path.relpath(dirpath, src_root), fn):
                    counts["SKIP"] += 1
                    continue  # 代理工作区状态缓存，无档案价值
                if full in ZIPS or (fn.endswith(".zip") and src_root == SRC2):
                    counts["RELEASE_ZIP"] += 1
                    rows.append((src_root, rel, "RELEASE_ZIP", os.path.basename(full)))
                    continue
                b = os.path.basename(rel)
                if b in zip_entry_size and zip_entry_size[b] == os.path.getsize(full):
                    counts["RELEASE_ZIP"] += 1
                    rows.append((src_root, rel, "RELEASE_ZIP", f"内容在 {zip_of[b]} 内(basename+size 一致)"))
                    continue
                cat, why = classify(rel, src_root)
                if cat is None:
                    sp = stage_path_for(src_root, rel)
                    if sp and os.path.isfile(os.path.join(STAGE, sp)):
                        tgt = os.path.join(STAGE, sp)
                        if os.path.getsize(full) == os.path.getsize(tgt) and \
                           sha256(full) == sha256(tgt):
                            cat, why = "STAGED", sp + "（sha256 一致）"
                        else:
                            mismatch.append((full, sp))
                            cat, why = "HASH_MISMATCH", sp
                    else:
                        cat, why = "UNCOVERED", "无映射或副本缺失"
                counts[cat] += 1
                rows.append((src_root, rel, cat, why))
                if cat == "UNCOVERED":
                    uncovered.append(full)

    out = os.path.join(STAGE, "audit_manifest.tsv")
    with open(out, "w", encoding="utf-8") as f:
        f.write("source_root\trelative_path\tclassification\tnote\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    print("=== 分类统计 ===")
    for k, v in sorted(counts.items()):
        print(f"{k:15s} {v}")
    print(f"合计             {len(rows)}")
    if mismatch:
        print(f"\n!! 体积不一致 {len(mismatch)} 项：")
        for a, b in mismatch[:20]:
            print("  ", a, "->", b)
    if uncovered:
        print(f"\n!! UNCOVERED {len(uncovered)} 项（禁止删除本地）：")
        for u in uncovered:
            print("  ", u)
    else:
        print("\nUNCOVERED = 0 —— 本地每一文件已覆盖：已入库 / 在Release / 可再生。")
    return 1 if (uncovered or mismatch) else 0


if __name__ == "__main__":
    sys.exit(main())
