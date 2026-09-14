# 猕猴年龄模板 Brainnetome 图谱：v4 预览版

本项目将成年猕猴 Brainnetome 皮层分区适配到六组年龄模板。提供两个独立附件：

- [v4 主图谱包](https://github.com/wenlii/Macaque-development-Brainnetome-Atlas/releases/download/v4-preview/macaque_age_atlas_v4_minimal.zip)：六组 32k 表面、原始 0.5 mm 模板空间体积、掩膜、QC 和验证脚本。六组均保留 248 个分区，**解剖 QC 状态仍为 HOLD**。
- [展示专用包](https://github.com/wenlii/Macaque-development-Brainnetome-Atlas/releases/download/v4-preview/macaque_age_atlas_v4_display_only.zip)：填洞、平滑后的外观标签、加密展示表面和图片，**仅供展示，不能用于后续计算**。
- [SHA-256 校验值](https://github.com/wenlii/Macaque-development-Brainnetome-Atlas/releases/download/v4-preview/SHA256SUMS.txt)、[发布清单](https://github.com/wenlii/Macaque-development-Brainnetome-Atlas/releases/download/v4-preview/release_manifest.json)。

请从 [Release 附件](https://github.com/wenlii/Macaque-development-Brainnetome-Atlas/releases/tag/v4-preview) 下载数据；GitHub 自动生成的 Source code ZIP 不是图谱数据包。

完整解压后，主图谱用每组 `atlas.wb.spec` 打开，展示版用 `DISPLAY_ONLY_dense.wb.spec` 打开。两套标签与掩膜不能混用。主图谱仍保留原有几何可疑位置；展示版的平滑和补洞不代表这些位置已经得到解剖修复。

[年龄与样本信息](metadata/age_groups.tsv)、[标签编号](metadata/volume_labels.tsv)、[逐组 QC](metadata/analysis_group_qc.tsv)、[逐 ROI QC](metadata/analysis_roi_qc.tsv)、[Group04 外观对比](figures/Group04_comparison_DISPLAY_ONLY.png) 和 [英文说明](README.md) 提供详细信息。

Group06 的实际年龄为 79–98 月龄，来自 3 只雄性动物的 7 次扫描；右 CG.RSr 在主图谱中仅占 8 个体素。分区数齐全不能替代独立解剖验证，模板量也不能视为个体测量均值。

作者、引用和现有许可说明见 CITATION.cff 与 LICENSE.md。两个 ZIP 中的本地构建回执保持原样，其“当时尚未上传”字段描述的是构建时状态；是否已公开发布以当前 GitHub Release 页面为准。
