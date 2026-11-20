# 结果保存与目录说明

源码与计算结果分开保存。本工作站统一使用 `D:/Desktop/codex-result`：

|流程|默认根目录|配置|
|---|---|---|
|小分子对接|`D:/Desktop/codex-result/easy_MD`|根目录 `config.yaml` 的 `paths.output_root`|
|蛋白互作|`D:/Desktop/codex-result/easy_PPI`|`easy_ppi/config.json` 的 `paths.output_root`|

这些是工作站默认值，新电脑先改成实际绝对路径。easy-PPI 的 `EASY_PPI_OUTPUT_ROOT` 可覆盖配置中的根目录，`EASY_PPI_CONFIG` 可指定其他 JSON 配置，单次运行显式 `--output` 优先。配置缺失或路径错误时报错，不静默创建项目内结果。新任务使用独立目录，拒绝覆盖已有目录。

## 推荐组织

```text
codex-result/
  easy_MD/<蛋白_小分子>/<运行目录>/
  easy_PPI/
    <蛋白A_蛋白B_日期>/<运行目录>/
    benchmarks/<基准名称>/<运行目录>/
    _archive/<开发或兼容性试验>/<运行目录>/
```

分类目录是组织约定，不是 CLI 根据任务自动推断的标签；用 `--output` 指定该次新运行的完整目录。省略 `--output` 时，PPI 在配置根目录创建 UTC 时间和随机标识组成的独立目录。开发试验和科学基准应明确指定 `_archive` 或 `benchmarks` 下的输出。

## 文件用途

|文件或目录|用途|
|---|---|
|`raw/`|RCSB 原始结构、实体、物种和组装元数据|
|`prepared/`|处理后的结构、原始编号及链映射|
|`analysis/`|Arpeggio 接触、FreeSASA SASA/BSA、PRODIGY 经验预测和 CSV/JSON|
|`docking/`|LightDock 原生运行数据、全部最终候选 PDB、评分和聚类|
|`validation/`|参考比较、负对照、面积不变性等验证|
|`visualization/`|可编辑 PSE、PNG 预览和视图检查|
|`logs/`|各外部命令的标准输出和错误日志|
|`provenance/`|本次计算使用的源码和配置快照|
|`run.json`|参数、命令、时间及阶段状态|
|`manifest.json`|文件大小及 SHA-256，用于完整性校验|

实验复合物直接分析时不需要进行 LightDock 对接。候选数量由采样参数决定；例如 1BRS 三种子基准保存 7,500 个最终候选，不能只保留评分第一名。运行成功、文件完整性通过和科学验证通过是不同结论。

## 本工作站历史迁移

2026-11-19 已将项目内历史结果迁移到上述根目录，保留全部计算文件内容和失败记录，并逐文件校验哈希。项目 `results/` 与 `ppi_results/` 只保留迁移说明。逐目录说明在桌面结果根目录，本地 `references/result-location-map.json` 提供旧路径到新路径映射；该机器专用清单不随源码发布。

历史日志、评分表和原始清单中的旧绝对路径保留原样，以维持原始审计记录和校验和；当前位置以本地迁移索引和新的导航/重开脚本为准。不要直接运行仍含旧硬编码路径的归档调查脚本。

GitHub 只包含源码、配置、测试和文档，不上传运行结构、候选、PSE/PNG、压缩交付包、虚拟环境、第三方二进制或本机迁移清单。
