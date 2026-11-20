# easy-MD / easy-PPI

面向 Windows 的蛋白质–小分子自动化对接工作流，以及独立的实验结构蛋白–蛋白互作模块，供 Claude Code、Codex 等可访问本地文件并执行命令的智能体使用。

下方首次使用流程针对 easy-MD 小分子对接；蛋白互作请阅读 [easy-PPI 使用说明](easy_ppi/README.zh-CN.md)。

## 公开下载

- 仓库：<https://github.com/wdc19981006-cell/easy-MD>
- 最新主分支 ZIP：[点击下载](https://github.com/wdc19981006-cell/easy-MD/archive/refs/heads/main.zip)
- Git：`git clone https://github.com/wdc19981006-cell/easy-MD.git easy-md`

本仓库根目录就是 Skill 根目录，包含 `SKILL.md`、`config.yaml`、`scripts/` 和测试。无需再进入旧版的内层 `easy-MD/` 文件夹。

## 新用户：复制这段给 Claude Code / Codex

```text
请从 https://github.com/wdc19981006-cell/easy-MD.git 下载最新 easy-MD 到当前工作区的 easy-md 文件夹。
先阅读 README.md、NEW_USER_GUIDE.md 和 SKILL.md，按说明完成首次配置。
我同意为这个工作流在独立环境安装缺少的依赖，并联网从官方来源获取所需软件。
请检查本机实际可用的 Python、Vina、OpenBabel、PLIP 和 PyMOL，
把实际程序路径和一个独立的结果输出目录写入配置，然后运行环境检查和必要的导入检查。
不要沿用原作者电脑的路径，不要覆盖已有配置或结果；已有同名目录时先检查内容。
本次只完成安装和检查，完成后告诉我哪些功能可用，以及下一步如何提交蛋白和小分子名称。
```

首次安装完成后，再输入：

```text
请使用 easy-md，准备 ROCK2 与 Belumosudil 的分子对接。
先展示蛋白、PDB、配体、结合位点依据、对接盒、重复次数、随机种子、计算强度、输出目录和预计耗时，等我确认后再运行。
```

智能体在准备阶段后暂停。用户批准本次具体参数后，才执行正式对接与后续分析。仅发送仓库 URL 不能保证智能体完成下载、安装和配置，推荐使用上面的完整指令。详细步骤见 [NEW_USER_GUIDE.md](NEW_USER_GUIDE.md)。

## 下载后的流程

1. **获取代码**：克隆仓库或解压 ZIP，并读取根目录说明。
2. **首次配置**：检查 Windows 环境，安装缺少的依赖，修改 `config.yaml` 的本机路径及输出目录。
3. **检查环境**：执行 `scripts/preflight.py`，另行验证 PLIP/OpenBabel Python 绑定和 PyMOL 导入。缺失项会报告给用户。
4. **提交任务**：输入蛋白和小分子名称；获取结构、准备受体/配体并生成对接盒。
5. **确认参数**：审阅结构、结合位点与计算参数，明确批准这一次计算。
6. **自动计算与交付**：默认运行 3 个随机种子的 Vina 对接、构象排序、PLIP 分析、PyMOL 可视化和结果校验。

没有足够结合位点依据或名称存在歧义时，需要用户提供信息后继续。

## 时间与环境要求

目前主要面向 Windows，原工作站使用 Python 3.12。需要 Vina、OpenBabel、`requirements.txt` 中的 Python 包，以及用于完整分析和可视化的 PLIP、PyMOL。软件安装包未包含在本仓库中。

`config.yaml` 保留原工作站的路径示例，**新电脑必须改为本机真实路径**，包括 `python`、`vina`、`obabel`、`pymol`、`pymol_python`、`output_root`，并按需要修改 `desktop` / `work_root`。

首次下载、安装和配置的耗时不计入单次对接耗时。已有配置的原工作站示例预计约 5–30 分钟，实际取决于 CPU、结构、配体、对接盒及参数；不能保证所有任务都在 6 分钟完成。尚未在一台全新的 Windows 电脑上验证自动安装全过程。

## 功能与结果

- 蛋白名称与基因符号规范化、人源 PDB 结构筛选、PubChem 配体获取。
- 受体与配体准备、对接盒生成、AutoDock Vina 多种子对接与构象排序。
- 可选原配体重对接验证。
- 仅使用 PLIP 进行相互作用分类，生成 XML、JSON 和文字报告。
- PyMOL 会话包含 `Protein`、`Contact_Residues`、`Ligand`、`Interactions` 四个可独立开关的对象。
- 运行状态记录、结构化输出和保存会话的几何校验。

结果写入配置的 `paths.output_root`。主要文件包括 `docking/scores.csv`、`structure/best_pose.pdb`、`analysis/interaction_summary.txt` 和 `visualization/interaction_visualization.pse`。标准流程不生成二维相互作用图。

可信的已有本地结果可用 `scripts/repackage.py OLD_RUN` 重新整理并生成简化 PyMOL 会话，无需重跑 Vina 或 PLIP，原目录保留。

## 文档与验证

- [SKILL.md](SKILL.md)：智能体工作流。
- [NEW_USER_GUIDE.md](NEW_USER_GUIDE.md)：首次安装、使用与常见问题。
- [Windows 环境规则](references/environment-windows.md)。
- [命令参考](references/workflow-commands.md)。
- [结构选择和验证](references/selection-and-validation.md)。
- [输出规范](references/outputs.md)。
- `tests/`：单元和集成测试；使用本机已配置的 Python 执行 `-B -m unittest discover -s tests -v`。

`SKILL.md` 与 `scripts/` 针对蛋白质–小分子对接；蛋白质–蛋白质互作使用独立的 `easy_ppi/`。两者均不执行分子动力学。标准 Vina 不模拟共价键形成；对接分数不能直接代表实验亲和力、药效或作用机制。

## easy-PPI：实验结构蛋白互作

只接受有实验来源的 PDB/mmCIF，不接入 AlphaFold、Boltz 或其他序列到结构预测模型。输入可为蛋白名称、PDB ID 或本地文件；检查物种、链、结构质量、缺失残基和生物学组装。已有可用实验复合物时直接分析实际界面；只有独立实验结构时才使用 LightDock 对接。

- Gemmi/pdb-tools 准备结构并保留原始文件和残基编号映射。
- PDBe Arpeggio 分类跨链接触；FreeSASA 计算 SASA/BSA；PRODIGY 提供明确标注的经验亲和力预测。
- 有实验参考时用 DockQ 验证；无参考时跳过，数值为 null。
- PyMOL 会话分层显示 Protein_A、Protein_B、Interface_Residues、Interactions，作用线使用虚线。
- 保留全部候选、参数、日志、评分及 CSV/JSON/PDB/PSE/PNG，不用程序无报错替代科学验证。

依赖使用独立环境，安装步骤和 Windows 兼容性见 [easy-PPI 使用说明](easy_ppi/README.zh-CN.md)；设计、软件选择、许可证来源和两阶段验收见 [架构文档](easy_ppi/ARCHITECTURE.zh-CN.md)。

**当前验证边界：**1BRS 实验界面分析通过；7,500 候选的重新对接未达到 DockQ≥0.23，保留 FAIL，不能宣称对接性能已通过验证。NRF2–KEAP1 使用 2FLU 的 NRF2 短肽–KEAP1 Kelch 实验片段分析，不代表全长构象。详见 [验证记录](easy_ppi/VALIDATION.zh-CN.md)。

## 默认结果位置

|流程|配置|本工作站默认输出|
|---|---|---|
|easy-MD 小分子对接|根目录 config.yaml → paths.output_root|D:/Desktop/codex-result/easy_MD|
|easy-PPI 蛋白互作|easy_ppi/config.json → paths.output_root|D:/Desktop/codex-result/easy_PPI|

新电脑须修改为本机绝对路径。PPI 支持 EASY_PPI_OUTPUT_ROOT / EASY_PPI_CONFIG；显式 --output 优先。配置错误时不会静默回退到项目 ppi_results。源码、计算结果、开发归档和基准数据分开保存，详见 [结果保存说明](references/result-storage.md)。
