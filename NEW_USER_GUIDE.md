# 新用户首次使用指南

## 1. 下载代码

公开地址：<https://github.com/wdc19981006-cell/easy-MD>。

可以让智能体使用 Git 克隆，也可下载 [main 分支 ZIP](https://github.com/wdc19981006-cell/easy-MD/archive/refs/heads/main.zip) 后解压。在解压目录中找到与 `scripts/` 同级的 `SKILL.md`，这个目录就是工作流根目录。

`main` 下载链接跟随主分支更新，它不是固定版本。需要复现时记录下载时的提交 SHA、配置和运行元数据。使用 Git 的用户可通过 `git rev-parse HEAD` 查看提交 SHA。

从仓库目录开始工作，让智能体明确读取 `SKILL.md`，即可按其步骤执行。若要在其他项目中复用，再按所用客户端的技能安装方式，将整个目录安装为名为 `easy-md` 的技能；单独复制 `SKILL.md` 会缺少脚本和配置。

## 2. 第一次会检查和安装什么

这不是自带科学软件的便携安装包，也没有一键安装脚本。README 中的首次配置指令会让智能体检查、安装和配置，具体安装方法取决于电脑上已有的 Python 和科学软件。

| 组件 | 用途 | 检查要点 |
| --- | --- | --- |
| Python | 执行工作流 | 确认实际可执行文件，原工作站使用 3.12；避开 Windows 商店占位程序 |
| requirements.txt 中的包 | 结构处理与分析 | 在同一个工作流 Python 环境中安装并检查导入 |
| AutoDock Vina | 正式对接 | `paths.vina` 指向真实可执行文件，检查 `--version` |
| OpenBabel | 格式转换与准备 | `paths.obabel` 指向真实可执行文件，检查 `-V` |
| PLIP 与 OpenBabel Python 绑定 | 相互作用分析 | 工作流 Python 能导入 `plip.structure.preparation` 和 `openbabel.openbabel` |
| PyMOL 与其 Python | 三维会话与校验 | `paths.pymol_python` 指向能导入 `pymol` 和 `numpy` 的解释器 |

安装 OpenBabel 命令行程序，并不意味着工作流 Python 已能导入 OpenBabel 绑定。`pip install -r requirements.txt` 也不负责安装 Vina 和 PyMOL 外部程序；若某个锁定版本在当前平台无法安装，智能体应报告具体失败项，再选择兼容安装方案，不能声称环境已准备完成。

建议将依赖放入独立环境。README 的复制指令明确授权为本工作流安装依赖；如果用户只要求下载，先完成下载和检查，在获得安装授权后再安装。

## 3. 写入本机配置

编辑 `config.yaml` 中的 `paths`。其中原作者的 `C:/Python312/...`、`D:/Interaction/...`、`D:/Desktop/...` 都是示例，不应直接用于新电脑。

- `python`：安装了工作流依赖的 Python 可执行文件。
- `vina` / `obabel`：本机 Vina 和 OpenBabel 的绝对路径。
- `pymol`：本机 PyMOL 启动程序。
- `pymol_python`：能运行 PyMOL 模块的 Python，可与工作流 Python 不同。
- `output_root`：独立的结果输出目录，推荐使用英文路径，避免科学软件的路径兼容问题。
- `desktop` / `work_root`：按本机情况调整；每个实际任务使用独立目录。

在执行任务前保存配置，并将选定的结果目录告知用户。不要用系统目录、用户主目录或源码目录本身作为单次任务工作目录。

## 4. 检查环境

使用第 3 步确认的 Python 运行；下面均为需要替换占位符的示意命令，不要原样粘贴。

PowerShell：

```powershell
& "<本机Python绝对路径>" "<工作流根目录>/scripts/preflight.py" --config "<工作流根目录>/config.yaml" --json
```

Git Bash：

```bash
"<Python的Git-Bash路径>" "C:/实际路径/easy-md/scripts/preflight.py" --config "C:/实际路径/easy-md/config.yaml" --json
```

Git Bash 中 Python 的执行路径可写成 `/c/.../python.exe`；传给 Windows 程序的脚本和数据路径应使用 `C:/...` / `D:/...` 格式，不要使用 `/tmp/...`。

预检 JSON 的 `ok` 表示核心对接依赖检查通过。`optional` 单独报告 PLIP 和 PyMOL 的能力；它检查包是否可发现及程序文件是否存在，不能替代真实导入和运行检查。

智能体还应使用工作流 Python 检查 `from plip.structure.preparation import PDBComplex`、`from openbabel import openbabel`，并使用配置的 PyMOL Python 检查 `import pymol, numpy`。完整报告和三维会话依赖这些组件，缺少它们时应在正式对接前报告，避免计算结束后才发现分析失败。

## 5. 提交蛋白和小分子

例如：“请使用 easy-md，准备 ROCK2 与 Belumosudil 的分子对接。”

智能体按 `SKILL.md` 运行 `pipeline.py --config CONFIG prepare --protein ROCK2 --ligand Belumosudil`。`CONFIG` 等占位符应替换为真实绝对路径，并使用配置的 Python 调用根目录 `scripts/` 下的脚本。

准备阶段会规范化名称、获取结构、准备受体和配体，并生成对接盒。若可调用生物数据库 MCP 工具，则优先使用；没有这些工具时，使用仓库内的 REST 获取模块，不要求新用户先安装作者私有的 MCP 服务。

若名称有歧义，用户选择正确的基因或化合物；若没有可靠的口袋依据，则提供有实验依据的对接盒中心和尺寸，不能自动改成盲对接。

## 6. 确认一次具体计算

智能体显示蛋白、PDB、配体、位点依据、盒中心和尺寸、随机种子、重复次数、计算强度、输出目录、预计耗时以及适用限制。

用户审阅并批准这些参数后，智能体运行 `pipeline.py execute --run-dir RUN_DIR --confirm`。首次安装授权不代表已经批准某一次分子对接。

## 7. 自动计算与结果交付

默认使用 3 个随机种子运行 Vina，排序所有构象，然后进行 PLIP 分析、PyMOL 可视化、结构化打包与校验；重对接仅在选择该选项时执行。PLIP 是唯一相互作用分类来源。

每个任务下包含 `structure/`、`docking/`、`analysis/`、`visualization/`、`run_status.json` 和 `metadata.json`。交付时提供：

- Vina 评分和最佳构象。
- PLIP 相互作用文字、JSON、XML 报告。
- PyMOL `.pse` 会话和 `.pml` 脚本。
- 各阶段状态、校验结果、配置和实际输出目录。

`.pse` 会话中可分别开关 `Protein`、`Contact_Residues`、`Ligand`、`Interactions`。标准流程不生成二维相互作用图。

失败时读取 `run_status.json` 与对应日志，报告失败阶段和已有产物；不要直接重复执行整个流程或覆盖原目录。状态检查使用 `pipeline.py status --run-dir RUN_DIR`，结果检查使用 `pipeline.py validate --run-dir RUN_DIR`。

## 时间预期

首次安装时间与单次计算时间分别说明。已有环境的原工作站配置预计约 5–30 分钟，不能据此承诺所有电脑和组合均在 6 分钟完成；网络、CPU、盒尺寸、配体自由度和参数都会影响时间。当前未在全新 Windows 电脑上完成自动安装全过程验证。
