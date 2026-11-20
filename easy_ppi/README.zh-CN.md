# easy-PPI v0.1

easy-MD 旁的独立蛋白–蛋白模块，只接受具有实验来源的 PDB/mmCIF。与原 `scripts/` 中的 Vina/PLIP 流程、依赖和配置独立。设计、软件许可证来源和验收标准见 [ARCHITECTURE.zh-CN.md](ARCHITECTURE.zh-CN.md)。

## 能力

- 名称 → MyGene 尝试基因规范化 → RCSB 实验候选 → 实体描述/基因、UniProt 和物种核对。名称歧义保留候选 JSON 并停止，不默认人源。融合构建体不自动选用。
- 支持 PDB ID、本地 PDB/mmCIF；检查来源、模型、分辨率、链、物种、缺失残基/主链/侧链、占有率、B-factor、修饰残基及生物学组装。无实验来源、主链缺失或非标准残基会阻断自动准备；未知物种明确标记。
- 已有实验生物学复合物优先直接分析。独立结构按 UniProt/精确序列检索共同实验条目，核对物种和构建体；检索失败、复杂复制链或构建体不兼容会停止，不能假称“没有复合物”。
- 名称输入没有指定某个构建体，因此可以选择同一 UniProt/物种的实际实验复合物，并在 `route.json/construct_review` 公开其与初始结构候选的序列差异与归档突变数；不默认宣称野生型。明确指定 PDB ID/文件时不会这样替换构建体。
- Gemmi 保留原始编号及 auth/label 映射；pdb-tools 检查格式。保留观察到的侧链，报告缺失原子，不自动补环、补侧链或使用结构预测。
- OpenBabel 原生 PDB 加氢，检查重原子身份和坐标不变、氢数量与命名唯一性；PDBe Arpeggio 负责全部作用分类。离子接触写作盐桥候选；不将几何距离自行判为氢键。
- FreeSASA 三次计算 A/B/AB，并给出逐残基 SASA、两侧埋藏量、总埋藏量与其一半。PRODIGY 明确标为经验预测；严重重叠构象不解释亲和力。
- LightDock 0.9.4：完整记录 setup/GSO/候选生成/原生聚类；默认使用其官方 FastDFIRE，可显式选 DFIRE。保存每个最终 glowworm 的 PDB。不同评分器结果不混排。
- DockQ 仅在有实验参考时运行。无参考写 `SKIPPED` 和 null，不产生伪分数。保留原生 CLI 输出，批量评价使用上游 API 全精度。按 LightDock 分数排名，DockQ 只做评价。
- PyMOL 的四个顶层项为 `Protein_A`、`Protein_B`、`Interface_Residues`、`Interactions`；作用线为虚线，只取 Arpeggio 的明确原子端点。平面/环接触保存在 JSON，未猜测原子端点。保存 PSE 后重开核对原子数、作用线的保存端点与虚线设置，PNG 另作预览。

## 环境

推荐独立 Linux/WSL2 环境：Python 3.11、OpenBabel、开源 PyMOL，并安装本目录 `requirements.txt`。不要把依赖装入 Vina 的环境。可使用本目录的 `environment-linux.yml`；先切换到本目录创建 Conda 环境，使其中的相对 requirements 路径可解析。WSL/Conda 安装本身不会由模块自动执行。

本工作站实际测试了 Windows 11 + Python 3.12 + Rtools GCC。普通 `pip install` 因缺少 MSVC 失败，已通过项目内部 `.ppi-venv` 和 `.ppi-tools` 构建，不修改全局 Python。当前测试环境继承原科学包，只在虚拟环境固定 NumPy 1.26.4；全局 NumPy 保持原版本。新机器建议使用完全独立环境。

Windows 的 LightDock 兼容处理均有记录：

1. 用已有 Rtools GCC 构建 FreeSASA/ProDy/LightDock。GCC 14+ 对旧 NumPy C 指针转换的诊断降为警告；未禁用运行时断言。
2. PyPI sdist 未包含 `.pyx`。从固定 release commit `eaf05c7c39a3d8fbf715cb74a2dfe2857778107f` 获取原始 `.pyx`，保存哈希，再用当前 Cython 重建 CPython 3.12 扩展。旧 `.c` 保留为 `.c.pre_cython`。
3. `readfp → read_file` 仅在科学子进程中处理 Python 配置 API 更名。
4. 上游 Windows 文本写入产生 CR-CR-LF。保留原始字节 gzip，逐文件记录 SHA256，只规范化初始化和 GSO 文件换行。
5. DockQ 使用上游不修改的官方 `operations_nocy` 脚本模式，避免缺失 MSVC。原生 CLI 的 NumPy JSON 序列化问题通过保存 CLI 文本和调用上游批量 API 处理；未重写 DockQ 公式。

开发构建辅助工具：`python -m easy_ppi.download_source NAME VERSION`、`python -m easy_ppi.fetch_lightdock_pyx`、`python easy_ppi/build_mingw.py ABSOLUTE_SOURCE_DIR`、`python -m easy_ppi.bootstrap_dockq`。源码下载目录必须尚不存在；工具保存来源与校验和。必须在项目根目录运行，并用已配置 GCC 的 PATH。PyMOL 可以使用独立解释器；本机使用已有 PyMOL 安装，不分发其二进制或授权文件。

## Git Bash 命令

下列均在 `D:/CodexProjects/easy-md` 运行。`python` 是说明性通称；本机实际使用指定解释器：

```bash
/d/CodexProjects/easy-md/.ppi-venv/Scripts/python.exe -m easy_ppi doctor

/d/CodexProjects/easy-md/.ppi-venv/Scripts/python.exe -m easy_ppi run \
  --protein-a 1BRS --chain-a A --protein-b 1BRS --chain-b D \
  --assembly-id 1 --pymol-python D:/Interaction/PyMol/python.exe

# 名称有歧义时输出候选 JSON；可改用明确的 PDB ID/链。
/d/CodexProjects/easy-md/.ppi-venv/Scripts/python.exe -m easy_ppi run \
  --protein-a Barnase --species-a 'Bacillus amyloliquefaciens' \
  --protein-b Barstar --species-b 'Bacillus amyloliquefaciens' \
  --pymol-python D:/Interaction/PyMol/python.exe

# 本地文件仍须具有实验方法/来源；未知来源的纯坐标不会被认作实验结构。
/d/CodexProjects/easy-md/.ppi-venv/Scripts/python.exe -m easy_ppi inspect D:/path/experimental.cif

# 仅在使用上游脚本模式时设置；正常 DockQ 安装无需此环境变量。
export EASY_PPI_DOCKQ_SCRIPT=D:/CodexProjects/easy-md/.ppi-tools/dockq-2.1.3/DockQ.py

/d/CodexProjects/easy-md/.ppi-venv/Scripts/python.exe -m easy_ppi benchmark \
  --seeds 42 2024 8888 --swarms 50 --glowworms 50 --steps 100 \
  --scoring fastdfire --cores 4 --pymol-python D:/Interaction/PyMol/python.exe

/d/CodexProjects/easy-md/.ppi-venv/Scripts/python.exe -m easy_ppi validate D:/path/to/new_run

PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 /d/CodexProjects/easy-md/.ppi-venv/Scripts/python.exe -m pytest easy_ppi/tests -q
```

`--output` 必须是不存在的新目录；默认保存到 **`D:/Desktop/codex-result/easy_PPI/<UTC时间>_<随机标识>/`**，由 `easy_ppi/config.json` 的 `paths.output_root` 控制，环境变量 `EASY_PPI_OUTPUT_ROOT` 可覆盖该根目录，`EASY_PPI_CONFIG` 可指定其他配置文件。默认根目录必须是绝对路径；配置缺失或错误时直接报错，不回退到项目目录。原 Vina 结果使用根目录 `config.yaml`，保存到 `D:/Desktop/codex-result/easy_MD`。

目录含 `raw/ prepared/ analysis/ docking/ validation/ visualization/ logs/`、`run.json`、`environment.json` 和 `manifest.json`。`validate` 只读核对清单及哈希，明确分开文件完整性和科学有效性。2026-11-19 之前的开发结果已迁移至桌面结果根目录；本工作站旧路径与新路径对照保存在本地 `references/result-location-map.json`（不上传）；通用目录规则见 [结果保存说明](../references/result-storage.md)。历史日志中的旧路径保持原样以保留审计证据。

已有外部实验参考时可添加 `--reference ID_OR_FILE --reference-chain-a CHAIN --reference-chain-b CHAIN`。参考先核对实验来源与组装，再交给 DockQ；从不进入 setup、优化或评分排序。生物学组装的复制链可用 `A@1`、`A@2` 显式选择；只有两个相同链拷贝时可自动成对。更复杂的多聚体、非标准残基和辅因子依赖界面需要人工审阅，不应把模块当作通用自动修复器。

## 科学边界

1BRS 实验界面分析和从同一复合物分离的 bound redocking 是两个独立测试。自比 DockQ=1 只是工具正对照，不是重新对接成功。小采样可正常生成 PDB 而 DockQ 全部失败；必须如实报告。三种子基准也不等于搜索收敛或任意蛋白上通用有效。

`sampling_validation` 表示是否找到 DockQ≥0.23 的候选；`ranking_validation` 表示评分 top10 是否命中；还需审阅严重重叠和界面化学。`oracle_best_for_sampling_only` 不能冒充实际预测的 top1。1BRS 原论文描述 Barstar C40A/C82A，不能只依据归档的 mutation 字段宣称野生型。实验结构缺失残基、侧链、水介导相互作用和推断质子化都限制结果解释。PRODIGY 不提供实验亲和力、疗效或致病证据。

模块状态、失败记录和所有候选都会保留；未通过科学验证的对接不能因为脚本无异常而标记成功。通用 CLI 返回 2 表示某阶段失败/受阻，基准返回 3 表示计算完成但科学验收未通过。本版本已完成测试的范围和失败结论见 [VALIDATION.zh-CN.md](VALIDATION.zh-CN.md)。GitHub 仅保存源码、配置、测试和说明；原始结构、全部候选、日志、PSE/PNG 与迁移清单保留在本地结果目录。
