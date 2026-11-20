# easy-PPI v0.1 架构与验收方案

## 范围和隔离

独立 Python 包 `easy_ppi`，独立配置、依赖和结果目录。不导入、不修改 `scripts/` 中的 Vina、PLIP 或配体准备逻辑。不接入结构预测模型。实验结构的刚体对接仍然是预测，不能证明两个蛋白存在生物学互作。

GitHub 发布包含独立 easy-PPI 源码、配置、文档和测试；不包含本地计算产物或第三方二进制，不修改 Vina 流程或历史结果。默认结果位于 `D:/Desktop/codex-result/easy_PPI/<UTC时间>_<随机标识>/`，由独立 `config.json` 配置；指定结果目录必须不存在。禁止静默回退到当前工作目录。2026-11-19 之前的项目内结果按用途迁移至桌面，保留原文件内容和历史路径记录，并提供逐目录迁移索引。

## 数据流

1. **doctor**：记录平台、解释器、发行包版本与许可证、可执行程序、WSL 状态。每个阶段独立报告可用性。Linux 优先用于 LightDock/Arpeggio/FreeSASA；Windows 原生环境必须实际验证，不能因为 pip 安装成功就认定可用。
2. **resolve**：名称先用 MyGene 解析（无默认人源限制，物种由用户可选指定），再查 RCSB 实验条目；使用实体的 UniProt/序列和物种检查同名误匹配。无唯一实体或链时输出候选并停止自动选择。PDB ID 下载原始 mmCIF，保存 REST 元数据；本地文件检查实验方法和来源，缺少实验来源时要求补充证据。
3. **inspect / prepare**：Gemmi 解析 PDB/mmCIF，检查实验方法、分辨率、R-work/R-free、生物学组装、物种、链、模型、缺失残基、主链原子、替代构象、占有率和修饰残基。第一版不自动补环；缺失主链、非标准残基、关键辅因子等不能静默修补或丢弃。保留原始文件和 auth/label 链、残基、插入码到规范 A/B 的映射。替代构象按整残基占有率确定主构象；记录处理决策。pdb-tools 验证 PDB 格式。
4. **route**：先查两实体共同实验条目及生物学组装。共同条目必须具有所选实体、同一组装且实际有界面；不能把同一晶胞内任意两链当作复合物。若不存在，经检索记录后进入 LightDock。检索失败是未知，不能记作没有实验复合物。独立结构需检查构建体、序列覆盖和物种。显式 `redock` 仅用于基准测试。
5. **analyze**：三项独立后端。
   - PDBe Arpeggio 在清理后的 mmCIF 上添加氢并分类；只提取 A–B 接触，保留完整 JSON。`ionic` 输出标为离子接触（盐桥候选），不升级为实验盐桥。芳环面接触单独保留；三维虚线只绘制有明确原子端点的接触，不自行猜测氢键。
   - FreeSASA 在相同重原子、参数、坐标下计算 SASA(A)、SASA(B)、SASA(AB)。`buried_total = A+B-AB`；`interface_area = buried_total/2`，单位 Å²。分别报告两侧埋藏量和逐残基面积。交换链与平移/旋转不应改变面积；远离构象应接近零。
   - PRODIGY 指定 A/B 和温度，保存原始输出；ΔG（kcal/mol）及 Kd（M）明确标为经验预测，不作为实验亲和力或 docking 排名依据。
6. **dock**：LightDock 原生 setup → GSO → 候选生成；记录随机种子、步数、glowworm 数、全部 swarms、评分函数及刚体/ANM 设置。保存每个最终候选及所有中间原生文件。刚体版本优先，不自动套用需要 AMBER 参数的氢/辅因子流程。候选按 LightDock 评分排序，保留聚类和碰撞信息；参考结构不参与排序。
7. **validate**：只在有独立标注的实验参考时调用 DockQ；链映射显式、序列核对，保留原始输出与 fnat、iRMSD、LRMSD、DockQ。无参考时 `SKIPPED/no_experimental_reference` 且数值 null。DockQ <0.23 不正确；0.23–0.49 可接受；0.49–0.80 中等；≥0.80 高质量。
8. **visualize**：PyMOL 创建 Protein_A、Protein_B、Interface_Residues、Interactions 四个顶层项。界面残基来源为 FreeSASA 埋藏量和 Arpeggio；原子作用线来源仅 Arpeggio，以虚线显示。会话保存、重新打开并核对对象、原子数及端点，PNG 仅为预览。
9. **audit**：阶段状态分别为 SUCCESS/FAILED/BLOCKED/SKIPPED，完整成功必须包括后端输出的解析与物理检查。清单记录 SHA256、参数、命令 argv、工作目录、开始结束时间、退出码、stdout/stderr、工具版本；JSON 和 CSV 中缺失值为空/null，不能写零代替未计算。

## 软件和许可证核验来源

下列为上游公开许可证；实际安装版本另记在环境报告中。不重新分发第三方二进制或修改上游授权声明。

|软件|用途|上游许可证 / 注意事项|来源|
|---|---|---|---|
|Gemmi|PDB/mmCIF、组装和映射|MPL-2.0|https://github.com/project-gemmi/gemmi/blob/master/LICENSE.txt|
|pdb-tools|导出格式检查|Apache-2.0，包含第三方声明|https://github.com/haddocking/pdb-tools/blob/master/LICENSE|
|LightDock|蛋白对接|GPL-3.0|https://lightdock.org/|
|PDBe Arpeggio|接触分类|GPL-3.0；需要 OpenBabel|https://github.com/PDBeurope/arpeggio/blob/master/setup.py|
|FreeSASA|SASA/BSA|MIT（应核对实际包及附带依赖）|https://github.com/mittinatten/freesasa/blob/master/LICENSE|
|PRODIGY|经验亲和力预测|Apache-2.0（以安装版 LICENSE 为准）|https://github.com/haddocking/prodigy|
|DockQ v2|实验参考构象比较|MIT（以安装版 LICENSE 为准；旧版需另查）|https://github.com/wallnerlab/DockQ/blob/master/LICENSE|
|开源 PyMOL|会话和图片|上游自定义宽松授权；商业 PyMOL 授权另论|https://github.com/schrodinger/pymol-open-source/blob/master/LICENSE|

GPL 工具通过独立进程使用，源码不复制到 Vina 流程；这不构成对未来打包分发方式的许可证法律结论。正式发布前确定项目自身许可证并附第三方声明。

## 1BRS 两阶段测试和验收

### 第一阶段：实验界面

- 获取 RCSB 1BRS 原始 mmCIF、条目/实体/组装证据与校验和，检查 2.0 Å X-ray、Barnase/Barstar 身份及 Bacillus 物种信息（数值以现场数据为准）。枚举非对称单元 A/B/C 和 D/E/F 的候选界面，选用生物学组装支持的链对；不得把六链整体当一个二聚体。
- 分离重原子结构，验证序列和残基编号；报告末端缺失残基，不能默补。计算 A/B/AB 和逐残基 SASA、BSA；核对两侧埋藏量之和、总面积定义及参数。
- Arpeggio 应有跨链接触，检查典型带电残基区的原子身份及距离；数量依赖加氢和定义，不能硬编码期待数来制造通过。对照文献/结构检查热点（例如 Barnase Arg59/Arg83/Arg87、Barstar Asp35/Asp39，逐一按原始编号核对）。
- 独立 PyMOL 重开 PSE、核对四项及作用线端点；验证 PNG 可读。
- DockQ 自比接近 1 是工具正对照，不宣称对接成功。
- Barstar 平移 200 Å 的负对照：BSA 接近零、无跨链接触、DockQ 应明显降低。

### 第二阶段：分离再对接

- 从同一实验链对分离单体，是 bound redocking，不能外推 unbound 对接成功。分别居中并确定性随机旋转 ligand，移除原复合物刚体变换；保留变换矩阵。
- 无参考界面 restraints；预先固定种子（42/2024/8888）、评分函数、采样量及步数。快速小采样只验证工具链，不能代替全局搜索验收。
- 本工作站先运行 10 swarms × 20 glowworms × 10 步的兼容性测试，再预先固定 50 swarms × 50 glowworms × 100 步 × 3 种子的 1BRS 基准。这个采样量只定义本次基准，不能保证任意蛋白的全局搜索收敛；通过与否仍以 DockQ 和物理检查为准。默认通用 CLI 保留更高采样参数。
- 性能检查发现 DFIRE 和上游 FastDFIRE 在五个 pilot 构象上的评分差约 0.15–0.42，不满足数值等价。正式 v0.1 基准选择 LightDock 官方 FastDFIRE（实测约快 10 倍），独立开始三种子测试，不混用两种评分或宣称结果等价。较慢的 DFIRE 运行保留为未完成的性能试验。评分器选择不依赖实验参考或最终 DockQ。
- 统计各种子候选数、LightDock 排名、聚类、碰撞、top1/top5/top10 的 DockQ 和 fnat/iRMSD/LRMSD。报告 oracle-best 仅作为采样覆盖指标，禁止用 DockQ 最佳替换 score-top1。
- 最少一个候选 DockQ≥0.23 表示该采样找到了可接受构象；top10 命中率单独判定排名能力。没有命中应明确 FAIL，并保留全部候选和失败原因，不能宣称科学验证通过。
- 选择 LightDock 评分首位候选重新运行第一阶段分析，比较实验界面与预测界面；亲和力差异只作经验预测对照。

### 自动化验证

单元测试覆盖误选链、实验来源缺失、替代构象一致性、编号插入码、跨链过滤、空后端结果、无参考 DockQ、面积计算约定、不覆盖结果和参考泄漏。集成测试需要真实工具，缺工具应 BLOCKED。Vina 原测试与源码哈希前后核对，确保原流程不变。
