# easy-MD — Easy Molecular Docking

**一键式分子对接 Claude Code Skill**，基于 AutoDock Vina 1.2.5 的全自动化流程。

[中文](#中文) | [English](#english)

---

## 中文

### 简介

easy-MD 是一个 Claude Code Skill，将分子对接全流程自动化：从靶标蛋白搜索、配体结构获取，到蛋白准备、对接执行、结果分析，只需一句话即可完成。

### 功能

- 🔬 自动搜索 PDB 蛋白结构（按分辨率排序，自动过滤非 X-RAY）
- 💊 自动搜索 PubChem 配体（SMILES 三级回退，解决新药数据缺失）
- 🧹 蛋白自动修复（pdbfixer：补缺失原子/残基、加氢 pH 7.4）
- 🎯 结合位点自动检测（基于共晶配体）
- ⚡ 三重独立对接（exhaustiveness=32, seeds 42/2024/8888）
- 📊 标准化输出（多模型 PDB + 分数表 + 完整报告）

### 依赖

#### 必需软件

| 软件 | 版本 | 说明 |
|------|------|------|
| Python | 3.12+ | 脚本执行环境 |
| AutoDock Vina | 1.2.5 | 分子对接引擎 |
| OpenBabel | 3.1.1 | 化学格式转换 |
| PyMOL | (可选) | 三维可视化 |

#### Python 库

```
rdkit pdbfixer openmm biopython numpy
```

安装：`pip install rdkit pdbfixer openmm biopython numpy`

### 安装

```bash
# 克隆到 Claude Code skills 目录
mkdir -p ~/.claude/skills
cd ~/.claude/skills
git clone https://github.com/<your-username>/easy-MD.git

# 重启 Claude Code 即可自动加载
```

### 使用

在 Claude Code 中直接说：

```
分子对接，ROCK2 与 Belumosudil，结果放桌面
```

Skill 自动触发，完成以下流程：

```
Workflow 0: 搜索蛋白 (RCSB PDB) + 配体 (PubChem) → 用户确认
Workflow 1: 蛋白准备 (download → pdbfixer → PDBQT)
Workflow 2: 配体准备 (SMILES → RDKit 3D → PDBQT)
Workflow 3: 对接盒子定义 (自动检测结合位点)
Workflow 4: Vina 对接 ×3 (parallel, exhaustiveness=32)
Workflow 5: 结果打包 (all_poses.pdb + scores.csv + summary.txt)
```

### 输出

```
{目标目录}/{蛋白}_{配体}/
├── all_poses.pdb         # 9个构象合并 (PyMOL state切换)
├── best_pose.pdb         # 最佳构象
├── scores.csv            # 结合能分数表
├── summary.txt           # 完整对接报告
├── protein_clean.pdb     # 处理后的蛋白
├── receptor.pdbqt        # Vina受体
├── ligand.pdbqt          # Vina配体
├── config.txt            # 对接参数
└── log.txt               # Vina完整日志
```

### MCP 工具

Skill 内置了 bio-server MCP 工具用于数据库搜索：

| 工具 | 数据库 | 功能 |
|------|--------|------|
| `search_pdb` | RCSB PDB | 蛋白结构搜索 |
| `search_ligand` | PubChem | 配体信息搜索 |
| `get_smiles` | PubChem | SMILES获取（三级回退） |
| `download_pdb` | RCSB PDB | PDB文件下载 |

### 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| Exhaustiveness | 32 | 搜索精度（发表级别） |
| Repeat | 3 | 独立运行次数 |
| Seeds | 42, 2024, 8888 | 随机种子 |
| num_modes | 9 | 输出构象数 |
| Energy range | 3 | 能量窗口 (kcal/mol) |

---

## English

### Overview

easy-MD is a Claude Code Skill that automates the complete molecular docking pipeline using AutoDock Vina 1.2.5. From target selection to publication-ready results — all in one command.

### Features

- 🔬 Auto PDB search (resolution-sorted, X-RAY filtered)
- 💊 Auto PubChem ligand search (3-level SMILES fallback)
- 🧹 Automatic protein preparation (pdbfixer)
- 🎯 Binding site auto-detection (co-crystal ligand based)
- ⚡ Triple independent docking (exhaustiveness=32)
- 📊 Standardized outputs (multi-model PDB + scores + report)

### Dependencies

| Software | Version |
|----------|---------|
| Python | 3.12+ |
| AutoDock Vina | 1.2.5 |
| OpenBabel | 3.1.1 |

Python packages: `rdkit pdbfixer openmm biopython numpy`

### Installation

```bash
mkdir -p ~/.claude/skills
cd ~/.claude/skills
git clone https://github.com/<your-username>/easy-MD.git
# Restart Claude Code
```

### Usage

Just say in Claude Code:

```
dock ROCK2 with Belumosudil
```

The skill triggers automatically and runs 6 workflows:

```
WF0: Protein + Ligand search → confirmation
WF1: Protein preparation (pdbfixer + OpenBabel)
WF2: Ligand preparation (RDKit + OpenBabel)
WF3: Docking box definition (auto-detect)
WF4: Vina docking ×3 (parallel, seeds 42/2024/8888)
WF5: Result packaging (all_poses.pdb + scores + summary)
```

### Output Structure

```
{target}/{protein}_{ligand}/
├── all_poses.pdb         # 9 poses merged
├── best_pose.pdb         # Best pose
├── scores.csv            # Affinity table
├── summary.txt           # Full report
├── protein_clean.pdb     # Prepared protein
├── receptor.pdbqt        # Vina receptor
├── ligand.pdbqt          # Vina ligand
├── config.txt            # Box config
└── log.txt               # Vina log
```

### License

MIT
