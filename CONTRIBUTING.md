# 开发约定

唯一正式分享地址是 [easy-MD](https://github.com/wdc19981006-cell/easy-MD)（Public）。
开发、调试、回归测试和科学验证在 [easy-md-test](https://github.com/wdc19981006-cell/easy-md-test)（Private）进行。
旧 `easy-md-skill` 是同一私有仓库重命名前的地址，不再用于新 remote 或安装说明。

## 本地工作区

本工作站 `D:/CodexProjects/easy-md` 是开发工作区：`origin` 指向 easy-md-test，`release` 指向 easy-MD。
默认 push 目标是 origin。先 `git fetch origin`，从最新 origin/main 创建 `feature/<主题>` 或 `fix/<主题>`；用 PR 合入测试仓库 main。
拉取采用 fast-forward，遇到分叉先检查，不用 reset --hard 或 force push 覆盖工作。

每个新 checkout 都要用**实际 Python 可执行文件**运行 `tools/install_hooks.py`。
原工作站 Git Bash 示例：

```bash
/c/Python312/python.exe D:/CodexProjects/easy-md/tools/install_hooks.py
/c/Python312/python.exe D:/CodexProjects/easy-md/tools/repository_guard.py
/c/Python312/python.exe -B -m unittest discover -s tests -v
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 /d/CodexProjects/easy-md/.ppi-venv/Scripts/python.exe -B -m pytest easy_ppi/tests -q -p no:cacheprovider
```

安装器仅修改当前仓库 Git 配置，保存解释器路径，不修改全局设置，也不覆盖已有其他 hooksPath。
正式仓库独立 checkout 的 origin 应指向正式仓库，并安装同样 hooks；不要在它的 main 直接开发。

## 保存什么

- 源码、单元/集成测试、配置**示例**、依赖说明、用户文档。
- 测试仓库额外保存 `development/methods/` 方法、`development/records/` 脱敏运行记录、`development/benchmarks/` CSV/JSON 指标摘要、`development/templates/` 记录模板。
- 每份运行记录写明输入身份/PDB ID、来源提交、算法/依赖版本、实际命令（隐去用户名和敏感路径）、种子、参数、耗时、错误与修复、程序状态、文件验证和科学验证结论。
- 失败、BLOCKED、SKIPPED/null 必须如实保存。记录不等于发布批准。

## 不能提交的内容

PDB/mmCIF/CIF/PDBQT/SDF/MOL/MOL2、全部候选/轨迹、中间结构、PSE、结构渲染图、压缩包、第三方程序、虚拟环境、缓存、个人配置、原始日志、密钥和凭据均保留本地。
不要把结构改名为 .txt/.json、嵌入 Markdown/源码或通过 LFS 绕过规则。合成测试在临时目录生成小型结构，不分发真实结构 fixture。
原始日志先人工脱敏，再把必要片段保存为 Markdown 或 UTF-8 文本；自动扫描不能替代人工检查邮箱、样本/患者信息和机构信息。

本地 hooks 检查暂存区全部文件、待推送的新提交；CI 检查当前树与新增提交。它们会拒绝被 `git add -f` 强行暂存的禁传文件。
hooks 可被本地用户绕过，GitHub Actions 也不能阻止一次已经发生的 push；正式 main 还需配置 PR 和必需检查，见发布流程。现有历史的单独审查不会因 .gitignore 自动消失。

Vina 结果：`D:/Desktop/codex-result/easy_MD`；PPI 结果：`D:/Desktop/codex-result/easy_PPI`。
保持独立配置、依赖和结果目录；不因整理源码移动、覆盖、删除或重新计算历史结果。

科学计算仍遵循 SKILL.md：Vina 的具体计算批准与依赖安装批准分开；PPI 使用实验来源与 Arpeggio。
后续发布见 [RELEASE_PROCESS.md](RELEASE_PROCESS.md)。
