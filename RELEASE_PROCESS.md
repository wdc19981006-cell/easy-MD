# 从开发到正式发布

## 仓库职责

| 仓库 | 可见性 | 用途 |
| --- | --- | --- |
| [easy-MD](https://github.com/wdc19981006-cell/easy-MD) | Public | 唯一正式 Skill、稳定源码、依赖/配置示例、用户指南、版本和发布证据 |
| [easy-md-test](https://github.com/wdc19981006-cell/easy-md-test) | Private | 功能开发、兼容性调试、测试、脱敏方法/错误记录、科学验证与指标摘要 |

仓库名不代表算法验证已通过。easy-PPI 当前 1BRS 重新对接的采样与排名结论仍为 FAIL；实验界面分析可用，详见 easy_ppi/VALIDATION.zh-CN.md。

## 验证与审核

1. 在 easy-md-test 的 feature/fix 分支修改，记录基线 SHA 与影响范围。
2. 分别运行 Vina 单元/集成测试、独立 PPI 测试及本次功能的真实计算案例。新计算必须遵循既有授权流程并写入新结果目录；仓库整理无需重跑科学任务。
3. 保存可复现命令、环境版本、种子/盒参数、评分、DockQ/RMSD、运行时间、脱敏错误与修复。结构、图、PSE、全部候选和原日志保留本地，仓库只记录本地 run ID 和摘要。
4. 分别审核 **程序完成**、**文件/几何验证**、**科学有效性**。单元测试或退出码 0 不能替代科学验证。没有实验参考的 DockQ 必须 SKIPPED/null。
5. 标注待发布的明确提交和审核结论。算法案例失败/未知时，停止该算法晋升；修复后重新验证，不能自动把探索代码同步到正式仓库。

## 可追溯同步

两个仓库有各自历史，不直接跨仓库整树覆盖或 merge 所有实验提交。选择经审核的提交，在正式仓库的新分支上使用 `git cherry-pick -x <已审核的开发SHA>`，保留来源；若历史不同造成冲突，逐文件审阅解决，不能用强制覆盖。

在**单独的正式仓库 checkout**执行（其 origin 为 easy-MD，dev 为 easy-md-test）：

```bash
git fetch origin
git fetch dev
git switch -c release/<版本> origin/main
git cherry-pick -x <已审核的开发SHA>
# 检查实际 diff，只含已批准的稳定代码；不要带入 development/、本机配置或结果。
git diff origin/main --stat
# 使用本机已配置的 Python 运行 guard 和两套测试。
git push origin HEAD
```

用正式仓库 PR 将该分支与 main 比较。PR 必须链接测试仓库提交和运行记录（私有链接仅授权人员可读），并在 `validation/releases/<版本>.md` 提供公开、脱敏的必要验证摘要。
更新 `VERSION`、`CHANGELOG.md`，记录影响、兼容性、测试通过/跳过项和科学边界。先检查所有新增提交，防止先提交结构再删除仍留在历史中。

人工审核代码和证据，必需 CI 通过后才合入正式 main。仓库治理、说明和提交检查可作为独立 PR 发布，保持科学源码哈希不变；不要求为文档更新重跑对接。
合入后再次核对 main 的 SHA 与文件，在确认可发布的提交创建**新**版本 tag 和 GitHub Release；已有 tag/版本不移动、不覆盖。当前已有历史版本在引入 VERSION 前没有统一标签，不补造历史测试或标签。

## CI、hooks 与权限

`Repository checks` 的三个检查名为 `repository-policy`、`vina-tests`、`ppi-tests`。
CI 使用测试依赖，不安装完整 Vina/PLIP/Arpeggio/LightDock 环境；真实 PyMOL 或实验 fixture 缺失会按原测试跳过，因此绿色 CI 只证明已执行的自动化测试通过。
本地完整测试和真实案例证据必须另行审阅；CI 不发布、同步代码或启动真实对接。

正式 main 建议要求 PR、三个必需检查、禁止 force push/删除，并对管理员生效；私人开发仓库能否启用同等保护取决于 GitHub 方案。
读取远端保护配置确认是否生效，不能仅依据请求无异常。若权限/方案不足，在管理记录中保留原因并由仓库所有者在 Settings → Branches/Rules 配置；不得声称已实现服务器端拦截。

2026-11-20 实际回读：Public easy-MD/main 已启用 PR、上述三个必需检查、strict、管理员约束、禁止强推和删除；未要求额外独立批准人数（当前为 0），人工审核证据仍须保留在 PR。
Private easy-md-test 的保护配置/回读返回 403：当前方案要求升级 GitHub Pro 才能保护私有分支。本次保持 Private，不改变方案、不降低正式仓库保护；测试仓库目前靠本地 hooks、CI 和开发约定，不能保证 GitHub 网页或其他客户端的直接 push 被服务器阻止。
需要所有者升级 GitHub Pro 后，在 Settings → Branches 对 main 配置同样 PR/必需检查规则并回读验证；无需向智能体提供或粘贴 Token。

## 历史审查与报告范围

2026-11-20 只读审查了重命名前两个完整 Git 镜像中所有现存 refs 的可达提交、分支和标签。
具体路径、分支清单和发现仅供所有者在私有测试仓库 `development/records/` 查看，不在公开发布材料披露私有历史细节。
本次不删除历史分支、不重写历史、不强制推送。若未来要清理，必须另行确认范围并保留备份；.gitignore 或删除当前文件不会清除历史对象。已删除/不可达的 GitHub 对象不在可审查范围内。
