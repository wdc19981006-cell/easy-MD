"""Concise, explicit reports; never turn incomplete calculations into zeros."""
from pathlib import Path
import json


def load(path):
    return json.loads(path.read_text(encoding='utf-8')) if path.is_file() else None


def write_report(directory, state, output=None):
    directory = Path(directory)
    lines = ['# easy-PPI 运行报告', '', f'运行目录：`{directory}`', '',
             '## 阶段状态', '', '|阶段|状态|说明|', '|---|---|---|']
    for name, stage in state['stages'].items():
        explanation = stage.get('error') or stage.get('reason') or ''
        lines.append(f'|{name}|{stage["status"]}|{explanation.replace(chr(10), " ").replace("|", "/")}|')
    lines += ['', '## 科学验证', '', '```json', json.dumps(state['scientific_validation'], ensure_ascii=False, indent=2), '```', '']
    sasa = load(directory / 'analysis/sasa.json')
    if sasa:
        lines += ['## 实验或当前构象界面', '',
                  f'- SASA(A/B/AB)：{sasa["sasa"]["A"]:.2f} / {sasa["sasa"]["B"]:.2f} / {sasa["sasa"]["AB"]:.2f} Å²。',
                  f'- 总埋藏 SASA = A+B−AB：{sasa["buried_total_A2"]:.2f} Å²。',
                  f'- 按一半定义的界面面积：{sasa["interface_area_half_A2"]:.2f} Å²。', '']
    contacts = load(directory / 'analysis/arpeggio/interactions.json')
    if contacts:
        lines += ['## PDBe Arpeggio', '', f'跨链原子接触记录：{len(contacts["atom_contacts"])}；多数 proximal 仅表示邻近，不等于特异化学作用。', '',
                  '标签可重叠，不能将计数相加当作独立相互作用数：', '', '```json', json.dumps(contacts['counts_by_label'], ensure_ascii=False, indent=2), '```', '',
                  'ionic 为离子接触/盐桥候选，未宣称实验验证的盐桥。', '']
    protonation = load(directory / 'analysis/arpeggio/protonation.json')
    if protonation:
        lines += [f'加氢：{protonation["added_hydrogens"]} 个 H；重原子最大位移 {protonation["heavy_atom_max_shift_A"]:.4g} Å。', '']
        if protonation.get('warnings'):
            lines += ['OpenBabel 原始警告须保留审阅：', '', '```text', '\n'.join(protonation['warnings']), '```', '']
    affinity = load(directory / 'analysis/affinity_empirical.json')
    if affinity:
        lines += ['## PRODIGY 经验预测', '', f'ΔG = {affinity["delta_g_kcal_mol"]} kcal/mol；Kd = {affinity["kd_M"]:.3g} M，温度 {affinity["temperature_C"]} °C。', '',
                  '这是结构接触模型的经验预测，不能代替实验亲和力。', '']
    else:
        lines += ['PRODIGY：未得到可解释的预测；无数值替代。', '']
    redocking = load(directory / 'validation/redocking.json')
    if redocking:
        top, best = redocking['score_top1'], redocking['oracle_best_for_sampling_only']
        lines += ['## LightDock 重新对接', '', f'共保存并评价 {redocking["candidate_count"]} 个最终候选。', '',
                  '|指标|LightDock 评分 top1|参考评价最佳（仅用于采样覆盖）|', '|---|---|---|',
                  *[f'|{key}|{top[key]:.4f}|{best[key]:.4f}|' for key in ('DockQ', 'fnat', 'iRMSD', 'LRMSD')],
                  f'|LightDock 排名|{top["rank"]}|{best["rank"]}|',
                  f'|<1.5 Å 严重重叠原子计数（A侧）|{top["severe_clash_atom_count_A_1_5A"]}|{best["severe_clash_atom_count_A_1_5A"]}|', '',
                  f'采样验证：{redocking["sampling_validation"]}；top10 排名验证：{redocking["ranking_validation"]}。', '',
                  'DockQ <0.23 为不正确；0.23–0.49 可接受；0.49–0.80 中等；≥0.80 高质量。DockQ 最佳构象不得冒充评分首位。', '',
                  '本测试为 bound redocking，不代表非结合态结构对接或普适性能；严重重叠需要进一步物理检查/精修。', '']
    lines += ['## 产物与边界', '', '原始结构与获取证据在 raw/；编号映射在 prepared/；全部评分、候选和聚类在 docking/；原生 stdout/stderr 在 logs/。', '',
              'PyMOL 会话和 PNG 在 visualization/；重对接首位构象的分析保存在 docking/top1_analysis/。manifest.json 用于校验文件完整性，文件完整性通过不代表科学验证通过。', '',
              '缺失残基/侧链未自动补齐；实验水分子没有进入蛋白直接接触模型；质子化是推断。详细局限见检查 JSON 和架构文档。', '']
    destination = Path(output) if output else directory / 'REPORT.zh-CN.md'
    destination.write_text('\n'.join(lines), encoding='utf-8')
