---
name: easy-MD
description: Easy Molecular Docking with AutoDock Vina. Protein preparation (pdbfixer), ligand preparation (RDKit+OpenBabel), Vina docking, result packaging. Use whenever user mentions 分子对接, docking, virtual screening, protein-ligand binding, 蛋白对接, 虚拟筛选, 药靶对接, 结合模式, 结合能, 结合姿态 — even in passing or as part of a larger task. Triggers on: 对接, Vina, dock, protein-ligand, target prediction, binding mode, easy-MD.
---

# Molecular Docking — AutoDock Vina Automated Pipeline

Automate the complete molecular docking workflow: from target selection to publication-ready visualization. Uses tools already installed on the user's Windows system — no additional dependencies needed.

## 🛠 Setup — Verify Known Paths (MANDATORY, ~1 sec)

**Tools are installed at fixed paths on this machine. Do NOT scan the filesystem — just verify the known paths exist.** Filesystem scans are slow and can trigger Windows Store popups when they stumble on `python3.exe` stubs.

### Step 0: Quick verification of known paths (parallel, ~1 sec)

Run all 4 checks in parallel (single message with multiple Bash calls):

```bash
"/c/Python312/python.exe" --version 2>&1
```

```bash
"/d/Interaction/vina/vina_1.2.5_win.exe" --version 2>&1
```

```bash
"/d/Interaction/OpenBabel-3.1.1/obabel.exe" -V 2>&1
```

```bash
ls "/d/Interaction/PyMol/PyMOLWin.exe" 2>&1
```

### Known Paths (this machine — ALWAYS use these)

| Tool | Path | Verify |
|------|------|--------|
| Python | `C:\Python312\python.exe` | `--version` |
| Vina | `D:\Interaction\vina\vina_1.2.5_win.exe` | `--version` |
| OpenBabel | `D:\Interaction\OpenBabel-3.1.1\obabel.exe` | `-V` |
| PyMOL | `D:\Interaction\PyMol\PyMOLWin.exe` | `ls` (check exists) |

### Step 1: Report findings as a table

```
Tool          Status        Path
──────────────────────────────────────────────
Python        ✅ 3.12.0     C:\Python312\python.exe
Vina          ✅ 1.2.5      D:\Interaction\vina\vina_1.2.5_win.exe
OpenBabel     ✅ 3.1.1      D:\Interaction\OpenBabel-3.1.1\obabel.exe
PyMOL         ✅ 已安装      D:\Interaction\PyMol\PyMOLWin.exe
```

### Step 2: If anything is missing, show install guide

**Python 缺失时：**
```
请安装 Python 3.12 (64-bit)：
  下载: https://www.python.org/ftp/python/3.12.0/python-3.12.0-amd64.exe
  安装时勾选 "Add Python to PATH"，或装到 C:\Python312\
  验证: python --version
```

**Vina 缺失时：**
```
请下载 AutoDock Vina 1.2.5 (单文件，无需安装)：
  下载: https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_win.exe
  放到任意目录，如 C:\Users\<用户名>\bin\
  验证: vina_1.2.5_win.exe --version
```

**OpenBabel 缺失时：**
```
请安装 OpenBabel 3.1.1：
  下载: https://github.com/openbabel/openbabel/releases/download/openbabel-3-1-1/OpenBabel-3.1.1-Windows-64bit.exe
  安装到 D:\Interaction\OpenBabel-3.1.1\
  验证: obabel.exe -V
```

### Step 3: Check Python packages (ONLY after Python is confirmed)

```bash
"$PYTHON" -c "
import sys
pkg_map = {
    'rdkit': 'rdkit',
    'pdbfixer': 'pdbfixer', 
    'simtk.openmm': 'openmm',
    'Bio': 'biopython',
    'numpy': 'numpy',
    'MDAnalysis': 'MDAnalysis',
    'matplotlib': 'matplotlib',
}
missing = []
for mod, pkg in pkg_map.items():
    try:
        __import__(mod)
        print(f'  ✅ {pkg}')
    except ImportError:
        print(f'  ❌ {pkg} — NOT INSTALLED')
        missing.append(pkg)
if missing:
    print(f'\n  Run: \"{sys.executable}\" -m pip install {\" \".join(missing)}')
    print('  ⚠️  openmm is ~200MB, may take several minutes to download')
"
```

### Step 4: Store found paths — use them for the entire session

```python
# These are the ONLY Python/Vina/OpenBabel paths you may use:
PYTHON_EXE = "<whatever was found>"    # NOT "python3" or "python"
VINA_EXE   = "<whatever was found>"
OBABEL_EXE = "<whatever was found>"
```

**⛔ NEVER use bare `python3` or `python` in Bash!** They resolve to Windows Store stubs (`/c/Users/*/AppData/Local/Microsoft/WindowsApps/`). When Claude Code's sandbox detects a Store process launching, it kills it → **exit code 49**. Each kill triggers a retry with backoff, which is **why the session feels slow** — what looks like one failed command was actually 3-5 retries × 2-second backoff = 10+ seconds wasted per attempt.

**⛔ NEVER pass `/tmp/...` paths to Windows executables!** Git Bash `/tmp/` = Unix path. Windows binaries (`python.exe`, `obabel.exe`, `vina.exe`) receive it as a literal string and can't resolve it → **exit code 2**. Write all scripts to `C:\Users\22241\.claude\tmp_docking\` and pass full Windows-native paths.

**Note on SSL:** PubChem may return SSL errors on Windows (CRYPT_E_REVOCATION_OFFLINE). Always use `curl -k` (--insecure) when querying PubChem REST API.

**Note on Chinese paths:** Windows Chinese pathnames (e.g., `桌面/测试`) can cause Python multiprocessing garbled output and OpenBabel encoding issues. When the user specifies a path containing Chinese characters, transparently remap to an ASCII-only temporary working directory:
```
C:\Users\22241\.claude\tmp_docking\
```
Process all steps there, then copy final results back to the user's target directory. This avoids encoding bugs in pdbfixer, OpenBabel, and Vina without the user needing to know.

**CRITICAL — Path handling rules (must follow to avoid silent data loss):**
1. **All Python scripts MUST start with `os.chdir(r'C:\Users\22241\.claude\tmp_docking')`** — NEVER rely on Bash `cd` to set the working directory for Python. Bash `cd` in a multi-command pipeline may not propagate to the heredoc subprocess, causing files to scatter into `~/.claude/`.
2. **Desktop path is ALWAYS `D:\Desktop`** — your Windows Desktop maps to `D:\Desktop`, NOT `D:\桌面`. The latter is a regular folder, not the system Desktop. Use `r'D:\Desktop\...'` for Windows native paths and `/d/Desktop/...` for Git Bash.
3. **Never use Bash `cp` for Chinese paths** — Bash's path encoding under Git Bash is unreliable with multi-byte characters.
4. **All `cd` commands in Bash should be on the same line as the action**, e.g., `cd /path && some_command`. Do not assume the `cd` persists across separate tool calls.
5. **Python print() MUST avoid non-ASCII characters** — Windows Chinese systems default to GBK encoding. `Å` → `Angstrom`, `μM` → `uM`. Otherwise: `UnicodeEncodeError: 'gbk' codec can't encode character`.

**CRITICAL — User Preferences (ALL docking runs):**
1. **Exhaustiveness = 32** (高精度, publication quality). Never use 4 or 8.
2. **Repeat = 3** independent runs with different random seeds, merge results, take best.
3. **No interaction analysis** — skip prolif, py3Dmol, 2D interaction plots. Just docking scores + poses.
4. **Multi-model PDB** — merge all 9 poses from the best run into a single `all_poses.pdb` with MODEL/ENDMDL records. PyMOL can switch states.
5. **Output subfolder** — always create `{PROTEIN}_{LIGAND}/` under the user-specified target directory (e.g., `ROCK2_Belumosudil/`).
6. **summary.txt** — always generate a plain-text summary with protein info, ligand info, docking parameters, and score table.
7. **Clean up temp files** — after copying results to target, delete the `tmp_docking/` working directory.

---

## Workflow 0: Confirmation Step (MANDATORY — must run before any docking)

**NEVER skip this step.** Always present the search results to the user and wait for explicit confirmation before proceeding.

### 🔌 MCP-First Policy (CRITICAL)

**ALWAYS use MCP `bio-server` tools as the primary method** for PDB and PubChem searches. The MCP tools handle SSL, encoding, error handling, and formatting automatically. Only fall back to curl/REST API when:

- MCP tool returns an error or no results
- You need to post-filter results (e.g., filter by X-RAY method, check resolution)
- `get_pdb_info` is called (⚠️ known bug: returns `'str' object has no attribute 'get'` — use REST API instead)

**MCP tools reference:**

| MCP Tool | Replaces | ⚠️ Known Issues & Fallback |
|----------|----------|---------------------------|
| `search_pdb(keyword, max_results)` | RCSB Search API curl | **偶发字段全 N/A/?**：并行 REST 未追到时，需手动 `curl + core/entry` 逐条补全 metadata |
| `search_ligand(name)` | PubChem REST API curl | **SMILES 可能为 `?`**（较新药物）：自动 fallback 到 2D record；若仍失败，手动 `curl + record_type=2d` |
| `get_smiles(name)` | PubChem canonical SMILES curl | **较新药物返回 "Not found"**：三级回退（name→CID→2D record），见下方 Ligand Search 章节 |
| `download_pdb(pdb_id, output_dir)` | biopython PDBList | Direct PDB file download |
| `get_pdb_info(pdb_id)` | ⚠️ **BROKEN** — 永远不要调用 | 返回 `'str' object has no attribute 'get'`。用 REST API：`curl -k -s "https://data.rcsb.org/rest/v1/core/entry/{PDB_ID}"` |

### Protein Search

**If user provides a PDB ID directly** (e.g., "4RWE"):
- Validate by fetching metadata from `https://data.rcsb.org/rest/v1/core/entry/{PDB_ID}`
- Show: PDB ID, title, resolution, method, organism, deposited date, ligand list

**If user provides a protein name/keyword** (e.g., "EGFR", "COX-2", "HIV-1 protease"):
- **Step 0 (MANDATORY): Resolve protein name → gene symbol** — use MyGene.info API if searching human proteins. Skip for viral/bacterial targets (use the name directly).
  - Human targets: `curl -k -s "https://mygene.info/v3/query?q={protein_name}&species=human&fields=symbol,name" | /c/Python312/python.exe -c "import sys,json; d=json.load(sys.stdin); hits=d.get('hits',[]); [print(h['symbol'], h.get('name','')) for h in hits[:5]]"`
  - Use the **gene symbol** (not common name) for human targets — common abbreviations are ambiguous (e.g. "COX-2" matched JMJD2D via full-text, but "PTGS2" found the actual cyclooxygenase-2)
  - **Viral/bacterial targets** (e.g., HIV-1 protease, SARS-CoV-2 Mpro): use the protein name directly — no gene symbol resolution needed
- **Step A: Search PDB via MCP `search_pdb()`** (PRIMARY — 1 call, ~2 sec):
  - Call `search_pdb("{gene_symbol or protein_name}", max_results=10)`
  - This returns PDB IDs sorted by resolution — much faster and more reliable than raw RCSB API
  - **Note:** `search_pdb` already filters for human proteins by default. For non-human targets (viral, bacterial), the search still works — just use the target name directly (e.g., `search_pdb("HIV-1 protease")`)
  - **⚠️ N/A Fallback:** 如果 `search_pdb` 结果中 Resolution / Method / Year / Ligands 全部显示 N/A 或 `?`（并行 REST 偶发失败），**立即用以下 one-liner 批量补全所有条目的 metadata**：
    ```bash
    /c/Python312/python.exe -c "
    import json,urllib.request,ssl
    ctx=ssl.create_default_context(); ctx.check_hostname=False; ctx.verify_mode=ssl.CERT_NONE
    for pid in ['8X92','7JNT','9EP8','8X8X','7JOV','8GDS','7P6O']:
        try:
            req=urllib.request.Request(f'https://data.rcsb.org/rest/v1/core/entry/{pid}')
            with urllib.request.urlopen(req,context=ctx,timeout=10) as r:
                d=json.loads(r.read())
            m=d.get('exptl',[{}])[0].get('method','?')
            r=d.get('rcsb_entry_info',{}).get('resolution_combined',['?'])
            t=d.get('struct',{}).get('title','?')[:60]
            ligs=[e.get('rcsb_chem_comp_info',{}).get('chem_comp_id','?') for e in d.get('rcsb_entry_container_identifiers',{}).get('non_polymer_entity_ids',[])[:3]]
            print(f'{pid} | {m:25s} | {r}A | Ligands: {ligs} | {t}')
        except Exception as e:
            print(f'{pid} | ERROR: {e}')
    "
    ```
    不要逐条等 `search_pdb` 重试 — 批量 REST 一次拿到所有 metadata。
- **Step B: Post-filter by experimental method** — MCP `search_pdb` may return NMR/Cryo-EM alongside X-RAY. Filter to keep only X-RAY DIFFRACTION:
  - For each candidate PDB ID, fetch `https://data.rcsb.org/rest/v1/core/entry/{PDB_ID}` (NOT `get_pdb_info` — it's broken)
  - Check `exptl[0].method` for "X-RAY DIFFRACTION"
  - Extract resolution, title, organism, ligands
  - Keep only X-RAY entries, sorted by resolution
- **⚠️ `get_pdb_info` is BROKEN** — returns `'str' object has no attribute 'get'`. Use REST API instead: `curl -k -s "https://data.rcsb.org/rest/v1/core/entry/{PDB_ID}"`
- **Fallback (only if MCP `search_pdb` fails or returns 0 results):** Use RCSB Search API directly:
  ```json
  {
    "query": {
      "type": "group",
      "logical_operator": "and",
      "nodes": [
        {"type": "terminal", "service": "full_text", "parameters": {"value": "{KEYWORD}"}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entity_source_organism.scientific_name", "operator": "exact_match", "value": "Homo sapiens"}}
      ]
    },
    "return_type": "entry",
    "request_options": {
      "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
      "paginate": {"start": 0, "rows": 10}
    }
  }
  ```
  - **CRITICAL:** Do NOT filter by `rcsb_entry_info.structure_determination_methodology` — this field always stores `"experimental"`, not `"X-RAY DIFFRACTION"`.
  - **CRITICAL:** Use `rcsb_entity_source_organism.scientific_name` (NOT `src_organism_ids.name` — that attribute does not exist).
- **CRITICAL — Empty header fallback:** Some PDB entries (e.g., 7JNT, 8X92) return empty header files at `files.rcsb.org/header/`. NEVER silently drop an entry for this — instead fetch `https://data.rcsb.org/rest/v1/core/entry/{PDB_ID}` and check `exptl[0].method`.
- Handle the case where no results match: relax organism filter, try broader keywords
- Present results as a table

### Ligand Search

**If user provides SMILES directly** → use as-is, skip PubChem search
**If user provides a molecule name** → search PubChem:
- **PRIMARY: Use MCP `search_ligand(name)`** — returns CID, formula, MW, SMILES, IUPAC in one call
- **For SMILES only:** Use MCP `get_smiles(name)` — returns Canonical SMILES string directly
- Parse JSON response, present formatted table

#### ⚠️ SMILES Fallback Chain (CRITICAL — 必须执行直到拿到有效 SMILES)

较新药物（如 Belumosudil/KD025, 2021年批准）的 `CanonicalSMILES` 在 PubChem property 端点常返回 `?` 或空值。**MCP 工具内置了三级回退**，但如果 MCP 仍然失败，手动执行以下 fallback：

| Level | Method | Command | 说明 |
|:-----:|--------|---------|------|
| 1 | MCP `search_ligand` / `get_smiles` | (primary, 已内置 fallback) | 大部分常见药物在这一级即可 |
| 2 | PubChem property by CID | `curl -k -s "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{CID}/property/CanonicalSMILES/JSON"` | 如果 Level 1 拿到 CID 但 SMILES 为空 |
| 3 | PubChem **2D record** | `curl -k -s "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{NAME}/record/JSON?record_type=2d"` | **最可靠**：从 `PC_Compounds[0].props` 中搜 `label=="SMILES"` |

**Level 3 提取 SMILES 的 one-liner：**
```bash
curl -k -s "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{NAME}/record/JSON?record_type=2d" | /c/Python312/python.exe -c "
import sys,json; d=json.load(sys.stdin)
for p in d['PC_Compounds'][0]['props']:
    if p['urn']['label'] in ('SMILES','Canonical SMILES'):
        print(p['value']['sval']); break
"
```

**备选名尝试**：如果原名查不到，依次尝试别名（如 Belumosudil → KD025 → SLx-2119）。PubChem 对商业名和研发代号的索引覆盖不同。

#### PubChem REST fallback (MCP fails entirely)

```bash
curl -k -s "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{NAME}/property/Title,MolecularFormula,MolecularWeight,CanonicalSMILES,IUPACName/JSON"
```

### Confirmation Prompt

Always end Workflow 0 with:
```
✅ Please confirm:
   Protein: {PDB_ID} — {Title} ({Resolution}Å)
   Ligand:  {Name} (CID: {CID})
   Parameters: exhaustiveness=32, repeat=3
   Proceed with docking? [yes / modify / cancel]
```

Do NOT continue to Workflow 1 until user explicitly confirms.

### ⚡ Speed Tip: Parallel Search

PDB and PubChem searches are independent — **always run them simultaneously** (two MCP calls in one step). MCP `search_pdb()` + `search_ligand()` run in parallel for ~2-second total wait time.

---

## ⏱️ Expected Runtime

**Default parameters: exhaustiveness=32, repeat=3**

| Stage | Typical Time | Bottleneck |
|-------|-------------|------------|
| Workflow 0 (search) | 3–8 sec | Network → parallel RCSB+PubChem |
| Workflow 1 (protein prep) | 10–60 sec | pdbfixer → scales with protein size |
| Workflow 2 (ligand prep) | 2–5 sec | RDKit → trivial |
| Workflow 3 (box definition) | 1–3 sec | Trivial |
| Workflow 4 (Vina docking ×3) | 20–40 min | CPU → exhaustiveness=32 ×3 repeats |
| Workflow 5 (result packaging) | 5–10 sec | File copy + summary |
| **Total** | **~25–45 min** | Vina ×3 dominates |

### Exhaustiveness & Repeat Guide

**Default: exhaustiveness=32, repeat=3** (严谨出图 / publication quality)

| exhaustiveness | ×1 Run | ×3 Runs | Use Case |
|:---:|--------|---------|----------|
| 4 | ~1 min | ~3 min | Quick test only |
| 8 | ~4 min | ~12 min | Draft / sanity check |
| 16 | ~10 min | ~30 min | Standard (if user explicitly asks) |
| **32** | **~25 min** | **~40 min** | **Publication quality (DEFAULT)** |

**Always use 32 unless the user explicitly asks for faster.** Never default to 8.

---

## Workflow 1: Protein Preparation

1. **Download PDB file** using biopython:
```python
import os
os.chdir(r'C:\Users\22241\.claude\tmp_docking')  # MUST be first — do NOT rely on Bash cd

from Bio.PDB import PDBList
pdbl = PDBList()
pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')
```

2. **Detect and select chain(s)** — MANDATORY. Most docking targets only need ONE chain. Multi-chain structures (dimers, complexes) must be reduced:
```python
from Bio.PDB import PDBParser, PDBIO, Select

class ChainSelect(Select):
    def __init__(self, chain_ids):
        self.chain_ids = set(chain_ids)
    def accept_chain(self, chain):
        return chain.id in self.chain_ids

# Parse and report all chains
parser = PDBParser(QUIET=True)
struct = parser.get_structure(pdb_id, pdb_file)
chains = [(c.id, len(list(c.get_residues()))) for c in struct[0].get_chains()]
print(f"Found {len(chains)} chain(s):")
for cid, nres in chains:
    print(f"  Chain {cid}: {nres} residues")

# Selection logic:
# - If 1 chain → use as-is
# - If multiple → default to chain A, inform user
# - If user specified a chain → use that one
selected_chain = 'A'  # default

# Extract selected chain
io = PDBIO()
io.set_structure(struct)
single_chain_pdb = os.path.join(output_dir, f'{pdb_id}_chain{selected_chain}.pdb')
io.save(single_chain_pdb, ChainSelect({selected_chain}))
```

3. **Prepare protein** with pdbfixer (on the single-chain PDB):
```python
import os
os.chdir(r'C:\Users\22241\.claude\tmp_docking')  # MUST be first

from pdbfixer import PDBFixer
from openmm.app import PDBFile

fixer = PDBFixer(filename=single_chain_pdb)
fixer.removeHeterogens(keepWater=False)
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.4)
PDBFile.writeFile(fixer.topology, fixer.positions, open(f'{output_dir}/protein_clean.pdb', 'w'))
```

4. **Convert to PDBQT** with OpenBabel:
```bash
"D:\Interaction\OpenBabel-3.1.1\obabel.exe" protein_clean.pdb -O receptor.pdbqt -xr --partialcharge gasteiger
```

---

## Workflow 2: Ligand Preparation

1. **SMILES → 3D** with RDKit:
```python
import os
os.chdir(r'C:\Users\22241\.claude\tmp_docking')  # MUST be first

from rdkit import Chem
from rdkit.Chem import AllChem

mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())

# MMFFOptimizeMolecule API varies by RDKit version:
#   Old API: returns int (0=converged, 1=not converged)
#   New API: returns tuple (status, energy)
result = AllChem.MMFFOptimizeMolecule(mol)
if isinstance(result, tuple):
    print(f"  MMFF: converged={result[0]==0}, energy={result[1]:.2f}")
else:
    print(f"  MMFF: status={result} (0=converged)")

writer = Chem.SDWriter('ligand.sdf')
writer.write(mol)
writer.close()
```

2. **Convert to PDBQT** with OpenBabel:
```bash
"D:\Interaction\OpenBabel-3.1.1\obabel.exe" ligand.sdf -O ligand.pdbqt --gen3d --partialcharge gasteiger
```

---

## Workflow 3: Docking Box Definition

Auto-detect the binding site from a reference ligand in the original PDB (if co-crystal exists). Otherwise ask user for center coordinates.

**Auto-detection code:**
```python
import numpy as np

def detect_binding_site(pdb_file):
    het_residues = {}
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('HETATM'):
                resn = line[17:20].strip()
                if resn in ('HOH', 'WAT', 'H2O', 'GOL', 'EDO', 'PEG', 'SO4', 'PO4'):
                    continue
                chain = line[21:22]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                key = f'{resn}_{chain}'
                if key not in het_residues:
                    het_residues[key] = []
                het_residues[key].append((x, y, z))
    
    if not het_residues:
        return None, None
    
    # Prefer chain A ligand, largest by atom count
    chain_a = {k: v for k, v in het_residues.items() if k.endswith('_A')}
    primary_key = max(chain_a, key=lambda k: len(chain_a[k])) if chain_a else max(het_residues, key=lambda k: len(het_residues[k]))
    coords = np.array(het_residues[primary_key])
    center = coords.mean(axis=0)
    
    min_xyz = coords.min(axis=0) - 8.0
    max_xyz = coords.max(axis=0) + 8.0
    box_size = max((max_xyz - min_xyz).max(), 22.0)
    
    return center, box_size

center, box_size = detect_binding_site(pdb_file)
```

**Write config.txt:**
```python
import os
os.chdir(r'C:\Users\22241\.claude\tmp_docking')  # MUST be first

config = f"""receptor = receptor.pdbqt
ligand = ligand.pdbqt

center_x = {center[0]:.2f}
center_y = {center[1]:.2f}
center_z = {center[2]:.2f}

size_x = {box:.0f}
size_y = {box:.0f}
size_z = {box:.0f}

num_modes = 9
energy_range = 3
exhaustiveness = 32
"""
with open('config.txt', 'w') as f:
    f.write(config)
```

---

## Workflow 4: Run Vina (Repeat 3× with different seeds)

**Always run Vina 3 times with different random seeds** to ensure the global minimum is found. Use seeds: 42, 2024, 8888. Run all 3 in parallel for speed.

**⚠️ CRITICAL: Use a two-step approach.** Git Bash can race `mkdir` return against `cd` in a `&` pipeline, causing `No such file or directory` for run2/run3. Always verify directories exist before launching Vina.

**Step A (sync):** Create directories and copy input files, then verify:
```bash
cd /c/Users/22241/.claude/tmp_docking && mkdir -p run1 run2 run3 && \
  cp receptor.pdbqt ligand.pdbqt config.txt run1/ && \
  cp receptor.pdbqt ligand.pdbqt config.txt run2/ && \
  cp receptor.pdbqt ligand.pdbqt config.txt run3/ && \
  ls run1/config.txt run2/config.txt run3/config.txt && echo "All directories ready"
```
The trailing `ls` confirms files are actually on disk before proceeding.

**Step B (async):** Launch all 3 runs simultaneously in background, then wait. **CRITICAL: Add `cd /` before `wait`** — this releases the shell's CWD from the run directories so Windows can delete them later:
```bash
cd /c/Users/22241/.claude/tmp_docking/run1 && "D:\Interaction\vina\vina_1.2.5_win.exe" --config config.txt --out docked.pdbqt --seed 42 2>&1 | tee log.txt &
cd /c/Users/22241/.claude/tmp_docking/run2 && "D:\Interaction\vina\vina_1.2.5_win.exe" --config config.txt --out docked.pdbqt --seed 2024 2>&1 | tee log.txt &
cd /c/Users/22241/.claude/tmp_docking/run3 && "D:\Interaction\vina\vina_1.2.5_win.exe" --config config.txt --out docked.pdbqt --seed 8888 2>&1 | tee log.txt &
cd / && wait && echo "=== All Vina runs complete ==="
```
The `cd /` before `wait` is essential: it moves the shell out of `run3/`, releasing the directory handle on Windows. Without it, `tmp_docking/` cannot be deleted.

**Note:** Use separate `run1/`, `run2/`, `run3/` directories (not files like `run1.pdbqt`) — this prevents Vina output file conflicts in parallel execution and keeps logs organized.

**After all 3 runs complete, compare best scores and select the best run:**

```python
import os, re
os.chdir(r'C:\Users\22241\.claude\tmp_docking')  # MUST be first

best_affinity = 999
best_run = None
best_seed = None
all_scores = []
seeds = {1: 42, 2: 2024, 3: 8888}

for run_num, seed in seeds.items():
    log_file = f'run{run_num}/log.txt'
    with open(log_file) as f:
        log_text = f.read()
    
    # Extract first (best) affinity from scoring table
    m = re.search(r'\s+1\s+([\-\d.]+)', log_text)
    if m:
        aff = float(m.group(1))
        all_scores.append((run_num, seed, aff))
        if aff < best_affinity:
            best_affinity = aff
            best_run = run_num
            best_seed = seed

print(f"Best run: #{best_run} (seed={best_seed}, affinity={best_affinity} kcal/mol)")
for r, s, a in all_scores:
    mark = " <-- BEST" if r == best_run else ""
    print(f"  Run {r} (seed={s}): {a:.3f} kcal/mol{mark}")

# Split best run's docked.pdbqt into individual PDB poses
# Use OpenBabel (see Workflow 5)
```

**Present the comparison table:**
```
⚡ Vina ×3 Complete:
┌────────┬──────┬──────────────────┐
│ Run    │ Seed │ Best Affinity    │
├────────┼──────┼──────────────────┤
│ #1     │ 42   │ -9.8 kcal/mol    │
│ #2     │ 2024 │ -10.2 kcal/mol ⭐│
│ #3     │ 8888 │ -9.6 kcal/mol    │
└────────┴──────┴──────────────────┘
Selected: Run #2 (seed=2024)
```

Also parse the full 9-mode score table from the best run's log and present it.

---

## Workflow 5: Result Packaging (NO interaction analysis)

### 5a. Split docked PDBQT and merge all 9 poses into single multi-model PDB

First split the best run's `docked.pdbqt` into individual PDB files using OpenBabel:
```bash
cd /c/Users/22241/.claude/tmp_docking && "D:\Interaction\OpenBabel-3.1.1\obabel.exe" run{BEST_RUN}/docked.pdbqt -O best_pose.pdb -m
```
This produces `best_pose1.pdb` through `best_pose9.pdb`.

Then merge them into a multi-model PDB with REMARK scores:

```python
import os, re
os.chdir(r'C:\Users\22241\.claude\tmp_docking')  # MUST be first

# Parse scores from the best run's log
scores = []
best_run = 2  # determined in Workflow 4
with open(f'run{best_run}/log.txt') as f:
    log_text = f.read()
for line in log_text.split('\n'):
    m = re.match(r'\s+(\d+)\s+([\-\d.]+)\s+([\-\d.]+)\s+([\-\d.]+)', line)
    if m:
        scores.append((int(m.group(1)), float(m.group(2)), float(m.group(3)), float(m.group(4))))

# Merge all 9 PDB poses into multi-model PDB
with open('all_poses.pdb', 'w') as outf:
    for i in range(1, 10):
        pose_file = f'best_pose{i}.pdb'
        if not os.path.exists(pose_file):
            continue
        aff = scores[i-1][1]
        rmsd_lb = scores[i-1][2]
        rmsd_ub = scores[i-1][3]
        outf.write(f'MODEL     {i}\n')
        outf.write(f'REMARK Vina Score: {aff:.3f} kcal/mol  RMSD l.b.: {rmsd_lb:.3f}  RMSD u.b.: {rmsd_ub:.3f}\n')
        with open(pose_file, 'r') as inf:
            for line in inf:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    outf.write(line)
        outf.write('ENDMDL\n')
```

### 5b. Save scores.csv

```python
import os
os.chdir(r'C:\Users\22241\.claude\tmp_docking')  # MUST be first

with open('scores.csv', 'w') as f:
    f.write('mode,affinity_kcal_mol,rmsd_lb,rmsd_ub\n')
    for mode, aff, lb, ub in scores:
        f.write(f'{mode},{aff:.3f},{lb:.3f},{ub:.3f}\n')
```

### 5c. Generate summary.txt

```python
import os
os.chdir(r'C:\Users\22241\.claude\tmp_docking')  # MUST be first

summary = f"""============================================================
  分子对接报告 -- {PROTEIN} + {LIGAND_NAME}
============================================================

[蛋白信息]
  基因名: {GENE_NAME}
  PDB ID: {PDB_ID}
  标题: {TITLE}
  物种: Homo sapiens
  方法: X-RAY DIFFRACTION
  分辨率: {RESOLUTION} A
  链: {CHAIN} ({N_RES} residues)
  共晶配体: {COCRYSTAL_LIGAND}

[配体信息]
  名称: {LIGAND_NAME} ({LIGAND_NAME_CN})
  PubChem CID: {CID}
  分子式: {FORMULA}
  分子量: {MW} g/mol
  SMILES: {SMILES}

[对接参数]
  软件: AutoDock Vina 1.2.5
  Exhaustiveness: 32
  重复次数: 3 (seeds: 42, 2024, 8888)
  Box center: {CX}, {CY}, {CZ}
  Box size: {BOX} x {BOX} x {BOX} A
  Number of modes: 9
  Energy range: 3 kcal/mol

[三组独立对接结果]
  Run      Seed     Best Affinity
  ----------------------------------------
  Run 1    42       {run1_best:.3f}       kcal/mol
  Run 2    2024     {run2_best:.3f}       kcal/mol{' <-- BEST' if best_run == 2 else ''}
  Run 3    8888     {run3_best:.3f}       kcal/mol{' <-- BEST' if best_run == 3 else ''}

[最佳对接结果 -- Run {best_run} (seed={best_seed})]
  Mode     Affinity         RMSD l.b.    RMSD u.b.
  --------------------------------------------------
"""
for mode, aff, lb, ub in scores:
    star = " <--" if mode == 1 else ""
    summary += f"  {mode:<8} {aff:<16.3f} {lb:<12.3f} {ub:.3f}{star}\n"

summary += f"""
[输出文件]
  all_poses.pdb  -- 多模型PDB (9个对接构象, PyMOL state切换)
  scores.csv     -- 结合能分数表
  summary.txt    -- 本报告
  receptor.pdbqt -- 受体 (Vina格式)
  ligand.pdbqt   -- 配体 (Vina格式)
  protein_clean.pdb -- 处理后的蛋白结构
  best_pose.pdb  -- 最佳构象 (Mode 1)
  config.txt     -- 对接盒子配置

============================================================
"""

with open('summary.txt', 'w', encoding='utf-8') as f:
    f.write(summary)
```

### 5d. Copy to target and clean up

**MANDATORY ORDER — follow exactly or risk data loss:**

**Step 1:** Create subfolder and copy files using Windows native paths in Python:

```python
import os, shutil
os.chdir(r'C:\Users\22241\.claude\tmp_docking')  # MUST be first

# Use Windows native path format — NEVER /d/.../ Unix-style for Chinese paths
target = r'{TARGET_DIR}\{PROTEIN}_{LIGAND}'  # e.g., r'D:\Desktop\测试\ROCK2_Belumosudil'
os.makedirs(target, exist_ok=True)

files_to_copy = [
    ('all_poses.pdb', 'all_poses.pdb'),
    ('protein_clean.pdb', 'protein_clean.pdb'),
    ('scores.csv', 'scores.csv'),
    ('summary.txt', 'summary.txt'),
    ('config.txt', 'config.txt'),
    ('receptor.pdbqt', 'receptor.pdbqt'),
    ('ligand.pdbqt', 'ligand.pdbqt'),
    ('ligand.sdf', 'ligand.sdf'),
    ('best_pose1.pdb', 'best_pose.pdb'),      # rename to best_pose.pdb
    ('run{best_run}/log.txt', 'log.txt'),      # best run's log
    ('4QTB.pdb', '{pdb_id}.pdb'),              # original PDB
]

for src, dst in files_to_copy:
    if os.path.exists(src):
        shutil.copy2(src, os.path.join(target, dst))
```

**Step 2: VERIFY files before deleting temp (CRITICAL — never skip this):**

```python
import os

target = r'{TARGET_DIR}\{PROTEIN}_{LIGAND}'
required = ['all_poses.pdb', 'summary.txt', 'scores.csv', 'best_pose.pdb',
            'receptor.pdbqt', 'ligand.pdbqt', 'protein_clean.pdb']

all_ok = True
for f in required:
    path = os.path.join(target, f)
    exists = os.path.exists(path)
    size = os.path.getsize(path) if exists else 0
    status = f'OK ({size} bytes)' if exists and size > 0 else 'MISSING or EMPTY!'
    print(f'  {f}: {status}')
    if not exists or size == 0:
        all_ok = False

if all_ok:
    print('All files verified. Safe to delete temp.')
else:
    print('WARNING: Some files missing or empty. DO NOT delete temp!')
    # Re-run the copy step to fix
```

**Step 3: Delete temp directory ONLY after verification passes.**

⚠️ **Windows file locking:** Vina subprocesses and `tee` hold directory handles (CWD) for seconds after exit. Before cleanup, `cd /` to release the CWD, then use Python with retry:

```python
import os, shutil, time
os.chdir('/')  # Release CWD from tmp_docking
d = r'C:\Users\22241\.claude\tmp_docking'
for i in range(10):
    try:
        shutil.rmtree(d)
        print('Temp cleaned')
        break
    except OSError:
        time.sleep(2)
        if i == 9:
            print('Note: temp not cleaned (files locked). Results are safe in target — harmless.')
```

**Key rules:**
- Always use Windows native paths (`r'D:\Desktop\...'`), NEVER Unix-style (`/d/桌面/...`)
- The Desktop is `D:\Desktop` (system path), NOT `D:\桌面` (regular folder)
- NEVER delete `tmp_docking/` before verifying files exist at target with correct sizes
- If verification fails, re-run Step 1 — do NOT proceed to deletion
- Bash `cp` is unreliable with Chinese paths — use Python `shutil.copy2` exclusively

### Output Directory Structure

```
{target_dir}/{PROTEIN}_{LIGAND}/
├── all_poses.pdb         # ★ 9个姿态合并 (PyMOL用state切换)
├── protein_clean.pdb     # 准备好的蛋白结构
├── ligand.sdf            # 配体3D结构
├── ligand.pdbqt          # Vina配体文件
├── receptor.pdbqt        # Vina受体文件
├── config.txt            # 对接盒子参数
├── log.txt               # Vina完整日志(最佳run)
├── scores.csv            # 分数表
└── summary.txt           # ★ 完整对接报告
```

**NO interaction analysis files** — no `2d_interaction.png`, `3d_view.html`, `interactions.txt`, or `best_pose*.pdb` individual files. All poses are in `all_poses.pdb`.

---

## Troubleshooting

| Problem | Likely Cause | Fix |
|---------|-------------|-----|
| Bash exit code 49 when using `python3` | `python3` is a Windows Store stub (`/c/Users/22241/AppData/Local/Microsoft/WindowsApps/python3`) that gets killed by sandbox | **Always use `"C:/Python312/python.exe"`** — never `python3` or `python` in any Bash command, heredoc, or pipe. |
| Bash exit code 2 "file not found" with `/tmp/` paths | Git Bash `/tmp/` is not recognized by Windows executables (`C:\Python312\python.exe`, `obabel.exe`, etc.) | **Always use Windows-native paths** (`C:/Users/22241/...`). Write scripts to `C:\Users\22241\.claude\tmp_docking\`, never `/tmp/`. |
| `TypeError: 'int' object is not subscriptable` from `MMFFOptimizeMolecule` | RDKit version mismatch — old API returns `int`, new API returns `tuple` | Use `isinstance(result, tuple)` to check before subscripting (see Workflow 2 updated code). |
| Protein too large / unnatural results | Multi-chain PDB, dimer/tetramer | **Always extract single chain first.** Dimer interfaces inflate scores artificially. |
| pdbfixer crashes | PDB has exotic residues | Add `fixer.findNonstandardResidues()` + `fixer.replaceNonstandardResidues()` |
| Vina "Parse error" | PDBQT malformed | Re-run obabel with `--partialcharge gasteiger` flag |
| No docking poses | Box too small / wrong location | Increase size_x/y/z to 25 or manually specify center |
| PubChem SSL error | Windows schannel CRYPT_E_REVOCATION_OFFLINE | Use `curl -k` flag |
| RCSB search returns 0 results | Wrong attribute names in query | Use `rcsb_entity_source_organism.scientific_name` for organism. Filter X-RAY by fetching PDB headers (`EXPDTA` record). |
| Empty search results | Keyword too specific | Try alternative names: gene symbol, full protein name, or broader terms. |
| PDB header file returns empty (e.g., 7JNT, 8X92) | Some newer PDB entries have empty header files at `files.rcsb.org/header/` | Fall back to REST API `data.rcsb.org/rest/v1/core/entry/{PDB_ID}` → `exptl[0].method`. Do NOT silently drop the entry. See Workflow 0 Step B updated code. |
| Chinese path encoding error | Python multiprocessing on Windows | Process all files in `C:\Users\22241\.claude\tmp_docking\` (ASCII-only), copy results back. |
| Chinese path copy fails silently | Using `/d/桌面/...` Unix paths in Python | Use Windows native paths `r'D:\Desktop\...'` in Python `shutil.copy2`. Desktop = D:\Desktop, NOT D:\桌面. |
| Files scattered in `~/.claude/` | Bash `cd` not propagated to Python heredoc | Every Python script must start with `os.chdir(r'C:\Users\22241\.claude\tmp_docking')`. Never trust Bash `cd`. |
| Vina results vary with same `--seed` | Multi-threading race in Vina 1.2.5 | Expected behavior. `--seed` sets initial seed but thread scheduling is non-deterministic. Running 3× and taking best is the correct strategy. |
| Temp files deleted before copy verified | Wrong order in Workflow 5 | ALWAYS verify files at target (Step 2) BEFORE deleting temp (Step 3). If verification fails, re-copy, never skip to deletion. |
| Background bash find+loop times out | Multi-drive find + version check loop exceeds timeout | Split into separate short commands (one per tool). Each takes 2-3 seconds. See Step 0 updated commands. |
| `AttributeError: module 'rdkit.Chem' has no attribute 'GetFormattedSmiles'` | RDKit version — this function does not exist | Use `Chem.MolToSmiles(mol)` instead. |
| `UnicodeEncodeError: 'gbk' codec can't encode character '\xc5'` | Windows Chinese system defaults to GBK, non-ASCII characters (e.g. Å) in print() crash | Avoid all non-ASCII in Python print(): use `Angstrom` for Å, `uM` for μM. See Path Handling Rule #5. |
| Vina run2/run3: `No such file or directory` | Git Bash `mkdir` + `&` race — directory not flushed to disk before `cd` | Use two-step approach: Step A creates dirs+verifies files (sync), Step B launches Vina (async). See Workflow 4. |
| `rm -rf tmp_docking`: Device or resource busy | Windows file locking — Vina/bash subprocess CWD still points into `tmp_docking/` | (1) Add `cd /` before `wait` in Vina Step B to release shell CWD. (2) Python cleanup: `os.chdir('/')` before `shutil.rmtree()`. (3) Retry 10× with 2s delay. Failure is harmless. |
