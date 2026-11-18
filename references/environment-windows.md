# Windows environment and path rules

## Contents

- Tool configuration
- Skill-root resolution
- Work-directory rules
- PowerShell invocation
- Git Bash invocation
- Compatibility notes

## Tool configuration

Read executable paths from `config.yaml`. The distributed defaults describe
the current workstation:

- Python: `C:/Python312/python.exe`
- Vina: `D:/Interaction/vina/vina_1.2.5_win.exe`
- OpenBabel: `D:/Interaction/OpenBabel-3.1.1/obabel.exe`
- PyMOL CLI: `D:/Interaction/PyMol/Scripts/pymol.exe`
- PyMOL Python: `D:/Interaction/PyMol/python.exe`

Run `scripts/preflight.py` before the workflow. Do not search the whole
filesystem when a configured path is missing; report the missing path and ask
whether the user wants to update `config.yaml`.

## Skill-root resolution

Resolve scripts relative to the loaded `SKILL.md`:

```text
<skill-root>/scripts/preflight.py
<skill-root>/scripts/protein_prepare.py
<skill-root>/scripts/pipeline.py
```

Do not embed `.claude`, `.codex`, `.agents`, or a user name in script paths.

Claude Code can expose `${CLAUDE_SKILL_DIR}`. Codex supplies the skill location
when it loads the skill. These are host details; bundled Python scripts must
remain host-neutral.

## Work-directory rules

- Create one dedicated directory per run.
- Use an ASCII-safe path for Windows scientific tools.
- Pass its absolute path through `--work-dir`.
- Never pass Git Bash `/tmp/...` to Windows executables.
- Never use a home directory, drive root, skill directory, or project root as
  the work directory.
- Existing non-empty directories without `.easy-md-run` are rejected.

## PowerShell invocation

Use PowerShell's call operator:

```powershell
& "C:\Python312\python.exe" "<skill-root>\scripts\preflight.py" `
  --config "<skill-root>\config.yaml"
```

Pass every argument as a separate quoted token. Do not translate commands into
background pipelines; `docking_run.py` manages the repeated processes.

## Git Bash invocation

Invoke Python with its Git Bash executable path, but pass Windows-native paths
to Python and Windows programs:

```bash
"/c/Python312/python.exe" "C:/path/to/easy-md/scripts/preflight.py" \
  --config "C:/path/to/easy-md/config.yaml"
```

Do not rely on an `export` from an earlier tool call. Pass `--work-dir`
explicitly to write-capable compatibility scripts. The `pipeline.py`
orchestrator instead accepts `--output-root`/`--run-dir` and passes the managed
work directory to its child scripts.

## Compatibility notes

The shared `SKILL.md` uses only standard `name` and `description` frontmatter.
Codex-specific UI metadata lives in `agents/openai.yaml`. Avoid adding
Claude-only frontmatter fields to the shared file if strict Codex validation is
required.
