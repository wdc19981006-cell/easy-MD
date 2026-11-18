# easy-MD repository guidance

Read `README.md`, `NEW_USER_GUIDE.md`, and `SKILL.md` before configuring or running this workflow. The repository root is the skill root; resolve scripts relative to it.

For a new user, check the operating system and actual Python environment first. This workflow targets Windows. The executable and output paths in `config.yaml` are examples from the original workstation: configure real local paths before invoking the pipeline. Do not assume the original author's programs or bio-server MCP are installed. Use bundled REST acquisition adapters when the MCP tools are unavailable.

A request to download alone authorizes downloading and inspection. Install dependencies when the user has authorized installation; the README first-use prompt provides that authorization for an isolated workflow environment. Preserve existing configurations and result directories. Check `scripts/preflight.py`, actual PLIP/OpenBabel imports, and the configured PyMOL Python before promising the full analysis and visualization workflow.

Follow the prepare / confirmation / execute sequence in `SKILL.md`. Dependency installation consent is not approval for Vina or redocking. Present the resolved inputs, site evidence, box, seeds, exhaustiveness, output directory, and expected runtime, then wait for approval of that concrete run. Report ambiguity or missing site evidence instead of silently choosing blind docking.

Classify interactions only with PLIP. Preserve the four-object PyMOL scene and its geometry validation. Keep output data, virtual environments, caches and credentials out of source commits. Report actual stage failures and validation results. Do not promise a fixed six-minute runtime or claim a fresh-machine installation test that has not been performed.
