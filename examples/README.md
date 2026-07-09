# examples/ — hackable Python workflows

These are **copy-and-edit starting points** for building custom TOPO simulation
workflows in Python. Unlike [`tutorials/`](../tutorials/) — which teach the
`.ini`-driven console commands (`topo-mdrun`, `topo-csp`) step by step — the
scripts here open up the *composition* those commands hide, so you can insert
your own logic (custom forces, temperature schedules, on-the-fly analysis)
between the building blocks.

Every script uses only TOPO's public API (`topo`, `topo.engine`); nothing here
reaches into private internals.

| Script | What it shows |
| --- | --- |
| [`custom_md.py`](custom_md.py) | The open version of `topo-mdrun`: `build_system → setup_simulation → attach_reporters → run_protocol → finalize_simulation`, with marked edit points for adding forces, defining a temperature schedule, and running in segments with a callback. |

## How to use

1. Copy a script somewhere of your own (don't edit it in place — that keeps the
   reference clean and update-friendly):

   ```bash
   cp examples/custom_md.py my_run.py
   ```

2. Edit the sections marked `# === EDIT: ... ===`.

3. Run it against a control file, exactly like the console command:

   ```bash
   python my_run.py -f md.ini
   ```

For the fully-featured runners (anneal/quench, restart, multi-copy plumbing),
read the canonical sources these examples are distilled from:
[`topo/mdrun/mdrun.py`](../topo/mdrun/mdrun.py) and
[`topo/engine.py`](../topo/engine.py).
