# enhancing-social-network-traversal-of-bmssp-algorithms-through-dijkstra-variants

BMSSP algorithms through Dijkstra variants evaluated on social-network traversals (LiveJournal, Orkut, soc-Pokec).

## Contents
- `bench_cli.cpp`, `benchmark_variants.cpp`, `preprocess.cpp` — C++ benchmark / preprocessing sources
- `bench_cli.exe` — prebuilt Windows CLI (rebuild from source if needed)
- `COLAB_REAL_RUN.py` — Colab reproduction script
- `BMSSP_Dijsktra_Simulation_UI.ipynb` — simulation UI notebook
- `ui_runs/` — per-run CSV outputs (`results_all_raw.csv`, `results_summary.csv`, `statistics_tests.csv`, per-variant logs) + slice `.npy` selectors (small enough to version)
- `*.csv` at root — cross-env summaries and statistics tests
- `Cloud Drive.txt` — Google Drive folder with full datasets

## Large datasets (not versioned)
Raw graphs and diffusion binaries exceed GitHub's 100 MB/file limit, so they are excluded via `.gitignore` and hosted on Drive:

- `livejournal.txt` (~1.0 GB), `livejournal_diffusion.bin` (~567 MB)
- `orkut.txt` (~1.7 GB), `orkut_diffusion.bin` (~949 MB)
- `soc_pokec.txt` (~423 MB), `soc_pokec_diffusion.bin` (~251 MB)
- `ui_runs/*slice*.bin` (91–528 MB each)
- `ABP - Video demonstration of code.mov` (~560 MB)

Download: https://drive.google.com/drive/folders/1yoDYA0uBHfwelBtRDgy9kd2slnP7XW36?usp=sharing

Place the files next to the executables/scripts to reproduce locally. See `COLAB_REAL_RUN.py` and the notebook for the expected layout.

## Run the Gradio web UI (`BMSSP_Dijsktra_Simulation_UI.ipynb`)

The notebook builds and launches a 4-tab Gradio app (`BMSSP vs Dijkstra Explorer`):

1. **How CSR works** — clickable CSR explainer + peek at a real graph.
2. **Run** — tick graphs + 14 algorithms (`Dijkstra`, `BMSSP_Org/SBBT/Bin/Fib_2025/2026`, `BMSSP_WC/CD/Expected`, `BMSSPy`, `RecDijkstra_AC`), optional statistical (Nx repeats, Mann-Whitney U + Holm) and line-graph (25/50/75% slices) runs. Each algorithm runs in its own process and is verified against Dijkstra.
3. **Watch it work** — narrated step-by-step animations on a small sample + comparison tab.
4. **Dashboard** — KPI cards, race, speed-up map, speed vs memory, memory over time, size line graph, statistics, details table.

Run these notebook cells **in order, once each**:

| Order | Cell | What it does |
|---|---|---|
| 1 | Step 0+1 Environment + Data Gathering | Detects Colab vs local, records CPU/RAM, sets `DATA_DIR` / `WORK`, downloads SNAP graphs to `DATA_DIR` if missing |
| 2 | `%%writefile preprocess.cpp` + compile/run guard | Writes preprocessor; skips if all 3 `*_diffusion.bin` exist |
| 3 | `%%writefile benchmark_variants.cpp` | Writes all algorithm implementations (must run before the UI build) |
| 4 | Install | `pip install -q 'gradio>=6,<7' plotly scipy` (plus `numpy/pandas` locally if missing) |
| 5 | Build CLI harness | Compiles `bench_cli` (`g++ -O3 -std=c++20`) and `csr_convert` (`g++ -O2 -std=c++17`) verbatim from `benchmark_variants.cpp`; needs Boost headers |
| 6 | UI core (helpers) | Defines CSR I/O, engine, checkpoint/resume, tracers, dashboard helpers |
| 7 | Gradio app (`build_app`) | Defines `build_app()` / `launch_app()` |
| 8 | Launch (`demo = build_app(); launch_app(demo)`) | Starts Gradio with `share=True` by default; open the `*.gradio.live` public link (full-screen) or local `http://127.0.0.1:7860` |

> Do not skip cells 1 and 3 — the UI reads `ENV`/`DATA_DIR` from cell 1 and the algorithm source from cell 3.

### A. Google Colab

1. Upload/open `BMSSP_Dijsktra_Simulation_UI.ipynb` in Colab (any runtime works; algorithms are CPU-only, T4/GPU not required).
2. Run cell 1 first. Approve the **Google Drive mount** prompt when it appears. This sets:
   - `DATA_DIR = /content/drive/MyDrive/GraphDatasets`
   - `WORK = /content/ui_runs`, checkpoints in `DATA_DIR/explorer_checkpoint` (survives browser reload; `Resume` survives runtime reset).
3. If the 3 `*_diffusion.bin` files are not in Drive yet, cell 1 downloads + cell 2 preprocesses the SNAP graphs (`soc_pokec`, `livejournal`, `orkut`). Or copy them from the Drive folder below into `GraphDatasets/` to skip this.
4. Run the remaining cells in table order. `ensure_boost()` auto-installs `libboost-all-dev` via `apt-get`; `g++` is preinstalled.
5. In the launch cell, wait for the `Running on public URL: https://xxxx.gradio.live` output and open it (allow popups). Keep the Colab tab open — closing the browser tab does not stop the run; reopening the link re-attaches via the 1.5 s poll timer.
6. If the Colab runtime disconnects: Runtime > Reconnect, re-run cells 1, 3–8, open the new link, press **Resume** on the Run tab (finished jobs load from `explorer_checkpoint`). **Start fresh** clears the saved session.
7. Large-graph note: Orkut + `BMSSP_CD`/`BMSSPy` need ~100 bytes/edge (~tens of GB). On low-RAM runtimes leave `Skip algorithms that need more RAM` ticked; skipped rows appear as `SKIPPED`, not failures.

### B. Visual Studio Code (local)

Best on **Linux / WSL2-Ubuntu / macOS**. Native Windows works for the UI but the per-algorithm sandbox (`resource` + `preexec_fn` memory cap, `/proc` sampling, `symlink`) is POSIX-only — use WSL2 for full benchmark behavior.

1. Prerequisites:
   - Python 3.10+ (notebook metadata: 3.13), VS Code + **Jupyter** extension, a `g++` toolchain (`g++ --version` must work).
   - Boost headers: `sudo apt-get install -y libboost-all-dev` (Ubuntu/WSL) or `brew install boost` (macOS).
   - VS Code selected interpreter / venv with:
     ```bash
     pip install jupyter ipykernel "gradio>=6,<7" plotly scipy pandas numpy
     python -m ipykernel install --user --name bmssp
     ```
2. Clone/open this repo folder in VS Code. Place datasets in `./GraphDatasets/` (default local `DATA_DIR`), or point elsewhere:
   ```bash
   # Windows PowerShell
   $env:GRAPH_DIR = "C:\path\to\GraphDatasets"
   # Linux/macOS/WSL
   export GRAPH_DIR=/path/to/GraphDatasets
   ```
   Local defaults when `GRAPH_DIR` is unset: `DATA_DIR = <repo>/GraphDatasets`, `WORK = <repo>/ui_runs`. Expected files: `soc_pokec_diffusion.bin`, `livejournal_diffusion.bin`, `orkut_diffusion.bin` (download from the Drive link below, or let cells 1–2 download + build them from SNAP).
3. Open `BMSSP_Dijsktra_Simulation_UI.ipynb` in VS Code, select the `bmssp` kernel, then run cells in the table order above (Run Below works; avoid running the launch cell twice — stop the first server first).
4. Build check: after cell 5 you must have `./bench_cli` (`./bench_cli.exe` on native Windows — prefer WSL so the UI's `./bench_cli` path resolves) and `./csr_convert`. If compile fails, the cell prints `stderr`; install Boost / fix `g++` and re-run.
5. Run the launch cell:
   ```python
   demo = build_app()
   launch_app(demo)              # share=True: local + *.gradio.live link
   launch_app(demo, share=False) # local-only: http://127.0.0.1:7860
   ```
   Open the local URL or public link in a browser (VS Code may offer Open in Simple Browser). Port `7860` must be free; forwarding is automatic in WSL2 (`localhost` maps to Windows).
6. Same resume model as Colab, but checkpoints live in `<DATA_DIR>/explorer_checkpoint` locally. Custom uploads (`.txt/.csv/.gz/.tar.gz/.bin`) are converted via `csr_convert` (Step 2 pipeline, diffusion weights `w=-log(p)*10000`, seed 42) and copied next to the checkpoint.

### Troubleshooting

- `benchmark_variants.cpp not found` → run cell 3 first.
- `bench_cli is not built yet` → re-run cell 5; check `g++` + Boost errors.
- `ModuleNotFoundError: gradio/scipy` → re-run cell 4 (VS Code: ensure the right kernel/venv).
- Windows `AttributeError: module 'os' has no attribute 'setsid'` / `import resource` failure → expected on native Windows; rerun from WSL2.
- No public link → set `share=True` and allow outbound HTTPS; or use `share=False` + local URL.
- Out-of-memory / timeout rows → normal on small machines; lower the graph/algorithm selection, keep RAM-skip on, raise `Timeout per algorithm` only on large-RAM machines.

## Reproduce (CLI only, no UI)
```bash
# preprocess (example)
g++ -O2 -std=c++17 preprocess.cpp -o preprocess
# benchmark
g++ -O2 -std=c++17 benchmark_variants.cpp -o benchmark_variants
```
