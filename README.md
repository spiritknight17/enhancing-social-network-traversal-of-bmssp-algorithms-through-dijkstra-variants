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

## Reproduce
```bash
# preprocess (example)
g++ -O2 -std=c++17 preprocess.cpp -o preprocess
# benchmark
g++ -O2 -std=c++17 benchmark_variants.cpp -o benchmark_variants
```
Or open `BMSSP_Dijsktra_Simulation_UI.ipynb`.
