# H-Line Pipeline — Setup & Command Reference

> The `--process` step currently reprocesses every `.dat` file through VIRGO to produce calibrated CSVs. In a future update, the observer scripts will output this directly on the Pi, so only the plotting step will be needed.

## Setup

### 1. Install VIRGO and dependencies

```bash
pip install astro-virgo matplotlib==3.7.5 scipy plotly
```

matplotlib 3.7.5 is required — newer versions cause a `GrouperView` error in VIRGO.

plotly is optional (needed only for `--format html` interactive 3D plots).

### 2. Apply the VIRGO patch

The pipeline expects a 5th column in the processed CSV (the median-filtered calibrated spectrum). Stock VIRGO only outputs 4 columns. Apply the patch:

```bash
# Find where virgo is installed
python3 -c "import virgo; print(virgo.__file__)"

# Replace virgo.py with the patched version
cp virgo_patch.py /path/to/virgo/virgo.py
```

Verify it worked — after processing, CSV files should have 5 columns:
```bash
head -1 ~/hydrogen_obs_pipeline/processed/loop_*/obs_*_spectra_filtered.csv
# Should show 5 comma-separated values per line
```

### 3. Calibration file

The pipeline **will not work without a calibration file**. You need to record one before processing:

```bash
# Record calibration (point dish at cold sky, away from galactic plane)
python3 h_observer.py    # option 3
# or
python3 h_quick.py       # option 3
```

This creates a file like `~/hydrogen_obs/cal_20260208_001639/calibration.dat`. You can either:

- Copy it to the pipeline's default location:
  ```bash
  mkdir -p ~/hydrogen_obs_pipeline/calibration
  cp ~/hydrogen_obs/cal_*/calibration.dat ~/hydrogen_obs_pipeline/calibration/
  ```
- Or pass it every time with `--cal`:
  ```bash
  python3 hline_pipeline.py --process loop_20260215 --cal ~/hydrogen_obs/cal_20260208_001639/calibration.dat
  ```

### 4. Match your SDR settings

Open `hline_pipeline.py` and edit the `obs_params` dict in `reprocess_loop()` (around line 109) to match the SDR and gains you used during observation:

```python
obs_params = {
    'dev_args': 'airspy=0,bias=1',   # change to match your SDR (e.g. 'rtl=0,bias=1')
    'rf_gain': 15,                    # match your h_observer.py / h_quick.py gains
    'if_gain': 15,
    'bb_gain': 15,
    'frequency': 1420405750.0,
    'bandwidth': 3e6,                 # match your SDR bandwidth
    'channels': 1024,
    't_sample': 1,
    'duration': 300,
    ...
}
```

These must match what you used when collecting data. If you used h_observer.py with an RTL-SDR at RF gain 49, set the same values here. Check your h_observer config with `python3 h_observer.py --config`.

### 5. Configure your observatory

Pass your location and dish pointing via command line:

```bash
python3 hline_pipeline.py --lat 40.7 --lon -74.0 --altitude 10 --az 180 --el 45 ...
```

Or edit the defaults at the top of `hline_pipeline.py`.

---

## Workflow

```
h_observer.py collects .dat files
        |
        v
hline_pipeline.py --process    reprocesses with VIRGO + calibration -> CSV
        |
        v
hline_pipeline.py --plot3d     generates plots from processed CSVs
```

1. Collect data with `h_observer.py` (produces `loop_*/obs_*_observation.dat`)
2. Record a calibration pointing at cold sky (produces `cal_*/calibration.dat`)
3. Process: `--process LOOP --cal /path/to/calibration.dat`
4. Plot: `--plot3d`, `--heatmap`, `--flipbook`, etc.

---

## Command Reference

### List loops

```bash
./hline_pipeline.py --list
./hline_pipeline.py --data-dir /path/to/obs --list
```

### Process

```bash
./hline_pipeline.py --process loop_20260215 --cal ~/hydrogen_obs/cal_20260208/calibration.dat
./hline_pipeline.py --process latest --cal /path/to/calibration.dat
```

### 3D plots (single loop)

```bash
./hline_pipeline.py --plot3d loop_20260215                          # matplotlib 3D
./hline_pipeline.py --plot3d loop_20260215 --format html            # plotly interactive
./hline_pipeline.py --plot3d loop_20260215 --format png --save out.png
./hline_pipeline.py --plot3d loop_20260215 --y-axis velocity        # velocity axis
./hline_pipeline.py --plot3d loop_20260215 --y-axis velocity --vlsr # VLSR corrected
./hline_pipeline.py --plot3d loop_20260215 --x-axis time            # time instead of RA
./hline_pipeline.py --plot3d loop_20260215 --view top               # top-down heatmap
./hline_pipeline.py --plot3d loop_20260215 --render solid           # solid surface
./hline_pipeline.py --plot3d loop_20260215 --clean                  # no grid/panes
./hline_pipeline.py --plot3d loop_20260215 --quality high           # full detail
```

### 3D plots (multi-loop combined)

```bash
./hline_pipeline.py --plot3d loop_20260215 loop_20260216 loop_20260217
./hline_pipeline.py --plot3d loop_20260215 loop_20260216 --format html
./hline_pipeline.py --plot3d loop_20260215 loop_20260216 --view top
./hline_pipeline.py --plot3d loop_20260215 loop_20260216 --y-axis velocity --vlsr
```

### Heatmap

```bash
./hline_pipeline.py --heatmap loop_20260215
./hline_pipeline.py --heatmap loop_20260215 --y-axis velocity
./hline_pipeline.py --heatmap loop_20260215 --y-axis velocity --vlsr
./hline_pipeline.py --heatmap loop_20260215 --x-axis time
./hline_pipeline.py --heatmap loop_20260215 --save heatmap.png
```

### Flipbook animation

```bash
./hline_pipeline.py --flipbook loop_20260215                        # single loop GIF
./hline_pipeline.py --flipbook loop_20260215 loop_20260216          # multi-loop GIF
./hline_pipeline.py --flipbook loop_20260215 --y-axis velocity
./hline_pipeline.py --flipbook loop_20260215 --save flipbook.gif
```

### Single observation

```bash
./hline_pipeline.py --hours loop_20260215                           # list available hours
./hline_pipeline.py --single 5                                      # plot hour 5 (latest loop)
./hline_pipeline.py --single 12.5 --loop loop_20260215             # specific loop
```

### Compare (overlay loops)

```bash
./hline_pipeline.py --compare loop_20260215 loop_20260216           # matplotlib overlay
./hline_pipeline.py --compare loop_20260215 loop_20260216 --format html
./hline_pipeline.py --compare loop_20260215 loop_20260216 --z-gap 5 # manual spacing
./hline_pipeline.py --compare loop_20260215 loop_20260216 --y-axis velocity --vlsr
```

### Stack (co-add by RA)

Stacking multiple loops by RA bin improves SNR. Best for fixed-pointing dishes where the same sky passes through each sidereal day.

```bash
./hline_pipeline.py --stack loop_20260215 loop_20260216 loop_20260217
./hline_pipeline.py --stack loop_20260215 loop_20260216 --format html
./hline_pipeline.py --stack loop_20260215 loop_20260216 --ra-bin 0.05  # finer bins
./hline_pipeline.py --stack loop_20260215 loop_20260216 --vlsr         # VLSR corrected
./hline_pipeline.py --stack loop_20260215 loop_20260216 --view top     # 2D view
```

### Export (Zenodo-ready zip)

```bash
./hline_pipeline.py --export loop_20260215 loop_20260216 loop_20260217
./hline_pipeline.py --export loop_20260215 loop_20260216 --name my_dataset
./hline_pipeline.py --export loop_20260215 --dest ~/Desktop
```

### Full pipeline (process + plot everything)

```bash
./hline_pipeline.py --all loop_20260215 --cal /path/to/calibration.dat
./hline_pipeline.py --all loop_20260215 loop_20260216 --cal /path/to/calibration.dat --format html
```

---

## Common options

| Option | Description |
|---|---|
| `--data-dir DIR` | Where your loop folders are (default: `~/hydrogen_obs`) |
| `--cal FILE` | Path to calibration.dat |
| `--lat`, `--lon`, `--altitude` | Observatory location |
| `--az`, `--el` | Dish azimuth and elevation |
| `--x-axis ra\|time` | X-axis: right ascension or time from start |
| `--y-axis freq\|velocity` | Y-axis: frequency (MHz) or Doppler velocity (km/s) |
| `--vlsr` | Apply VLSR correction (requires astropy) |
| `--format show\|png\|html` | Output format |
| `--view 3d\|top` | 3D surface or top-down heatmap |
| `--render grid\|solid` | Surface style |
| `--quality fast\|medium\|high` | Matplotlib detail level |
| `--clean` | Remove grid lines and background panes |
| `--no-baseline` | Skip baseline correction |
| `--save PATH` | Save to file instead of showing |

---

Copyright (c) 2026 Ayushman Tripathi · MIT License
