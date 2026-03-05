# H-Line Observer

> Automated hydrogen line data collection for Raspberry Pi · Requires [VIRGO](https://virgo.readthedocs.io/en/latest/examples.html) + GNU Radio

---

## Prerequisites

- **Python 3.7+** — `sudo apt install python3`
- **GNU Radio** — required by VIRGO for SDR data acquisition
- **[VIRGO](https://virgo.readthedocs.io/en/latest/examples.html)** — the spectrometer engine both scripts use under the hood

```bash
sudo apt install python3 gnuradio
pip install virgo --break-system-packages
```

---

## Which Script?

| | `h_quick.py` | `h_observer.py` |
|---|---|---|
| **Purpose** | Manual testing tool | Full automated pipeline |
| **Use case** | Run it once, pick an option, observe. Good for checking your setup works before going automated. | Runs 24/7 as a background service, creates daily folders, auto-resumes after power cuts. |

---

## `h_quick.py` — Manual

1. **Check gains at the top of the file**
   Defaults to **Airspy Mini / R2** — if that's your SDR you don't need to change anything. Otherwise edit `DEVICE`, `RF_GAIN`, `IF_GAIN`, `BB_GAIN` to match your hardware.

2. **Run it**
   ```bash
   python3 h_quick.py
   ```

3. **Record a calibration first**
   Choose option 3, point dish at cold sky, press Enter.

4. **Take an observation**
   Option 1 for single, option 2 for continuous loop. Select your saved calibration when asked.

> **Note:** Data saves to `~/hydrogen_obs/` — edit `OUTPUT_BASE` at the top to change this.

---

## `h_observer.py` — Automated

1. **Install and set up the service**
   ```bash
   python3 h_observer.py --install
   ```

2. **Set your SDR and output folder**
   ```bash
   python3 h_observer.py --config
   ```
   Pick your hardware preset (Airspy Mini, RTL-SDR, HackRF, etc.) and where to save data.

3. **Record a calibration**
   ```bash
   python3 h_observer.py
   ```
   Choose option 3, point dish at cold sky, press Enter.

4. **Start the service**
   ```bash
   systemctl --user start h_observer
   systemctl --user status h_observer
   ```
   The service starts automatically on every login from now on.

> **Note:** To watch live logs: `journalctl --user -u h_observer -f` · To stop: `systemctl --user stop h_observer`

---

## Output

Each loop creates a folder: `~/hydrogen_obs/loop_YYYYMMDD_HHMMSS/`

Inside: one `.dat` file and one `.png` plot per observation, plus a `loop_info.txt` with settings.

Calibrations are saved separately under `~/hydrogen_obs/cal_*/` and indexed in `.calibrations.json`.

---

Copyright (c) 2026 Ayushman Tripathi · MIT License
