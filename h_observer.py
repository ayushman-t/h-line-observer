#!/usr/bin/python3
# Copyright (c) 2026 Ayushman Tripathi
# Released under MIT License
#
# H-Line Observer — automated hydrogen line data collection for Raspberry Pi
# Requires VIRGO (pip install virgo) and an Airspy or RTL-SDR
#
# Usage:
#   python3 h_observer.py              interactive menu
#   python3 h_observer.py --auto       headless loop (for systemd)
#   python3 h_observer.py --install    first-time setup
#   python3 h_observer.py --config     view/change settings

import os
import sys
import json
import time
import signal
import shutil
import argparse
import subprocess
from pathlib import Path
from datetime import datetime, date

DEFAULTS = {
    'device':        'airspy=0,bias=1',
    'rf_gain':       15,
    'if_gain':       15,
    'bb_gain':       15,
    'frequency':     1420405750,
    'bandwidth':     3000000,
    'channels':      1024,
    'n_filter':      20,
    'm_filter':      35,
    'duration':      300,
    'sleep_between': 5,
    'output_base':   str(Path.home() / 'hydrogen_obs'),
}

CONFIG_FILE = Path.home() / '.config' / 'h_observer.conf'

PRESETS = {
    '1': {'name': 'Airspy Mini',           'device': 'airspy=0,bias=1', 'rf_gain': 15, 'if_gain': 15, 'bb_gain': 15},
    '2': {'name': 'Airspy R2',             'device': 'airspy=0,bias=1', 'rf_gain': 15, 'if_gain': 15, 'bb_gain': 15},
    '3': {'name': 'RTL-SDR (bias tee on)', 'device': 'rtl=0,bias=1',   'rf_gain': 49, 'if_gain': 0,  'bb_gain': 0},
    '4': {'name': 'RTL-SDR (no bias tee)', 'device': 'rtl=0',          'rf_gain': 49, 'if_gain': 0,  'bb_gain': 0},
    '5': {'name': 'HackRF One',            'device': 'hackrf=0',        'rf_gain': 0,  'if_gain': 40, 'bb_gain': 62},
    '6': {'name': 'HackRF One (+amp)',     'device': 'hackrf=0',        'rf_gain': 14, 'if_gain': 40, 'bb_gain': 62},
}


def load_config():
    cfg = dict(DEFAULTS)
    if CONFIG_FILE.exists():
        try:
            with open(CONFIG_FILE) as f:
                cfg.update(json.load(f))
        except Exception as e:
            print(f"[WARN] Could not read config: {e} — using defaults")
    return cfg


def save_config(cfg):
    CONFIG_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(CONFIG_FILE, 'w') as f:
        json.dump(cfg, f, indent=2)


def print_config(cfg):
    print(f"\n  Device   : {cfg['device']}")
    print(f"  Gains    : RF={cfg['rf_gain']}  IF={cfg['if_gain']}  BB={cfg['bb_gain']}")
    print(f"  Duration : {cfg['duration']}s")
    print(f"  Output   : {cfg['output_base']}")
    print(f"  Virgo    : {find_virgo() or 'NOT FOUND'}")


def cal_index_path(cfg):
    return Path(cfg['output_base']) / '.calibrations.json'


def load_calibrations(cfg):
    path = cal_index_path(cfg)
    if path.exists():
        try:
            with open(path) as f:
                return json.load(f)
        except Exception:
            pass
    return []


def save_cal_entry(cfg, entry):
    cals = load_calibrations(cfg)
    cals.append(entry)
    path = cal_index_path(cfg)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as f:
        json.dump(cals, f, indent=2)


def find_virgo():
    v = shutil.which('virgo')
    if v:
        return v
    for c in [
        Path.home() / '.local' / 'bin' / 'virgo',
        Path('/usr/local/bin/virgo'),
        Path('/usr/bin/virgo'),
    ]:
        if c.exists():
            return str(c)
    return None


def _pip_install(package):
    for extra in [[], ['--break-system-packages']]:
        try:
            r = subprocess.run(
                [sys.executable, '-m', 'pip', 'install', package] + extra,
                capture_output=True, text=True
            )
            if r.returncode == 0:
                return True
        except Exception:
            pass
    return False


def build_virgo_cmd(cfg, output_file, plot_file=None, cal_file=None, duration=None):
    virgo = find_virgo()
    if not virgo:
        raise RuntimeError("virgo not found — run: python3 h_observer.py --install")
    dur = duration if duration is not None else cfg['duration']
    cmd = [
        virgo,
        '-da', cfg['device'],
        '-f',  str(cfg['frequency']),
        '-b',  str(cfg['bandwidth']),
        '-c',  str(cfg['channels']),
        '-t',  '1',
        '-d',  str(dur),
        '-rf', str(cfg['rf_gain']),
        '-if', str(cfg['if_gain']),
        '-bb', str(cfg['bb_gain']),
        '-o',  output_file,
    ]
    if plot_file:
        cmd += ['-p', plot_file, '-n', str(cfg['n_filter']), '-m', str(cfg['m_filter'])]
    if cal_file:
        cmd += ['-C', cal_file]
    return cmd


def record_calibration(cfg, duration=60):
    print("\n" + "="*55)
    print("  RECORD CALIBRATION")
    print("="*55)
    print_config(cfg)
    print("\n  Point dish at cold sky (away from galactic plane)")
    input("  Press Enter when ready...")

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    out_dir   = Path(cfg['output_base']) / f'cal_{timestamp}'
    out_dir.mkdir(parents=True, exist_ok=True)
    cal_file  = str(out_dir / 'calibration.dat')

    print(f"\n  Recording {duration}s...")
    try:
        result = subprocess.run(
            build_virgo_cmd(cfg, cal_file, duration=duration),
            capture_output=True, text=True
        )
        if result.returncode != 0:
            print(f"  x {result.stderr.strip()}")
            return None
    except RuntimeError as e:
        print(f"  x {e}")
        return None

    save_cal_entry(cfg, {'name': f'Cal {timestamp}', 'file': cal_file, 'timestamp': timestamp})
    print(f"  + Saved: {cal_file}")
    return cal_file


def select_calibration(cfg):
    cals = load_calibrations(cfg)
    print("\n  Calibration:")
    print("    0. Skip")
    print("    1. Record new")
    if cals:
        print("\n    Saved:")
        for i, cal in enumerate(cals):
            status = "ok" if Path(cal['file']).exists() else "missing"
            print(f"    {i+2}. {cal['name']}  [{status}]")

    choice = input("\n  Select [0]: ").strip() or "0"
    if choice == "0":
        return None
    if choice == "1":
        dur = input("  Duration [60]: ").strip()
        return record_calibration(cfg, duration=int(dur) if dur.isdigit() else 60)
    try:
        idx = int(choice) - 2
        if 0 <= idx < len(cals):
            f = cals[idx]['file']
            if Path(f).exists():
                return f
            print("  file missing")
    except ValueError:
        pass
    return None


def run_observation(cfg, output_dir, obs_id, cal_file=None):
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    obs_file  = str(Path(output_dir) / f'{obs_id}_{timestamp}_observation.dat')
    plot_file = str(Path(output_dir) / f'{obs_id}_{timestamp}_plot.png')

    print(f"  [{datetime.now().strftime('%H:%M:%S')}] {obs_id} ({cfg['duration']}s)...",
          end=' ', flush=True)
    try:
        result = subprocess.run(
            build_virgo_cmd(cfg, obs_file, plot_file, cal_file),
            capture_output=True, text=True
        )
        if result.returncode == 0:
            print("ok")
            return True
        print(f"failed  {result.stderr.strip()[:80]}")
        return False
    except Exception as e:
        print(f"error  {e}")
        return False


def run_auto(cfg):
    running    = True
    output_dir = None
    count      = 0

    def _stop(signum, frame):
        nonlocal running
        print(f"\n[{datetime.now().strftime('%H:%M:%S')}] Stopping after current observation...")
        running = False

    signal.signal(signal.SIGTERM, _stop)
    signal.signal(signal.SIGINT,  _stop)

    cal_file = None
    for c in reversed(load_calibrations(cfg)):
        if Path(c['file']).exists():
            cal_file = c['file']
            break

    print("="*50)
    print("  H-LINE AUTO OBSERVER")
    print("="*50)
    print(f"  Device   : {cfg['device']}")
    print(f"  Gains    : RF={cfg['rf_gain']} IF={cfg['if_gain']} BB={cfg['bb_gain']}")
    print(f"  Duration : {cfg['duration']}s  Sleep: {cfg['sleep_between']}s")
    print(f"  Output   : {cfg['output_base']}")
    print(f"  Cal      : {Path(cal_file).parent.name if cal_file else 'none'}")
    print("="*50 + "\n")

    current_date = None

    while running:
        today = date.today()

        if today != current_date:
            current_date = today
            base         = Path(cfg['output_base'])
            existing     = sorted(base.glob(f'loop_{today.strftime("%Y%m%d")}*'))

            if existing:
                output_dir = str(existing[-1])
                count      = len(list(Path(output_dir).glob('*_observation.dat')))
                print(f"\n=== Resuming: {Path(output_dir).name} ({count} obs) ===\n")
            else:
                ts         = datetime.now().strftime('%Y%m%d_%H%M%S')
                output_dir = str(base / f'loop_{ts}')
                Path(output_dir).mkdir(parents=True, exist_ok=True)
                with open(Path(output_dir) / 'loop_info.txt', 'w') as f:
                    f.write(f"Started  : {datetime.now().isoformat()}\n")
                    f.write(f"Device   : {cfg['device']}\n")
                    f.write(f"Gains    : RF={cfg['rf_gain']} IF={cfg['if_gain']} BB={cfg['bb_gain']}\n")
                    f.write(f"Duration : {cfg['duration']}s\n")
                    f.write(f"Cal      : {cal_file or 'none'}\n")
                count = 0
                print(f"\n=== New loop: {Path(output_dir).name} ===\n")

        count += 1
        run_observation(cfg, output_dir, f'obs_{count:04d}', cal_file=cal_file)

        for _ in range(cfg['sleep_between']):
            if not running:
                break
            time.sleep(1)

    if output_dir:
        print(f"\nStopped. {count} observations in {Path(output_dir).name}")


def run_config(cfg):
    while True:
        print("\n" + "="*55)
        print("  CONFIGURATION")
        print("="*55)
        print_config(cfg)
        print("\n  1. Change SDR / gains")
        print("  2. Change output folder")
        print("  3. Change observation duration")
        print("  0. Back")
        choice = input("\n  Select: ").strip()

        if choice == "0":
            break

        elif choice == "1":
            print("\n  Select your SDR:")
            for k, p in PRESETS.items():
                print(f"    {k}. {p['name']}")
            print("    0. Manual")
            hw = input("\n  Select: ").strip()
            if hw in PRESETS:
                p = PRESETS[hw]
                cfg['device']  = p['device']
                cfg['rf_gain'] = p['rf_gain']
                cfg['if_gain'] = p['if_gain']
                cfg['bb_gain'] = p['bb_gain']
            elif hw == "0":
                v = input(f"  Device string [{cfg['device']}]: ").strip()
                if v: cfg['device'] = v
                v = input(f"  RF gain [{cfg['rf_gain']}]: ").strip()
                if v.isdigit(): cfg['rf_gain'] = int(v)
                v = input(f"  IF gain [{cfg['if_gain']}]: ").strip()
                if v.isdigit(): cfg['if_gain'] = int(v)
                v = input(f"  BB gain [{cfg['bb_gain']}]: ").strip()
                if v.isdigit(): cfg['bb_gain'] = int(v)
            save_config(cfg)
            print("  Saved.")

        elif choice == "2":
            v = input(f"  Output folder [{cfg['output_base']}]: ").strip()
            if v:
                cfg['output_base'] = str(Path(v).expanduser())
                Path(cfg['output_base']).mkdir(parents=True, exist_ok=True)
                save_config(cfg)
                print("  Saved.")

        elif choice == "3":
            v = input(f"  Duration [{cfg['duration']}s]: ").strip()
            if v.isdigit():
                cfg['duration'] = int(v)
                save_config(cfg)
                print("  Saved.")


def run_menu(cfg):
    print("\n" + "="*55)
    print("  H-LINE OBSERVER")
    print("="*55)
    print_config(cfg)

    while True:
        print("\n  1. Single observation")
        print("  2. Continuous loop")
        print("  3. Record calibration")
        print("  4. Config")
        print("  0. Exit")
        choice = input("\n  Select: ").strip()

        if choice == "0":
            break
        elif choice == "4":
            run_config(cfg)
        elif choice == "3":
            dur = input("\n  Duration [60]: ").strip()
            record_calibration(cfg, duration=int(dur) if dur.isdigit() else 60)
            input("\n  Press Enter...")
        elif choice in ("1", "2"):
            print_config(cfg)
            cal_file   = select_calibration(cfg)
            ts         = datetime.now().strftime('%Y%m%d_%H%M%S')
            label      = 'obs' if choice == '1' else 'loop'
            output_dir = str(Path(cfg['output_base']) / f'{label}_{ts}')
            Path(output_dir).mkdir(parents=True, exist_ok=True)

            if choice == "1":
                run_observation(cfg, output_dir, 'obs_0001', cal_file=cal_file)
                input("\n  Press Enter...")
            else:
                count = 0
                print("\n  Running — Ctrl+C to stop\n")
                try:
                    while True:
                        count += 1
                        run_observation(cfg, output_dir, f'obs_{count:04d}', cal_file=cal_file)
                        time.sleep(cfg['sleep_between'])
                except KeyboardInterrupt:
                    print(f"\n  Stopped after {count} observations")
                    input("\n  Press Enter...")


def run_install():
    print("\n" + "="*55)
    print("  H-LINE OBSERVER — SETUP")
    print("="*55)

    if sys.version_info < (3, 7):
        print(f"  Python 3.7+ required (you have {sys.version})")
        sys.exit(1)
    print(f"  Python {sys.version_info.major}.{sys.version_info.minor} ok")

    virgo = find_virgo()
    if virgo:
        print(f"  virgo: {virgo}")
    else:
        print("  virgo not found — installing...")
        if _pip_install('virgo'):
            virgo = find_virgo()
            if virgo:
                print(f"  virgo installed: {virgo}")
            else:
                print("  virgo installed but not in PATH")
                print("  Add to PATH: export PATH=$PATH:~/.local/bin")
        else:
            print("  pip install failed")
            print("  Try: pip install virgo --break-system-packages")

    script_path  = Path(os.path.abspath(__file__))
    service_dir  = Path.home() / '.config' / 'systemd' / 'user'
    service_dir.mkdir(parents=True, exist_ok=True)
    service_file = service_dir / 'h_observer.service'

    service_file.write_text(
        "[Unit]\n"
        "Description=H-Line Auto Observer\n"
        "After=network.target\n\n"
        "[Service]\n"
        f"ExecStart={sys.executable} {script_path} --auto\n"
        "Restart=always\n"
        "RestartSec=10\n\n"
        "[Install]\n"
        "WantedBy=default.target\n"
    )
    print(f"  Service file written")

    try:
        subprocess.run(['systemctl', '--user', 'daemon-reload'], check=True)
        subprocess.run(['systemctl', '--user', 'enable', 'h_observer'], check=True)
        print("  Service enabled (not started yet)")
    except Exception as e:
        print(f"  systemctl error: {e}")
        print("  Run manually:")
        print("    systemctl --user daemon-reload")
        print("    systemctl --user enable h_observer")

    print("\n" + "="*55)
    print("  Setup complete!")
    print("="*55)
    print(f"\n  Next steps:")
    print(f"    python3 {script_path.name} --config   set your SDR and output folder")
    print(f"    python3 {script_path.name}             record calibration, then start observing")
    print(f"\n  When ready to start the auto service:")
    print(f"    systemctl --user start h_observer")
    print(f"    systemctl --user status h_observer")
    print(f"    journalctl --user -u h_observer -f")


def check_deps():
    issues = []
    if not find_virgo():
        issues.append(('virgo not found', 'pip install virgo  (or with --break-system-packages)'))
    if not shutil.which('gnuradio-companion'):
        issues.append(('gnuradio missing', 'sudo apt install gnuradio'))
    if issues:
        print("\n  Dependency issues:")
        for name, cmd in issues:
            print(f"    {name:20s}  ->  {cmd}")
        print()


def main():
    parser = argparse.ArgumentParser(
        description='H-Line Observer — automated hydrogen line data collection'
    )
    parser.add_argument('--auto',    action='store_true', help='Headless mode for systemd')
    parser.add_argument('--install', action='store_true', help='First-time setup')
    parser.add_argument('--config',  action='store_true', help='View/change settings')
    args = parser.parse_args()
    check_deps()

    if args.install:
        run_install()
    elif args.config:
        cfg = load_config()
        run_config(cfg)
    elif args.auto:
        run_auto(load_config())
    else:
        run_menu(load_config())


if __name__ == '__main__':
    main()
