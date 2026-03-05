#!/usr/bin/python3
# Copyright (c) 2026 Ayushman Tripathi
# Released under MIT License
#
# H-Line Quick — manual observation tool for testing
# For automated 24/7 operation use h_observer.py instead

import os
import json
import subprocess
from datetime import datetime

DEVICE    = 'airspy=0,bias=1'
RF_GAIN   = 15
IF_GAIN   = 15
BB_GAIN   = 15
FREQUENCY = 1420405750
BANDWIDTH = 3000000
CHANNELS  = 1024
N_FILTER  = 20
M_FILTER  = 35

OUTPUT_BASE = os.path.expanduser('~/hydrogen_obs')
CAL_INDEX   = os.path.join(OUTPUT_BASE, '.calibrations.json')


def find_virgo():
    import shutil
    v = shutil.which('virgo')
    if v:
        return v
    for p in [
        os.path.expanduser('~/.local/bin/virgo'),
        '/usr/local/bin/virgo',
        '/usr/bin/virgo',
    ]:
        if os.path.exists(p):
            return p
    return None


def load_calibrations():
    if os.path.exists(CAL_INDEX):
        with open(CAL_INDEX) as f:
            return json.load(f)
    return []


def save_calibration(cal_info):
    cals = load_calibrations()
    cals.append(cal_info)
    os.makedirs(os.path.dirname(CAL_INDEX), exist_ok=True)
    with open(CAL_INDEX, 'w') as f:
        json.dump(cals, f, indent=2)


def select_calibration():
    cals = load_calibrations()
    print("\n  Calibration:")
    print("    0. Skip")
    print("    1. Record NEW")
    if cals:
        print("\n    Saved:")
        for i, cal in enumerate(cals):
            status = "ok" if os.path.exists(cal['file']) else "missing"
            print(f"    {i+2}. {cal['name']}  [{status}]")
    choice = input("\n  Select [0]: ").strip() or "0"
    if choice == "0": return None
    if choice == "1": return "NEW"
    try:
        idx = int(choice) - 2
        if 0 <= idx < len(cals) and os.path.exists(cals[idx]['file']):
            return cals[idx]['file']
    except ValueError:
        pass
    return None


def record_calibration():
    print("\n" + "="*55)
    print("  RECORD CALIBRATION")
    print("="*55)
    dur = input("\n  Duration [60]: ").strip()
    duration = int(dur) if dur.isdigit() else 60
    print("\n  Point dish at cold sky")
    input("  Press Enter when ready...")

    timestamp  = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir = os.path.join(OUTPUT_BASE, f'cal_{timestamp}')
    os.makedirs(output_dir, exist_ok=True)
    cal_file   = os.path.join(output_dir, 'calibration.dat')

    print(f"\n  Recording {duration}s...")
    virgo = find_virgo()
    if not virgo:
        print("  virgo not found")
        return None

    subprocess.run([
        virgo,
        '-da', DEVICE, '-f', str(FREQUENCY), '-b', str(BANDWIDTH),
        '-c', str(CHANNELS), '-t', '1', '-d', str(duration),
        '-rf', str(RF_GAIN), '-if', str(IF_GAIN), '-bb', str(BB_GAIN),
        '-o', cal_file
    ])

    save_calibration({'name': f'Cal {timestamp}', 'file': cal_file, 'timestamp': timestamp})
    print(f"\n  Saved: {cal_file}")
    return cal_file


def run_observation(duration, cal_file, output_dir, obs_id):
    print(f"\n{'='*55}\n  {obs_id}\n{'='*55}")
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    obs_file  = os.path.join(output_dir, f'{obs_id}_{timestamp}_observation.dat')
    plot_file = os.path.join(output_dir, f'{obs_id}_{timestamp}_plot.png')

    virgo = find_virgo()
    if not virgo:
        print("  virgo not found")
        return

    print(f"\n  Recording {duration}s...")
    cmd = [
        virgo,
        '-da', DEVICE, '-f', str(FREQUENCY), '-b', str(BANDWIDTH),
        '-c', str(CHANNELS), '-t', '1', '-d', str(duration),
        '-rf', str(RF_GAIN), '-if', str(IF_GAIN), '-bb', str(BB_GAIN),
        '-o', obs_file, '-p', plot_file,
        '-n', str(N_FILTER), '-m', str(M_FILTER)
    ]
    if cal_file:
        cmd += ['-C', cal_file]
    subprocess.run(cmd)
    print(f"\n  Done: {os.path.basename(plot_file)}")


def main():
    print("\n" + "="*55 + "\n  H-LINE QUICK\n" + "="*55)
    os.makedirs(OUTPUT_BASE, exist_ok=True)

    virgo = find_virgo()
    print(f"  Virgo  : {virgo or 'NOT FOUND — pip install virgo'}")
    print(f"  Output : {OUTPUT_BASE}")

    while True:
        print("\n  1. Single observation")
        print("  2. Continuous loop")
        print("  3. Record calibration")
        print("  0. Exit")
        choice = input("\n  Select: ").strip()

        if choice == "0":
            break

        elif choice == "3":
            record_calibration()
            input("\n  Press Enter...")

        elif choice in ("1", "2"):
            dur = input("\n  Duration [60]: ").strip()
            duration = int(dur) if dur.isdigit() else 60

            cal_choice = select_calibration()
            cal_file   = record_calibration() if cal_choice == "NEW" else cal_choice

            timestamp  = datetime.now().strftime('%Y%m%d_%H%M%S')
            label      = 'obs' if choice == '1' else 'loop'
            output_dir = os.path.join(OUTPUT_BASE, f'{label}_{timestamp}')
            os.makedirs(output_dir, exist_ok=True)

            if choice == "1":
                run_observation(duration, cal_file, output_dir, 'obs_0001')
                input("\n  Press Enter...")
            else:
                import time
                count = 0
                print("\n  Running — Ctrl+C to stop\n")
                try:
                    while True:
                        count += 1
                        run_observation(duration, cal_file, output_dir, f'obs_{count:04d}')
                        time.sleep(5)
                except KeyboardInterrupt:
                    print(f"\n  Stopped after {count} observations")
                    input("\n  Press Enter...")


if __name__ == "__main__":
    main()
