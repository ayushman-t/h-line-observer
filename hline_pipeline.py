#!/usr/bin/env python3
"""
H-Line Pipeline
================
Processes local radio observatory data with Virgo and generates:
- 3D plots with RA axis (supports multi-day combined plots)
- Flipbook animations
- Single observation plots
- ezRA pipeline integration

Usage:
    python3 hline_pipeline.py --help
    python3 hline_pipeline.py --data-dir ~/hydrogen_obs --list
    python3 hline_pipeline.py --process loop_20260209 --cal /path/to/calibration.dat
    python3 hline_pipeline.py --plot3d
    python3 hline_pipeline.py --all loop_20260209 --cal /path/to/calibration.dat
"""

import os
import sys
import glob
import csv
import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
try:
    import plotly.graph_objects as go
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
from datetime import datetime, timedelta

# ============== CONFIGURATION ==============

# Local paths
LOCAL_BASE = os.path.expanduser("~/hydrogen_obs_pipeline")
LOCAL_RAW = os.path.join(LOCAL_BASE, "raw")
LOCAL_PROCESSED = os.path.join(LOCAL_BASE, "processed")
LOCAL_OUTPUT = os.path.join(LOCAL_BASE, "output")
LOCAL_CAL = os.path.join(LOCAL_BASE, "calibration")

# Calibration file
CAL_FILE = os.path.join(LOCAL_CAL, "calibration.dat")

# SDR type
SDR_TYPE = "Airspy"

# Observatory location — update to your location
OBS_LAT = 0.0
OBS_LON = 0.0
OBS_ALT = 0  # meters

# Dish pointing — update to match your setup
DISH_AZ = 180
DISH_EL = 45

# Frequency range for H-line (matching ezRA working format)
FREQ_MIN = 1419.5
FREQ_MAX = 1421.0
F_REST = 1420.405751  # MHz
C_KMS = 299792.458

# IST offset
IST_OFFSET = timedelta(hours=5, minutes=30)

# ============== LOCAL DATA FUNCTIONS ==============
def list_local_loops(data_dir, limit=10):
    """List available loop folders in local data directory"""
    loops = sorted(glob.glob(os.path.join(data_dir, "loop_*")), reverse=True)
    if loops:
        loops = loops[:limit]
        print(f"\nFound {len(loops)} loop folders in {data_dir} (showing {len(loops)}):")
        for loop in loops:
            dat_files = glob.glob(os.path.join(loop, "*_observation.dat"))
            print(f"  {os.path.basename(loop)} ({len(dat_files)} .dat files)")
        return [os.path.basename(l) for l in loops]
    else:
        print(f"No loop folders found in {data_dir}")
        return []

# ============== PROCESSING FUNCTIONS ==============
def reprocess_loop(loop_path, cal_file=CAL_FILE):
    """Reprocess .dat files using Virgo to get median filtered data"""
    import virgo
    
    # Check if cal file exists locally, if not fetch it
    if not os.path.exists(cal_file):
        print(f"ERROR: Calibration file not found: {cal_file}")
        print("  Use --cal /path/to/calibration.dat to specify your calibration file")
        return None
    
    loop_name = os.path.basename(loop_path)
    output_dir = os.path.join(LOCAL_PROCESSED, loop_name)
    os.makedirs(output_dir, exist_ok=True)
    
    obs_files = sorted(glob.glob(f"{loop_path}/obs_*_observation.dat"))
    
    if not obs_files:
        print(f"ERROR: No observation .dat files found in {loop_path}")
        return None
    
    print(f"Processing {len(obs_files)} observations with Virgo...")
    print(f"Cal file: {cal_file}")
    print(f"SDR: {SDR_TYPE}")
    
    # Observation parameters to pass directly (avoids header parsing issues)
    obs_params = {
        'dev_args': 'airspy=0,bias=1',
        'rf_gain': 15,
        'if_gain': 15,
        'bb_gain': 15,
        'frequency': 1420405750.0,
        'bandwidth': 3e6,
        'channels': 1024,
        't_sample': 1,
        'duration': 300,
        'loc': (OBS_LAT, OBS_LON, OBS_ALT),
        'ra_dec': '',
        'az_alt': ''
    }
    
    success_count = 0
    for i, obs_file in enumerate(obs_files):
        basename = os.path.basename(obs_file).replace("_observation.dat", "")
        out_csv = f"{output_dir}/{basename}_spectra_filtered.csv"
        out_plot = f"{output_dir}/{basename}_plot.png"
        
        if os.path.exists(out_csv):
            success_count += 1
            continue
        
        print(f"[{i+1}/{len(obs_files)}] {basename}...", end=" ", flush=True)
        
        try:
            virgo.plot(
                obs_parameters=obs_params,
                obs_file=obs_file,
                cal_file=cal_file,
                spectra_csv=out_csv,
                plot_file=out_plot,
                n=20, m=35,
                f_rest=1420.4057517667e6,
                slope_correction=True,
                dB=True
            )
            print("OK")
            success_count += 1
        except Exception as e:
            print(f"ERROR: {e}")
    
    print(f"\nProcessed {success_count}/{len(obs_files)} observations to: {output_dir}")
    return output_dir

# ============== COORDINATE FUNCTIONS ==============
def datetime_from_filename(filename):
    """Extract datetime from filename like obs_0001_20260209_091534_..."""
    parts = os.path.basename(filename).split('_')
    # obs_0001_20260209_091534_observation.dat
    date_str = parts[2]  # 20260209
    time_str = parts[3]  # 091534
    dt = datetime.strptime(f"{date_str}_{time_str}", "%Y%m%d_%H%M%S")
    return dt  # This is IST

def ist_to_utc(dt_ist):
    """Convert IST datetime to UTC"""
    return dt_ist - IST_OFFSET

def calculate_lst(dt_utc, lon_deg):
    """Calculate Local Sidereal Time in hours"""
    # Julian Date
    jd = 367 * dt_utc.year - int(7 * (dt_utc.year + int((dt_utc.month + 9) / 12)) / 4) + \
         int(275 * dt_utc.month / 9) + dt_utc.day + 1721013.5 + \
         (dt_utc.hour + dt_utc.minute / 60 + dt_utc.second / 3600) / 24
    
    # Days since J2000.0
    d = jd - 2451545.0
    
    # Greenwich Mean Sidereal Time
    gmst = 18.697374558 + 24.06570982441908 * d
    gmst = gmst % 24
    
    # Local Sidereal Time
    lst = gmst + lon_deg / 15.0
    lst = lst % 24
    
    return lst

def altaz_to_radec(alt_deg, az_deg, lat_deg, lst_hours):
    """Convert Alt/Az to RA/DEC"""
    import math
    
    alt = math.radians(alt_deg)
    az = math.radians(az_deg)
    lat = math.radians(lat_deg)
    
    # Calculate DEC
    sin_dec = math.sin(alt) * math.sin(lat) + math.cos(alt) * math.cos(lat) * math.cos(az)
    dec = math.asin(sin_dec)
    
    # Calculate Hour Angle
    cos_ha = (math.sin(alt) - math.sin(lat) * sin_dec) / (math.cos(lat) * math.cos(dec))
    cos_ha = max(-1, min(1, cos_ha))  # Clamp
    ha = math.acos(cos_ha)
    
    if math.sin(az) > 0:
        ha = 2 * math.pi - ha
    
    # Calculate RA
    ra = lst_hours * 15 - math.degrees(ha)
    ra = ra % 360
    
    return (ra / 15.0) % 24, math.degrees(dec)  # RA in hours (0-24), DEC in degrees

def get_observation_coords(filename, az=DISH_AZ, el=DISH_EL):
    """Get RA/DEC for an observation"""
    dt_ist = datetime_from_filename(filename)
    dt_utc = ist_to_utc(dt_ist)
    lst = calculate_lst(dt_utc, OBS_LON)
    ra, dec = altaz_to_radec(el, az, OBS_LAT, lst)
    return ra, dec, dt_ist, dt_utc

# ============== VLSR CORRECTION ==============

def calculate_vlsr(ra_h, dec_deg, dt_utc, lat=OBS_LAT, lon=OBS_LON, alt=OBS_ALT):
    """Calculate VLSR correction in km/s to ADD to topocentric velocity to get LSR.

    astropy's radial_velocity_correction only supports 'barycentric'/'heliocentric',
    so we:
      1. Get barycentric correction from astropy
      2. Add solar motion projected onto the observation direction (IAU standard:
         solar apex at RA=18.0h, DEC=+30.0 deg, speed=20.0 km/s)

    Returns 0.0 if astropy unavailable.
    """
    try:
        from astropy.time import Time
        from astropy.coordinates import SkyCoord, EarthLocation
        import astropy.units as u
        import math

        # Step 1: barycentric correction (observer → barycentre)
        location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=alt * u.m)
        obstime  = Time(dt_utc)
        sc       = SkyCoord(ra=ra_h * 15 * u.deg, dec=dec_deg * u.deg, frame='icrs')
        v_bary   = sc.radial_velocity_correction(kind='barycentric',
                                                  obstime=obstime,
                                                  location=location).to(u.km / u.s).value

        # Step 2: solar motion correction (barycentre → LSR)
        # IAU standard solar apex: RA=18.0h, DEC=+30.0 deg, v_sun=20.0 km/s
        ra_apex  = math.radians(18.0 * 15.0)   # 270 deg
        dec_apex = math.radians(30.0)
        ra_obs   = math.radians(ra_h * 15.0)
        dec_obs  = math.radians(dec_deg)

        # Projection of solar motion onto observation direction
        cos_angle = (math.sin(dec_obs) * math.sin(dec_apex) +
                     math.cos(dec_obs) * math.cos(dec_apex) *
                     math.cos(ra_obs - ra_apex))
        v_solar  = 20.0 * cos_angle   # km/s

        return v_bary + v_solar

    except ImportError:
        return 0.0
    except Exception as e:
        print(f"  VLSR warning: {e}")
        return 0.0


_ASTROPY_CHECKED = None
def check_astropy():
    global _ASTROPY_CHECKED
    if _ASTROPY_CHECKED is not None:
        return _ASTROPY_CHECKED
    try:
        import astropy
        _ASTROPY_CHECKED = True
        return True
    except ImportError:
        print("WARNING: astropy not found -- VLSR disabled.  pip install astropy")
        _ASTROPY_CHECKED = False
        return False


def build_vlsr_grid(data, freq, use_vlsr=True):
    """Shift each obs onto LSR velocity grid and interpolate onto common axis.
    Returns (vel_grid, snr_vlsr_2d, vlsr_vals).
    """
    from scipy.interpolate import interp1d
    vel_topo  = C_KMS * (F_REST - freq) / F_REST
    vlsr_vals = []
    for d in data:
        v = calculate_vlsr(d["ra"], d["dec"], d["dt_utc"]) if use_vlsr else 0.0
        vlsr_vals.append(v)
        d["vlsr"] = v
    vlsr_vals = np.array(vlsr_vals)
    all_lsr   = np.array([vel_topo + v for v in vlsr_vals])
    vel_grid  = np.linspace(all_lsr.max(), all_lsr.min(), len(freq))
    snr_grid  = np.full((len(data), len(vel_grid)), np.nan)
    for i, d in enumerate(data):
        lsr_i    = vel_topo + vlsr_vals[i]
        sidx     = np.argsort(lsr_i)
        fi       = interp1d(lsr_i[sidx], d["snr"][sidx], kind="linear",
                            bounds_error=False, fill_value=np.nan)
        snr_grid[i] = fi(vel_grid)
    return vel_grid, snr_grid, vlsr_vals



def baseline_correct_spectrum(freq, arr):
    """Subtract baseline using median of explicitly line-free channels only.

    The naive np.median(arr) includes H-line channels — when H-line is bright
    it pulls the median up, making adjacent channels go artificially negative
    (the dark band artifact in heatmaps/3D plots at high-emission RA).

    Standard radio-astronomy fix: compute the baseline from channels that are
    clearly outside the H-line emission window only (line-free channels).
    We exclude F_REST ± 1.5 MHz which covers the full galactic H-line width.
    This is a scalar (0th-order) subtraction — robust to RFI in any channel.
    """
    # Exclude the H-line emission window from the median estimate.
    # Band is 1419.8-1421.5 MHz, F_REST=1420.406.
    # ±0.5 MHz exclusion leaves ~240 line-free bins (low side: 1419.8-1419.9, high side: 1420.9-1421.5).
    # This prevents the H-line pulling the median up at bright-emission RA positions,
    # which would cause artificially negative values in surrounding channels (dark band).
    lo_bound = 1420.006  # F_REST - 0.4 MHz (~84 km/s, covers galactic emission)
    hi_bound = 1420.806  # F_REST + 0.4 MHz
    line_free = (freq < lo_bound) | (freq > hi_bound)

    if line_free.sum() >= 10:
        baseline = np.median(arr[line_free])
    else:
        baseline = np.median(arr)  # fallback (shouldn't happen with our band)

    return arr - baseline

# ============== PLOTTING FUNCTIONS ==============
def load_processed_data(processed_dir, baseline_correct=True):
    """Load all processed CSV files"""
    files = sorted(glob.glob(f"{processed_dir}/*_spectra_filtered.csv"))
    
    if not files:
        print(f"ERROR: No processed CSV files found in {processed_dir}")
        return [], None
    
    print(f"Loading {len(files)} processed files...")
    
    all_data = []
    common_freq = None
    
    for f in files:
        freq, filtered = [], []
        with open(f, 'r') as fp:
            for row in csv.reader(fp):
                if len(row) >= 5:
                    fr = float(row[0])
                    if FREQ_MIN <= fr <= FREQ_MAX:
                        freq.append(fr)
                        filtered.append(float(row[4]))  # Column 4 = red line
        
        if not freq:
            continue
        
        if common_freq is None:
            common_freq = np.array(freq)
        
        arr = np.array(filtered)
        if baseline_correct:
            arr = baseline_correct_spectrum(np.array(freq), arr)
        
        ra, dec, dt_ist, dt_utc = get_observation_coords(f)
        hours_from_start = None  # Will calculate later
        
        all_data.append({
            'file': f,
            'freq': common_freq,
            'snr': arr,
            'ra': ra,
            'dec': dec,
            'dt_ist': dt_ist,
            'dt_utc': dt_utc
        })
    
    # Calculate hours from start
    if all_data:
        start_time = all_data[0]['dt_ist']
        for d in all_data:
            d['hours'] = (d['dt_ist'] - start_time).total_seconds() / 3600
    
    return all_data, common_freq

# ============== SETI FUNCTIONS ==============
def load_raw_data(loop_path):
    """Load raw observation.dat files (no VIRGO processing, no calibration)
    This is the pure FFT output from the SDR - for SETI narrowband searches
    Uses same H-line frequency region to avoid bandpass edges
    """
    dat_files = sorted(glob.glob(f"{loop_path}/obs_*_observation.dat"))
    
    if not dat_files:
        print(f"ERROR: No raw observation.dat files found in {loop_path}")
        return [], None
    
    print(f"Loading {len(dat_files)} raw .dat files (no calibration)...")
    
    # Check first file to determine format
    first_file = dat_files[0]
    first_power = np.fromfile(first_file, dtype=np.float32)
    
    # SDR settings (from observation parameters)
    center_freq = 1420.405750  # MHz
    bandwidth = 3.0  # MHz
    
    # Determine number of channels from file
    # If file has N values, it could be:
    # - N channels of one spectrum
    # - Multiple spectra of 1024 channels each
    if len(first_power) == 1024:
        channels = 1024
    elif len(first_power) % 1024 == 0:
        channels = 1024
        print(f"File contains {len(first_power)//1024} spectra of 1024 channels")
    else:
        channels = len(first_power)
        print(f"DEBUG: Using {channels} channels from file")
    
    # Calculate frequency axis
    freq = np.linspace(center_freq - bandwidth/2, center_freq + bandwidth/2, channels)
    
    # Filter to H-line region
    mask = (freq >= FREQ_MIN) & (freq <= FREQ_MAX)
    freq_filtered = freq[mask]
    
    print(f"Frequency range: {freq_filtered[0]:.3f} - {freq_filtered[-1]:.3f} MHz ({len(freq_filtered)} channels)")
    
    all_data = []
    
    for f in dat_files:
        try:
            # Raw .dat file contains only power values (float32)
            power = np.fromfile(f, dtype=np.float32)
            
            # Check if it matches expected channels
            if len(power) != channels:
                # Some files might have time-averaged data (multiple spectra)
                if len(power) % channels == 0:
                    # Average all spectra
                    n_spectra = len(power) // channels
                    power = power.reshape(n_spectra, channels).mean(axis=0)
                else:
                    print(f"Warning: {f} has {len(power)} values, expected {channels}")
                    continue
            
            # Apply frequency filter
            power_filtered = power[mask]
            
            ra, dec, dt_ist, dt_utc = get_observation_coords(f)
            
            all_data.append({
                'file': f,
                'freq': freq_filtered,
                'power': power_filtered,  # Raw power, not SNR
                'ra': ra,
                'dec': dec,
                'dt_ist': dt_ist,
                'dt_utc': dt_utc
            })
        except Exception as e:
            print(f"Error loading {f}: {e}")
            continue
    
    # Calculate hours from start
    if all_data:
        start_time = all_data[0]['dt_ist']
        for d in all_data:
            d['hours'] = (d['dt_ist'] - start_time).total_seconds() / 3600
        print(f"Loaded {len(all_data)} observations")
    
    return all_data, freq_filtered

def load_calibrated_data(processed_dir):
    """Load VIRGO processed data WITHOUT median subtraction
    This is calibrated but not baseline-corrected - better for SETI than raw,
    as it removes bandpass shape but keeps all signals
    """
    return load_processed_data(processed_dir, baseline_correct=False)

def plot_seti_waterfall(data_path, data_type='calibrated', view='2d', format='show', x_axis='ra', save_path=None):
    """Generate SETI waterfall plot (time vs frequency, intensity as color)
    
    data_type: 'raw' = raw observation.dat (no processing)
               'calibrated' = VIRGO processed (no median subtraction)
    view: '2d'/'top' = flat heatmap, '3d' = 3D surface
    format: 'show' = matplotlib, 'html' = plotly
    x_axis: 'ra' = Right Ascension, 'time' = hours from start
    """
    # Handle 'top' as alias for '2d'
    if view == 'top':
        view = '2d'
    
    if data_type == 'raw':
        data, freq = load_raw_data(data_path)
        power_key = 'power'
        title_suffix = "Raw FFT (No Calibration)"
    else:
        data, freq = load_calibrated_data(data_path)
        power_key = 'snr'
        title_suffix = "Calibrated (No Baseline)"
    
    if not data:
        print("No data to plot!")
        return
    
    print(f"Plotting SETI waterfall: {len(data)} observations")
    
    # Build waterfall array
    all_power = []
    all_x = []
    
    for d in data:
        all_power.append(d[power_key])
        if x_axis == 'ra':
            all_x.append(d['ra'])
        else:
            all_x.append(d['hours'])
    
    all_power = np.array(all_power)
    all_x = np.array(all_x)
    
    # Sort by x-axis
    sort_idx = np.argsort(all_x)
    all_power = all_power[sort_idx]
    all_x = all_x[sort_idx]
    
    x_label = 'Right Ascension (hours)' if x_axis == 'ra' else 'Time (hours)'
    
    if view == '3d' and format == 'html':
        # Plotly 3D
        if not PLOTLY_AVAILABLE:
            print("ERROR: Plotly not installed! Run: pip install plotly")
            return
        
        X, Y = np.meshgrid(all_x, freq, indexing='ij')
        
        fig = go.Figure(data=[go.Surface(
            x=X, y=Y, z=all_power,
            colorscale='Viridis',
            colorbar=dict(title=dict(text='Power', side='right')),
            lighting=dict(ambient=0.95, diffuse=0.8, specular=0.05, fresnel=0.05),
        )])
        
        fig.update_layout(
            title=dict(text=f'SETI 3D - {title_suffix}', x=0.5, xanchor='center'),
            scene=dict(
                xaxis=dict(title=x_label),
                yaxis=dict(title='Frequency (MHz)'),
                zaxis=dict(title='Power'),
                camera=dict(eye=dict(x=1.5, y=1.5, z=0.8)),
                aspectmode='manual',
                aspectratio=dict(x=1, y=1, z=1),
            ),
            width=1200, height=800,
        )
        
        if save_path is None:
            save_path = os.path.join(LOCAL_OUTPUT, f"seti_3d_{data_type}.html")
        fig.write_html(save_path, include_plotlyjs=True, full_html=True)
        print(f"Saved: {save_path}")
        print("Open in browser for interactive 3D!")
        
    elif view == '3d':
        # Matplotlib 3D
        from mpl_toolkits.mplot3d import Axes3D
        
        X, Y = np.meshgrid(all_x, freq, indexing='ij')
        
        fig = plt.figure(figsize=(14, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        # Decimate for live view
        if save_path:
            rstride, cstride = 1, 1
        else:
            rstride = max(1, len(all_x) // 60)
            cstride = max(1, len(freq) // 60)
        
        ax.plot_surface(X, Y, all_power, cmap='viridis', alpha=0.9,
                        rstride=rstride, cstride=cstride, linewidth=0, edgecolor='none', antialiased=True)
        
        ax.set_xlabel(x_label)
        ax.set_ylabel('Frequency (MHz)')
        ax.set_zlabel('Power')
        ax.set_title(f'SETI 3D - {title_suffix}')
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved: {save_path}")
        else:
            plt.show()
    
    else:
        # 2D waterfall (default)
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # Use pcolormesh for waterfall
        im = ax.pcolormesh(freq, all_x, all_power, shading='auto', cmap='viridis')
        
        ax.set_xlabel('Frequency (MHz)')
        ax.set_ylabel(x_label)
        ax.set_title(f'SETI Waterfall - {title_suffix}')
        
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Power (dB)' if data_type == 'calibrated' else 'Raw Power')
        
        # Mark H-line rest frequency
        ax.axvline(x=F_REST, color='red', linestyle='--', alpha=0.5, label=f'H-line ({F_REST} MHz)')
        ax.legend(loc='upper right')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved: {save_path}")
        else:
            plt.show()

def seti_find_candidates(data_path, data_type='calibrated', sigma=5):
    """Find narrowband signal candidates (>sigma threshold)
    
    Returns list of (time, frequency, power, snr_above_mean) tuples
    """
    if data_type == 'raw':
        data, freq = load_raw_data(data_path)
        power_key = 'power'
    else:
        data, freq = load_calibrated_data(data_path)
        power_key = 'snr'
    
    if not data:
        print("No data!")
        return []
    
    print(f"Searching for narrowband candidates (>{sigma} sigma)...")
    
    candidates = []
    
    for d in data:
        power = d[power_key]
        obs_freq = d.get('freq', freq)  # Use observation's freq if available, else common
        mean_p = np.mean(power)
        std_p = np.std(power)
        threshold = mean_p + sigma * std_p
        
        # Find peaks above threshold
        peaks = np.where(power > threshold)[0]
        
        for idx in peaks:
            candidates.append({
                'time': d['hours'],
                'dt_ist': d['dt_ist'],
                'ra': d['ra'],
                'dec': d['dec'],
                'freq_mhz': obs_freq[idx] if idx < len(obs_freq) else freq[idx],
                'power': power[idx],
                'snr': (power[idx] - mean_p) / std_p,
                'file': d['file']
            })
    
    print(f"Found {len(candidates)} candidates above {sigma} sigma")
    
    # Filter out H-line region (likely galactic hydrogen, not ET)
    hline_min = F_REST - 0.5  # MHz
    hline_max = F_REST + 0.5
    
    non_hline = [c for c in candidates if not (hline_min <= c['freq_mhz'] <= hline_max)]
    print(f"Excluding H-line region: {len(non_hline)} candidates remain")
    
    # Sort by SNR (highest first)
    non_hline.sort(key=lambda x: x['snr'], reverse=True)
    
    return non_hline

def seti_report(candidates, top_n=20):
    """Print top SETI candidates"""
    if not candidates:
        print("No candidates to report!")
        return
    
    print(f"\n=== Top {min(top_n, len(candidates))} SETI Candidates ===")
    print(f"{'#':<3} {'Time(h)':<8} {'Freq(MHz)':<12} {'SNR':<8} {'RA':<8} {'DEC':<8}")
    print("-" * 55)
    
    for i, c in enumerate(candidates[:top_n]):
        print(f"{i+1:<3} {c['time']:<8.2f} {c['freq_mhz']:<12.6f} {c['snr']:<8.1f} {c['ra']:<8.2f} {c['dec']:<8.1f}")
    
    print(f"\nNote: H-line region ({F_REST-0.5:.2f} - {F_REST+0.5:.2f} MHz) excluded")
    print("High SNR at same freq across multiple times = likely RFI")
    print("Signal that drifts in freq over time = interesting!")

def plot_3d(processed_dir, x_axis='ra', y_axis='freq', baseline_correct=True, render='grid', view='3d', clean=False, quality='fast', save_path=None, use_vlsr=False):
    """Generate 3D plot with selectable X-axis (ra, time) and Y-axis (freq, velocity)"""
    data, freq = load_processed_data(processed_dir, baseline_correct=baseline_correct)
    
    if not data:
        print("No data to plot!")
        return
    
    all_x = []
    all_snr = []
    all_ra = []
    all_dec = []
    
    for d in data:
        if x_axis == 'ra':
            all_x.append(d['ra'])
        else:
            all_x.append(d['hours'])
        all_snr.append(d['snr'])
        all_ra.append(d['ra'])
        all_dec.append(d['dec'])
    
    all_x = np.array(all_x)
    all_snr = np.array(all_snr)
    all_ra = np.array(all_ra)
    all_dec = np.array(all_dec)
    
    # Sort by x-axis to avoid surface interpolation artifacts
    sort_idx = np.argsort(all_x)
    all_x = all_x[sort_idx]
    all_snr = all_snr[sort_idx]
    all_ra = all_ra[sort_idx]
    all_dec = all_dec[sort_idx]
    
    # Y-axis: frequency or velocity
    if y_axis == 'velocity':
        if use_vlsr and check_astropy():
            sorted_data = [data[sort_idx[i]] for i in range(len(sort_idx))]
            y_data, all_snr, vlsr_vals = build_vlsr_grid(sorted_data, freq)
            mean_v = float(np.nanmean(vlsr_vals))
            y_label = f'Velocity LSR (km/s)  [VLSR corr: {mean_v:+.1f} km/s]'
        else:
            y_data = C_KMS * (F_REST - freq) / F_REST
            y_label = 'Velocity (km/s)'
    else:
        y_data = freq
        y_label = 'Frequency (MHz)'

    X, Y = np.meshgrid(all_x, y_data, indexing='ij')
    
    # Top view uses 2D pcolormesh (no polygon artifacts)
    if view == 'top':
        fig = plt.figure(figsize=(14, 8))
        ax = fig.add_subplot(111)
        
        mesh = ax.pcolormesh(X, Y, all_snr, cmap='viridis', shading='auto')
        plt.colorbar(mesh, ax=ax, label='SNR (dB)')
        
        if x_axis == 'ra':
            ax.set_xlabel('Right Ascension (hours)')
            ax.set_xlim(0, 24)
        else:
            ax.set_xlabel('Time (hours)')
        
        ax.set_ylabel(y_label)
        
        # Date range
        start_date = data[0]['dt_ist'].strftime('%Y-%m-%d')
        end_date = data[-1]['dt_ist'].strftime('%Y-%m-%d')
        date_str = start_date if start_date == end_date else f"{start_date} to {end_date}"
        dec = data[0]['dec']
        ax.set_title(f'H-Line Observatory - {date_str}\nDEC: {dec:.1f}°')
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved: {save_path}")
            plt.close()
        else:
            plt.show()
        return
    
    # 3D view
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Decimate for live view based on quality setting
    if save_path:
        rstride, cstride = 1, 1  # Full resolution for saved image
    else:
        # Quality levels: fast=60 (smooth rotation), medium=200, high=350 (detailed)
        quality_map = {'fast': 60, 'medium': 200, 'high': 350}
        target = quality_map.get(quality, 60)
        rstride = max(1, len(all_x) // target)
        cstride = max(1, len(y_data) // target)
    
    if render == 'solid':
        # Solid render - no grid lines
        ax.plot_surface(X, Y, all_snr, cmap='viridis', alpha=0.9, 
                        rstride=rstride, cstride=cstride, linewidth=0, edgecolor='none', antialiased=True)
    else:
        # Grid render - also no grid lines now
        ax.plot_surface(X, Y, all_snr, cmap='viridis', alpha=0.8, 
                        rstride=rstride, cstride=cstride, linewidth=0, edgecolor='none', antialiased=True)
    
    # Remove panes and grids if clean mode
    if clean:
        # Make panes transparent
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        # Make grid lines transparent
        ax.xaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
        ax.yaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
        ax.zaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    
    # Find and label peaks (observations with SNR > 5 above median)
    median_snr = np.median(all_snr)
    threshold = median_snr + 5
    
    for i in range(len(all_x)):
        peak_idx = np.argmax(all_snr[i])
        peak_val = all_snr[i, peak_idx]
        
        if peak_val > threshold:
            peak_freq = freq[peak_idx]
            # Calculate Doppler velocity
            doppler = C_KMS * (F_REST - peak_freq) / F_REST
            
            # Add text annotation (only for strong peaks, every 10th to avoid clutter)
            if i % 10 == 0:
                label = f'({all_ra[i]:.1f}h, {all_dec[i]:.0f}°, {doppler:.0f}km/s)'
                ax.text(all_x[i], y_data[peak_idx], peak_val + 0.5, label, 
                        fontsize=6, ha='center', alpha=0.8)
    
    if x_axis == 'ra':
        ax.set_xlabel('Right Ascension (hours)')
        ax.set_xlim(0, 24)
    else:
        ax.set_xlabel('Time (hours)')
        ax.set_xlim(all_x.min(), all_x.max())
    
    ax.set_ylabel(y_label)
    ax.set_ylim(y_data.min(), y_data.max())
    
    if baseline_correct:
        ax.set_zlabel('SNR (baseline corrected)')
    else:
        ax.set_zlabel('SNR (raw)')
    
    # Date range (single loop can span multiple days)
    start_date = data[0]['dt_ist'].strftime('%Y-%m-%d')
    end_date = data[-1]['dt_ist'].strftime('%Y-%m-%d')
    date_str = start_date if start_date == end_date else f"{start_date} to {end_date}"
    dec = data[0]['dec']  # DEC is constant for fixed dish
    ax.set_title(f'H-Line Observatory - {date_str}\nDEC: {dec:.1f}°')
    
    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"Saved: {save_path}")
        plt.close()
    else:
        plt.show()

def plot_3d_plotly(processed_dir, x_axis='ra', y_axis='freq', baseline_correct=True, save_path=None, use_vlsr=False):
    """Generate interactive 3D plot using Plotly (saves as HTML)"""
    if not PLOTLY_AVAILABLE:
        print("ERROR: Plotly not installed! Run: pip install plotly")
        return
    
    data, freq = load_processed_data(processed_dir, baseline_correct=baseline_correct)
    
    if not data:
        print("No data to plot!")
        return
    
    all_x = []
    all_snr = []
    all_ra = []
    all_dec = []
    
    for d in data:
        if x_axis == 'ra':
            all_x.append(d['ra'])
        else:
            all_x.append(d['hours'])
        all_snr.append(d['snr'])
        all_ra.append(d['ra'])
        all_dec.append(d['dec'])
    
    all_x = np.array(all_x)
    all_snr = np.array(all_snr)
    all_ra = np.array(all_ra)
    all_dec = np.array(all_dec)
    
    # Sort by x-axis
    sort_idx = np.argsort(all_x)
    all_x = all_x[sort_idx]
    all_snr = all_snr[sort_idx]
    all_ra = all_ra[sort_idx]
    all_dec = all_dec[sort_idx]
    
    print(f"Plotly grid: {len(all_x)} x {len(freq)} points (full resolution)")
    

    # Y-axis: frequency or velocity
    if y_axis == 'velocity':
        if use_vlsr and check_astropy():
            sorted_data = [data[sort_idx[i]] for i in range(len(sort_idx))]
            y_data, all_snr, vlsr_vals = build_vlsr_grid(sorted_data, freq)
            mean_v = float(__import__('numpy').nanmean(vlsr_vals))
            y_label = f'Velocity LSR (km/s)  [VLSR corr: {mean_v:+.1f} km/s]'
        else:
            y_data = C_KMS * (F_REST - freq) / F_REST
            y_label = 'Velocity (km/s)'
    else:
        y_data = freq
        y_label = 'Frequency (MHz)'

    X, Y = np.meshgrid(all_x, y_data, indexing='ij')
    doppler_arr = C_KMS * (F_REST - freq) / F_REST
    
    # Build customdata for fast hover
    ra_grid = np.broadcast_to(all_ra[:, np.newaxis], all_snr.shape)
    dec_grid = np.broadcast_to(all_dec[:, np.newaxis], all_snr.shape)
    vel_grid = np.broadcast_to(doppler_arr[np.newaxis, :], all_snr.shape)
    freq_grid = np.broadcast_to(freq[np.newaxis, :], all_snr.shape)
    customdata = np.stack([ra_grid, dec_grid, vel_grid, freq_grid], axis=-1)
    
    x_label = 'Right Ascension (hours)' if x_axis == 'ra' else 'Time (hours)'
    z_label = 'SNR (dB)'
    
    # Date range (single loop can span multiple days)
    start_date = data[0]['dt_ist'].strftime('%Y-%m-%d')
    end_date = data[-1]['dt_ist'].strftime('%Y-%m-%d')
    date_str = start_date if start_date == end_date else f"{start_date} to {end_date}"
    dec = data[0]['dec']
    
    fig = go.Figure(data=[go.Surface(
        x=X, y=Y, z=all_snr,
        colorscale='Viridis',
        colorbar=dict(title=dict(text='SNR (dB)', side='right')),
        customdata=customdata,
        hovertemplate='<b>RA:</b> %{customdata[0]:.2f}h | <b>DEC:</b> %{customdata[1]:.1f}°<br>' +
                      '<b>Freq:</b> %{customdata[3]:.3f} MHz<br>' +
                      '<b>Velocity:</b> %{customdata[2]:.0f} km/s<br>' +
                      '<b>SNR:</b> %{z:.2f} dB<extra></extra>',
        lighting=dict(ambient=0.95, diffuse=0.8, specular=0.05, fresnel=0.05),
        contours=dict(
            z=dict(show=False),
            x=dict(show=False),
            y=dict(show=False),
        ),
    )])
    
    fig.update_layout(
        title=dict(
            text=f'<b>H-Line Observatory</b><br>{date_str} | DEC: {dec:.1f}°',
            x=0.5, xanchor='center'
        ),
        scene=dict(
            xaxis=dict(title=x_label, gridcolor='rgba(128,128,128,0.3)', showbackground=True, backgroundcolor='rgb(250,250,250)'),
            yaxis=dict(title=y_label, gridcolor='rgba(128,128,128,0.3)', showbackground=True, backgroundcolor='rgb(250,250,250)'),
            zaxis=dict(title=z_label, gridcolor='rgba(128,128,128,0.3)', showbackground=True, backgroundcolor='rgb(250,250,250)'),
            camera=dict(eye=dict(x=1.5, y=1.5, z=0.8)),
            aspectmode='manual',
            aspectratio=dict(x=1, y=1, z=1),  # Equal cube aspect
        ),
        width=1200,
        height=800,
        margin=dict(l=10, r=10, t=60, b=10),
    )
    
    if save_path is None:
        save_path = os.path.join(LOCAL_OUTPUT, "hline_3d_interactive.html")
    
    fig.write_html(save_path, include_plotlyjs=True, full_html=True)
    print(f"Saved: {save_path}")
    print("Open in browser for interactive 3D!")

def plot_heatmap(processed_dir, x_axis='ra', y_axis='freq', baseline_correct=True, save_path=None, use_vlsr=False):
    """Generate 2D heatmap - cleaner than 3D for survey data"""
    data, freq = load_processed_data(processed_dir, baseline_correct=baseline_correct)
    
    if not data:
        print("No data to plot!")
        return
    
    all_x = []
    all_snr = []
    
    for d in data:
        if x_axis == 'ra':
            all_x.append(d['ra'])
        else:
            all_x.append(d['hours'])
        all_snr.append(d['snr'])
    
    all_x = np.array(all_x)
    all_snr = np.array(all_snr)
    
    # Sort by x-axis
    sort_idx = np.argsort(all_x)
    all_x = all_x[sort_idx]
    all_snr = all_snr[sort_idx]
    
    # Y-axis: frequency or velocity
    if y_axis == 'velocity':
        if use_vlsr and check_astropy():
            sorted_data = [data[sort_idx[i]] for i in range(len(sort_idx))]
            y_data, all_snr, vlsr_vals = build_vlsr_grid(sorted_data, freq)
            mean_v = float(__import__('numpy').nanmean(vlsr_vals))
            y_label = f'Velocity LSR (km/s)  [VLSR corr: {mean_v:+.1f} km/s]'
        else:
            y_data = C_KMS * (F_REST - freq) / F_REST
            y_label = 'Velocity (km/s)'
    else:
        y_data = freq
        y_label = 'Frequency (MHz)'
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # pcolormesh for clean rendering
    im = ax.pcolormesh(all_x, y_data, all_snr.T, shading='auto', cmap='viridis')
    
    cbar = plt.colorbar(im, ax=ax, label='SNR (dB)')
    
    if x_axis == 'ra':
        ax.set_xlabel('Right Ascension (hours)')
        ax.set_xlim(0, 24)
    else:
        ax.set_xlabel('Time (hours)')
    
    ax.set_ylabel(y_label)
    
    # Date range (single loop can span multiple days)
    start_date = data[0]['dt_ist'].strftime('%Y-%m-%d')
    end_date = data[-1]['dt_ist'].strftime('%Y-%m-%d')
    date_str = start_date if start_date == end_date else f"{start_date} to {end_date}"
    dec = data[0]['dec']  # DEC is constant for fixed dish
    baseline_str = "baseline corrected" if baseline_correct else "raw"
    ax.set_title(f'H-Line Observatory - {date_str} ({baseline_str})\nDEC: {dec:.1f}°')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"Saved: {save_path}")
        plt.close()
    else:
        plt.show()

def plot_flipbook(processed_dir, y_axis='freq', baseline_correct=True, save_path=None, fps=5, use_vlsr=False):
    """Generate flipbook animation"""
    from matplotlib.animation import FuncAnimation, PillowWriter
    
    data, freq = load_processed_data(processed_dir, baseline_correct=baseline_correct)
    
    if not data:
        print("No data!")
        return
    
    # Y-axis: frequency or velocity
    if y_axis == 'velocity':
        if use_vlsr and check_astropy():
            # Pre-compute per-obs LSR velocity (shifted inline during animation)
            _vlsr_shifts = [calculate_vlsr(d['ra'], d['dec'], d['dt_utc']) for d in data]
        y_data = C_KMS * (F_REST - freq) / F_REST
        y_label = 'Velocity LSR (km/s)' if (use_vlsr and check_astropy()) else 'Velocity (km/s)'
    else:
        _vlsr_shifts = None
        y_data = freq
        y_label = 'Frequency (MHz)'
    
    # Find y limits for SNR
    all_snr = np.array([d['snr'] for d in data])
    snr_min, snr_max = all_snr.min() - 0.5, all_snr.max() + 0.5
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    def animate(i):
        ax.clear()
        ax.plot(y_data, data[i]['snr'], 'r-', lw=1.5)
        if y_axis == 'velocity':
            ax.axvline(x=0, color='g', linestyle='--', alpha=0.5, label='0 km/s (rest)')
        else:
            ax.axvline(x=F_REST, color='g', linestyle='--', alpha=0.5, label=f'{F_REST:.2f} MHz (rest)')
        ax.set_xlim(y_data[-1], y_data[0])
        ax.set_ylim(snr_min, snr_max)
        ax.set_xlabel(y_label)
        ax.set_ylabel('SNR')
        ax.set_title(f'H-Line Observatory | RA: {data[i]["ra"]:.2f}h | DEC: {data[i]["dec"]:.1f}° | {data[i]["dt_ist"].strftime("%Y-%m-%d %H:%M")} IST')
        ax.legend(loc='upper right')
        ax.grid(True, alpha=0.3)
    
    anim = FuncAnimation(fig, animate, frames=len(data), interval=1000//fps)
    
    if save_path is None:
        save_path = os.path.join(LOCAL_OUTPUT, "hline_flipbook.gif")
    
    save_dir = os.path.dirname(save_path)
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    
    anim.save(save_path, writer=PillowWriter(fps=fps))
    print(f"Saved: {save_path}")
    plt.close()

def plot_single(processed_dir, hour):
    """Plot single observation closest to specified hour"""
    data, freq = load_processed_data(processed_dir)
    
    if not data:
        print("No data!")
        return
    
    # Find closest to requested hour
    closest = min(data, key=lambda d: abs(d['hours'] - hour))
    
    velocity = C_KMS * (F_REST - freq) / F_REST
    
    plt.figure(figsize=(12, 6))
    plt.plot(velocity, closest['snr'], 'r-', lw=1.5)
    plt.axvline(x=0, color='g', linestyle='--', alpha=0.5, label='0 km/s (rest)')
    plt.xlabel('Velocity (km/s)')
    plt.ylabel('SNR')
    plt.title(f'H-Line Observatory | RA: {closest["ra"]:.2f}h | DEC: {closest["dec"]:.1f}° | Hour {closest["hours"]:.1f} | {closest["dt_ist"].strftime("%Y-%m-%d %H:%M")} IST')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.gca().invert_xaxis()
    plt.show()

def list_hours(processed_dir):
    """List available observation hours"""
    data, _ = load_processed_data(processed_dir, baseline_correct=False)
    
    if not data:
        print("No data available!")
        return
    
    print(f"\nAvailable observations ({len(data)} total):")
    print("-" * 60)
    print(f"{'Hour':>6} | {'RA (h)':>8} | {'Time (IST)':>20}")
    print("-" * 60)
    
    for d in data[::10]:  # Show every 10th
        print(f"{d['hours']:>6.1f} | {d['ra']:>8.2f} | {d['dt_ist'].strftime('%Y-%m-%d %H:%M')}")
    
    print("-" * 60)
    print(f"Range: Hour 0 to {data[-1]['hours']:.1f}")

# ============== MULTI-DAY FUNCTIONS ==============
def load_multiday_data(processed_dirs, baseline_correct=True):
    """Load and combine data from multiple processed directories
    
    Returns:
        all_data: list of observation dicts with 'loop_idx' and 'loop_name' added
        common_freq: frequency array
        loop_info: list of dicts with start_hour, end_hour, date, name for each loop
    """
    all_data = []
    common_freq = None
    loop_info = []
    
    for loop_idx, processed_dir in enumerate(processed_dirs):
        data, freq = load_processed_data(processed_dir, baseline_correct)
        if data:
            loop_name = os.path.basename(processed_dir)
            loop_date = data[0]['dt_ist'].strftime('%Y-%m-%d')
            
            # Add loop tracking to each observation
            for d in data:
                d['loop_idx'] = loop_idx
                d['loop_name'] = loop_name
            
            all_data.extend(data)
            if common_freq is None:
                common_freq = freq
    
    if not all_data:
        return [], None, []
    
    # Sort by datetime
    all_data.sort(key=lambda d: d['dt_ist'])

    # Assign hours SEQUENTIALLY — loops are placed back-to-back with no gaps.
    # This prevents blank stretches and pcolormesh smearing when days are missing.
    # Each loop's observations are spaced by their real internal timing,
    # but the gap between loops is collapsed to zero.
    seq_hour = 0.0
    for loop_idx, processed_dir in enumerate(processed_dirs):
        loop_obs = sorted(
            [d for d in all_data if d['loop_idx'] == loop_idx],
            key=lambda d: d['dt_ist']
        )
        if not loop_obs:
            continue
        loop_start = loop_obs[0]['dt_ist']
        for d in loop_obs:
            d['hours'] = seq_hour + (d['dt_ist'] - loop_start).total_seconds() / 3600
        # Advance seq_hour by the actual duration of this loop
        loop_duration = (loop_obs[-1]['dt_ist'] - loop_start).total_seconds() / 3600
        seq_hour += loop_duration + 0.1  # tiny gap so dashed boundary lines show

    # Build loop_info with hour boundaries
    for loop_idx, processed_dir in enumerate(processed_dirs):
        loop_obs = [d for d in all_data if d['loop_idx'] == loop_idx]
        if loop_obs:
            loop_info.append({
                'idx': loop_idx,
                'name': os.path.basename(processed_dir),
                'date': loop_obs[0]['dt_ist'].strftime('%Y-%m-%d'),
                'start_hour': loop_obs[0]['hours'],
                'end_hour': loop_obs[-1]['hours'],
                'start_ra': min(d['ra'] for d in loop_obs),
                'end_ra':   max(d['ra'] for d in loop_obs),
                'n_obs': len(loop_obs)
            })
    
    print(f"Combined {len(all_data)} observations from {len(processed_dirs)} loops")
    print(f"Time span: {all_data[0]['dt_ist']} to {all_data[-1]['dt_ist']}")
    for info in loop_info:
        print(f"  {info['date']}: {info['n_obs']} obs, RA {info['start_ra']:.1f}h - {info['end_ra']:.1f}h")
    
    return all_data, common_freq, loop_info

def plot_3d_multiday(processed_dirs, y_axis='freq', baseline_correct=True, render='grid', view='3d', clean=False, quality='fast', save_path=None, use_vlsr=False):
    """Generate 3D plot for multi-day data with date+RA labels"""
    data, freq, loop_info = load_multiday_data(processed_dirs, baseline_correct=baseline_correct)
    
    if not data:
        print("No data to plot!")
        return
    
    # Build hours/snr arrays with NaN sentinel rows at loop boundaries.
    # Without this, matplotlib draws surface triangles between the last obs of
    # loop N (RA ~20h, H-line at one velocity) and first obs of loop N+1
    # (RA ~5h, H-line at different velocity), creating comb-like spike artifacts.
    all_hours = []
    all_snr = []
    n_freq = len(freq)
    prev_loop_idx = data[0]['loop_idx']

    for d in data:
        if d['loop_idx'] != prev_loop_idx:
            # Insert a NaN row to break the surface at loop boundary
            all_hours.append(all_hours[-1] + 1e-6)  # tiny time step
            all_snr.append(np.full(n_freq, np.nan))
            prev_loop_idx = d['loop_idx']
        all_hours.append(d['hours'])
        all_snr.append(d['snr'])

    all_hours = np.array(all_hours)
    all_snr = np.array(all_snr)

    # Y-axis: frequency or velocity
    if y_axis == 'velocity':
        if use_vlsr and check_astropy():
            y_data, all_snr, vlsr_vals = build_vlsr_grid(data, freq)
            mean_v = float(np.nanmean(vlsr_vals))
            y_label = f'Velocity LSR (km/s)  [VLSR corr: {mean_v:+.1f} km/s]'
        else:
            y_data = C_KMS * (F_REST - freq) / F_REST
            y_label = 'Velocity (km/s)'
    else:
        y_data = freq
        y_label = 'Frequency (MHz)'

    X, Y = np.meshgrid(all_hours, y_data, indexing='ij')

    # Build date legend string
    def _ra_label(info):
        span = info['end_ra'] - info['start_ra']
        if span > 20:
            return f"{info['date']} (full sky)"
        return f"{info['date']} (RA {info['start_ra']:.1f}-{info['end_ra']:.1f}h)"
    date_legend = " | ".join([_ra_label(info) for info in loop_info])

    # Top view — build fresh arrays directly from data (no NaN rows from 3D fix)
    # Use shading='nearest' so each cell is colored independently, no interpolation
    if view == 'top':
        top_hours = np.array([d['hours'] for d in data])
        top_snr   = np.array([d['snr']   for d in data])
        fig = plt.figure(figsize=(16, 8))
        ax = fig.add_subplot(111)
        mesh = ax.pcolormesh(top_hours, y_data, top_snr.T, cmap='viridis', shading='nearest')
        plt.colorbar(mesh, ax=ax, label='SNR (dB)')
        
        ax.set_xlabel('Time (hours from start)')
        ax.set_ylabel(y_label)
        
        # Add small tick marks + date labels at top only — no full dividing lines
        y_min, y_max = y_data.min(), y_data.max()
        tick_height = (y_max - y_min) * 0.06  # top 6% of plot only
        for info in loop_info[1:]:  # Skip first loop
            ax.plot([info['start_hour'], info['start_hour']],
                    [y_max - tick_height, y_max],
                    color='white', linestyle='-', alpha=0.8, linewidth=1.2)
            d = info['date']
            ax.text(info['start_hour'] + (loop_info[1]['start_hour'] - loop_info[0]['start_hour']) * 0.01,
                    y_max, f" {d}",
                    color='white', fontsize=7, va='top', ha='left', rotation=90)
        
        # Title with all dates
        ax.set_title(f'H-Line Observatory | DEC: {data[0]["dec"]:.1f}°\n{date_legend}', fontsize=10)
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved: {save_path}")
            plt.close()
        else:
            plt.show()
        return
    
    # 3D view — plot EACH LOOP as a separate surface so matplotlib never
    # tries to draw triangles between loops, eliminating edge comb artifacts.
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(111, projection='3d')

    alpha = 0.9 if render == 'solid' else 0.8

    for info in loop_info:
        loop_obs = [d for d in data if d['loop_idx'] == info['idx']]
        if not loop_obs:
            continue
        loop_hours = np.array([d['hours'] for d in loop_obs])
        loop_snr   = np.array([d['snr']   for d in loop_obs])

        if y_axis == 'velocity' and use_vlsr and check_astropy():
            # snr already vlsr-corrected into all_snr; slice out this loop's rows
            # find indices of this loop in the original data list
            loop_indices = [i for i, d in enumerate(data) if d['loop_idx'] == info['idx']]
            loop_snr = all_snr[loop_indices]

        if save_path:
            rstride, cstride = 1, 1
        else:
            quality_map = {'fast': 60, 'medium': 200, 'high': 350}
            target = quality_map.get(quality, 60)
            rstride = max(1, len(loop_hours) // target)
            cstride = max(1, len(y_data) // target)

        Xl, Yl = np.meshgrid(loop_hours, y_data, indexing='ij')
        ax.plot_surface(Xl, Yl, loop_snr, cmap='viridis', alpha=alpha,
                        rstride=rstride, cstride=cstride,
                        linewidth=0, edgecolor='none', antialiased=False)
    
    # Remove panes and grids if clean mode
    if clean:
        # Make panes transparent
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        # Make grid lines transparent
        ax.xaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
        ax.yaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
        ax.zaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    
    ax.set_xlabel('Time (hours from start)')
    ax.set_ylabel(y_label)
    
    if baseline_correct:
        ax.set_zlabel('SNR (baseline corrected)')
    else:
        ax.set_zlabel('SNR (raw)')
    
    # Title with date range
    start_date = data[0]['dt_ist'].strftime('%Y-%m-%d')
    end_date = data[-1]['dt_ist'].strftime('%Y-%m-%d')
    if start_date == end_date:
        title_date = start_date
    else:
        title_date = f"{start_date} to {end_date}"
    
    ax.set_title(f'H-Line Observatory\n{title_date} | RA: {data[0]["ra"]:.1f}h to {data[-1]["ra"]:.1f}h | DEC: {data[0]["dec"]:.1f}°')
    
    # Add date+RA annotations on the plot
    # Print loop info at bottom
    print("\nLoop breakdown:")
    for info in loop_info:
        print(f"  {info['date']}: Hours {info['start_hour']:.1f}-{info['end_hour']:.1f}, RA {info['start_ra']:.1f}-{info['end_ra']:.1f}h")
    
    ax.set_title(f'H-Line Observatory | DEC: {data[0]["dec"]:.1f}°\n{date_legend}', fontsize=10)
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {save_path}")
        plt.close()
    else:
        plt.show()

def plot_3d_plotly_multiday(processed_dirs, y_axis='freq', baseline_correct=True, save_path=None, use_vlsr=False):
    """Generate interactive 3D plot for multi-day data using Plotly"""
    if not PLOTLY_AVAILABLE:
        print("ERROR: Plotly not installed! Run: pip install plotly")
        return
    
    data, freq, loop_info = load_multiday_data(processed_dirs, baseline_correct=baseline_correct)
    
    if not data:
        print("No data to plot!")
        return
    
    # Insert NaN sentinel rows at loop boundaries to prevent spike artifacts
    all_hours = []
    all_snr = []
    all_ra = []
    all_dec = []
    n_freq = len(freq)
    prev_loop_idx = data[0]['loop_idx']

    for d in data:
        if d['loop_idx'] != prev_loop_idx:
            all_hours.append(all_hours[-1] + 1e-6)
            all_snr.append(np.full(n_freq, np.nan))
            all_ra.append(all_ra[-1])
            all_dec.append(all_dec[-1])
            prev_loop_idx = d['loop_idx']
        all_hours.append(d['hours'])
        all_snr.append(d['snr'])
        all_ra.append(d['ra'])
        all_dec.append(d['dec'])

    all_hours = np.array(all_hours)
    all_snr = np.array(all_snr)
    all_ra = np.array(all_ra)
    all_dec = np.array(all_dec)

    print(f"Plotly grid: {len(all_hours)} x {len(freq)} points")
    
    # Y-axis: frequency or velocity
    if y_axis == 'velocity':
        if use_vlsr and check_astropy():
            y_data, all_snr, vlsr_vals = build_vlsr_grid(data, freq)
            mean_v = float(np.nanmean(vlsr_vals))
            y_label = f'Velocity LSR (km/s)  [VLSR corr: {mean_v:+.1f} km/s]'
        else:
            y_data = C_KMS * (F_REST - freq) / F_REST
            y_label = 'Velocity (km/s)'
    else:
        y_data = freq
        y_label = 'Frequency (MHz)'

    X, Y = np.meshgrid(all_hours, y_data, indexing='ij')
    doppler_arr = C_KMS * (F_REST - freq) / F_REST
    
    # Build customdata for fast hover
    ra_grid = np.broadcast_to(all_ra[:, np.newaxis], all_snr.shape)
    dec_grid = np.broadcast_to(all_dec[:, np.newaxis], all_snr.shape)
    vel_grid = np.broadcast_to(doppler_arr[np.newaxis, :], all_snr.shape)
    freq_grid = np.broadcast_to(freq[np.newaxis, :], all_snr.shape)
    hours_grid = np.broadcast_to(all_hours[:, np.newaxis], all_snr.shape)
    customdata = np.stack([ra_grid, dec_grid, vel_grid, freq_grid, hours_grid], axis=-1)
    
    z_label = 'SNR (dB)'
    
    # Build date legend
    def _ra_label(info):
        span = info['end_ra'] - info['start_ra']
        if span > 20:
            return f"{info['date']} (full sky)"
        return f"{info['date']} (RA {info['start_ra']:.1f}-{info['end_ra']:.1f}h)"
    date_legend = " | ".join([_ra_label(info) for info in loop_info])

    fig = go.Figure(data=[go.Surface(
        x=X, y=Y, z=all_snr,
        colorscale='Viridis',
        colorbar=dict(title=dict(text='SNR (dB)', side='right')),
        customdata=customdata,
        hovertemplate='<b>Time:</b> %{customdata[4]:.1f}h | <b>RA:</b> %{customdata[0]:.2f}h | <b>DEC:</b> %{customdata[1]:.1f}°<br>' +
                      '<b>Freq:</b> %{customdata[3]:.3f} MHz | <b>Vel:</b> %{customdata[2]:.0f} km/s<br>' +
                      '<b>SNR:</b> %{z:.2f} dB<extra></extra>',
        lighting=dict(ambient=0.9, diffuse=0.9, specular=0.1, fresnel=0.1),
    )])
    
    fig.update_layout(
        title=dict(
            text=f'<b>H-Line Observatory</b> | DEC: {data[0]["dec"]:.1f}°<br>{date_legend}',
            x=0.5, xanchor='center'
        ),
        scene=dict(
            xaxis=dict(title='Time (hours)', gridcolor='rgba(128,128,128,0.3)', showbackground=True, backgroundcolor='rgb(250,250,250)'),
            yaxis=dict(title=y_label, gridcolor='rgba(128,128,128,0.3)', showbackground=True, backgroundcolor='rgb(250,250,250)'),
            zaxis=dict(title=z_label, gridcolor='rgba(128,128,128,0.3)', showbackground=True, backgroundcolor='rgb(250,250,250)'),
            camera=dict(eye=dict(x=1.5, y=1.5, z=0.8)),
            aspectmode='manual',
            aspectratio=dict(x=1, y=1, z=1),  # Equal cube aspect
        ),
        width=1200,
        height=800,
        margin=dict(l=10, r=10, t=80, b=10),  # More top margin for longer title
    )
    
    # Print loop breakdown
    print("\nLoop breakdown:")
    for info in loop_info:
        print(f"  {info['date']}: Hours {info['start_hour']:.1f}-{info['end_hour']:.1f}, RA {info['start_ra']:.1f}-{info['end_ra']:.1f}h")
    
    if save_path is None:
        save_path = os.path.join(LOCAL_OUTPUT, "hline_3d_multiday_interactive.html")
    
    fig.write_html(save_path, include_plotlyjs=True, full_html=True)
    print(f"Saved: {save_path}")
    print("Open in browser for interactive 3D!")

def plot_flipbook_multiday(processed_dirs, y_axis='freq', baseline_correct=True, save_path=None, fps=5, use_vlsr=False):
    """Generate flipbook animation for multi-day data"""
    from matplotlib.animation import FuncAnimation, PillowWriter
    
    data, freq, loop_info = load_multiday_data(processed_dirs, baseline_correct=baseline_correct)
    
    if not data:
        print("No data!")
        return
    
    # Y-axis: frequency or velocity
    if y_axis == 'velocity':
        y_data = C_KMS * (F_REST - freq) / F_REST
        y_label = 'Velocity (km/s)'
    else:
        y_data = freq
        y_label = 'Frequency (MHz)'
    
    # Find y limits for SNR
    all_snr = np.array([d['snr'] for d in data])
    snr_min, snr_max = all_snr.min() - 0.5, all_snr.max() + 0.5
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    def animate(i):
        ax.clear()
        if y_axis == 'velocity' and _vlsr_shifts:
            frame_vel = y_data + _vlsr_shifts[i]
            ax.plot(frame_vel, data[i]['snr'], 'r-', lw=1.5)
            ax.axvline(x=0, color='g', linestyle='--', alpha=0.5, label=f'0 km/s LSR  (VLSR={_vlsr_shifts[i]:+.1f})')
            ax.set_xlim(frame_vel[-1], frame_vel[0])
        elif y_axis == 'velocity':
            ax.plot(y_data, data[i]['snr'], 'r-', lw=1.5)
            ax.axvline(x=0, color='g', linestyle='--', alpha=0.5, label='0 km/s (rest)')
            ax.set_xlim(y_data[-1], y_data[0])
        else:
            ax.plot(y_data, data[i]['snr'], 'r-', lw=1.5)
            ax.axvline(x=F_REST, color='g', linestyle='--', alpha=0.5, label=f'{F_REST:.2f} MHz (rest)')
            ax.set_xlim(y_data[-1], y_data[0])
        ax.set_ylim(snr_min, snr_max)
        ax.set_xlabel(y_label)
        ax.set_ylabel('SNR')
        ax.set_title(f'H-Line Observatory | {data[i]["dt_ist"].strftime("%Y-%m-%d %H:%M")} IST | RA: {data[i]["ra"]:.2f}h | DEC: {data[i]["dec"]:.1f}°')
        ax.legend(loc='upper right')
        ax.grid(True, alpha=0.3)
    
    anim = FuncAnimation(fig, animate, frames=len(data), interval=1000//fps)
    
    if save_path is None:
        save_path = os.path.join(LOCAL_OUTPUT, "hline_flipbook_combined.gif")
    
    save_dir = os.path.dirname(save_path)
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    
    anim.save(save_path, writer=PillowWriter(fps=fps))
    print(f"Saved: {save_path}")
    plt.close()

# ============== COMPARE / OVERLAY FUNCTIONS ==============

# Per-loop colors for compare mode (matplotlib colormaps and hex colors)
LOOP_COLORS_MPL  = ['Blues', 'Reds', 'Greens', 'Purples', 'Oranges']
LOOP_COLORS_RGBA = [
    'rgba(30,120,220,{})',    # blue
    'rgba(220,50,50,{})',     # red
    'rgba(40,180,80,{})',     # green
    'rgba(160,60,200,{})',    # purple
    'rgba(230,140,30,{})',    # orange
]
LOOP_COLORS_HEX  = ['#1e78dc', '#dc3232', '#28b450', '#a03cc8', '#e68c1e']

def _build_ra_grid(all_data, loop_info, n_pts=None):
    """Return common RA grid and per-loop interpolated SNR arrays.

    Each loop is interpolated onto the same RA axis so the surfaces
    line up when overlaid.  Observations outside a loop's RA range
    are set to NaN (matplotlib/plotly handle NaN as transparent gaps).
    """
    from scipy.interpolate import interp1d

    # Global RA range from all loops
    all_ra = [d['ra'] for d in all_data]
    ra_min, ra_max = min(all_ra), max(all_ra)

    # Choose resolution: default ~300 pts so plot stays fast
    if n_pts is None:
        n_pts = 300

    ra_grid = np.linspace(ra_min, ra_max, n_pts)

    loop_snr_grids = []
    for info in loop_info:
        obs = [d for d in all_data if d['loop_idx'] == info['idx']]
        obs.sort(key=lambda d: d['ra'])

        loop_ra  = np.array([d['ra']  for d in obs])
        loop_snr = np.array([d['snr'] for d in obs])   # shape (n_obs, n_freq)

        n_freq = loop_snr.shape[1]
        snr_interp = np.full((n_pts, n_freq), np.nan)

        # Interpolate each freq bin separately across RA
        for fi in range(n_freq):
            # Only interpolate within this loop's RA coverage
            mask = (ra_grid >= loop_ra[0]) & (ra_grid <= loop_ra[-1])
            if mask.sum() < 2:
                continue
            f_interp = interp1d(loop_ra, loop_snr[:, fi], kind='linear',
                                bounds_error=False, fill_value=np.nan)
            snr_interp[mask, fi] = f_interp(ra_grid[mask])

        loop_snr_grids.append(snr_interp)

    return ra_grid, loop_snr_grids


def plot_compare(processed_dirs, y_axis='freq', baseline_correct=True,
                 z_gap=None, alpha=0.65, save_path=None, quality='fast', use_vlsr=False):
    """3D overlay / compare plot: each loop stacked on same RA axis.

    - X axis : RA (right ascension hours) — aligned across loops
    - Y axis : Frequency or velocity
    - Z axis : SNR + loop offset (so surfaces don't sit on top of each other)
    - Each loop gets a distinct colour, slightly transparent
    - NaN gaps where a loop doesn't cover a particular RA

    z_gap: vertical spacing between loops in SNR units.
            None = auto (30% of peak-to-peak SNR range)
    """
    try:
        from scipy.interpolate import interp1d
    except ImportError:
        print("ERROR: scipy required for compare mode.  pip install scipy")
        return

    data, freq, loop_info = load_multiday_data(processed_dirs,
                                                baseline_correct=baseline_correct)
    if not data:
        print("No data to plot!")
        return

    if len(loop_info) < 2:
        print("Need at least 2 loops for compare mode.")
        return

    # Build common RA grid
    quality_pts = {'fast': 200, 'medium': 300, 'high': 500}
    n_pts = quality_pts.get(quality, 200)
    # VLSR: pre-correct each obs's SNR onto common LSR velocity grid before RA interpolation
    if y_axis == 'velocity' and use_vlsr and check_astropy():
        vel_lsr, snr_corrected, vlsr_vals = build_vlsr_grid(data, freq)
        for i, d in enumerate(data):
            d['snr'] = snr_corrected[i]
        freq = vel_lsr   # now "freq" axis is actually LSR velocity
        mean_v = float(np.nanmean(vlsr_vals))
        print(f"  VLSR applied: mean correction {mean_v:+.1f} km/s")
    ra_grid, loop_snr_grids = _build_ra_grid(data, loop_info, n_pts=n_pts)

    # Y axis
    if y_axis == 'velocity':
        y_data  = C_KMS * (F_REST - freq) / F_REST
        y_label = 'Velocity (km/s)'
    else:
        y_data  = freq
        y_label = 'Frequency (MHz)'

    # Auto z_gap: 30% of peak-to-peak range across all valid data
    if z_gap is None:
        all_vals = np.concatenate([g[~np.isnan(g)] for g in loop_snr_grids])
        z_gap = (all_vals.max() - all_vals.min()) * 0.30
        z_gap = max(z_gap, 0.5)   # at least 0.5 dB

    print(f"\nCompare mode: {len(loop_info)} loops, z_gap={z_gap:.2f} dB, RA {ra_grid[0]:.2f}-{ra_grid[-1]:.2f}h")

    fig = plt.figure(figsize=(16, 10))
    ax  = fig.add_subplot(111, projection='3d')

    X, Y = np.meshgrid(ra_grid, y_data, indexing='ij')   # (n_pts, n_freq)

    legend_patches = []
    for li, (info, snr_grid) in enumerate(zip(loop_info, loop_snr_grids)):
        Z = snr_grid + li * z_gap   # vertical offset per loop

        # Build a solid-color facecolor array with transparency
        cmap_name  = LOOP_COLORS_MPL[li % len(LOOP_COLORS_MPL)]
        base_color = plt.get_cmap(cmap_name)(0.65)   # mid-range of colormap
        face_rgba  = np.full((*Z.shape, 4), base_color)
        face_rgba[..., 3] = alpha
        # Dim where data is NaN (make transparent)
        face_rgba[np.isnan(Z), 3] = 0.0
        Z_plot = np.where(np.isnan(Z), np.nanmin(Z[~np.isnan(Z)]) if (~np.isnan(Z)).any() else 0, Z)

        ax.plot_surface(X, Y, Z_plot, facecolors=face_rgba,
                        linewidth=0, antialiased=True,
                        rstride=max(1, n_pts // 80),
                        cstride=max(1, len(y_data) // 80))

        # Build legend patch
        import matplotlib.patches as mpatches
        patch = mpatches.Patch(color=base_color,
                               label=f"{info['date']}  RA {info['start_ra']:.1f}-{info['end_ra']:.1f}h  (+{li * z_gap:.1f} dB offset)")
        legend_patches.append(patch)

        print(f"  Loop {li+1}: {info['date']}  {info['n_obs']} obs  offset={li * z_gap:.1f} dB")

    ax.legend(handles=legend_patches, loc='upper left', fontsize=8)

    ax.set_xlabel('RA (hours)')
    ax.set_ylabel(y_label)
    ax.set_zlabel('SNR + offset (dB)')
    ax.set_title(
        f'H-Line Observatory — Loop Comparison\n'
        f'DEC: {data[0]["dec"]:.1f}°  |  '
        + "  vs  ".join(i["date"] for i in loop_info),
        fontsize=10
    )

    # Add H-line marker
    if y_axis == 'velocity':
        ax.plot([ra_grid[0], ra_grid[-1]], [0, 0],
                [ax.get_zlim()[0]] * 2, 'w--', alpha=0.4, lw=0.8)
    else:
        ax.plot([ra_grid[0], ra_grid[-1]], [F_REST, F_REST],
                [ax.get_zlim()[0]] * 2, 'w--', alpha=0.4, lw=0.8)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {save_path}")
        plt.close()
    else:
        plt.show()


def plot_compare_plotly(processed_dirs, y_axis='freq', baseline_correct=True,
                        z_gap=None, alpha=0.55, save_path=None, use_vlsr=False):
    """Interactive Plotly overlay compare plot: each loop on same RA axis.

    Same concept as plot_compare() but fully interactive HTML output.
    Each loop is a separate Surface trace with its own colourscale.
    """
    if not PLOTLY_AVAILABLE:
        print("ERROR: Plotly not installed!  Run: pip install plotly")
        return

    try:
        from scipy.interpolate import interp1d
    except ImportError:
        print("ERROR: scipy required for compare mode.  pip install scipy")
        return

    data, freq, loop_info = load_multiday_data(processed_dirs,
                                                baseline_correct=baseline_correct)
    if not data:
        print("No data to plot!")
        return

    if len(loop_info) < 2:
        print("Need at least 2 loops for compare mode.")
        return

    ra_grid, loop_snr_grids = _build_ra_grid(data, loop_info, n_pts=350)

    # Y axis
    if y_axis == 'velocity':
        y_data  = C_KMS * (F_REST - freq) / F_REST
        y_label = 'Velocity (km/s)'
    else:
        y_data  = freq
        y_label = 'Frequency (MHz)'

    # Auto z_gap
    if z_gap is None:
        all_vals = np.concatenate([g[~np.isnan(g)] for g in loop_snr_grids])
        z_gap = (all_vals.max() - all_vals.min()) * 0.30
        z_gap = max(z_gap, 0.5)

    print(f"\nCompare mode (Plotly): {len(loop_info)} loops, z_gap={z_gap:.2f} dB")

    X, Y = np.meshgrid(ra_grid, y_data, indexing='ij')

    # Plotly single-color colourscales (name -> [[0, rgba], [1, rgba]])
    def solid_colorscale(hex_col, opacity):
        r, g, b = int(hex_col[1:3], 16), int(hex_col[3:5], 16), int(hex_col[5:7], 16)
        c = f'rgba({r},{g},{b},{opacity})'
        return [[0, c], [1, c]]

    traces = []
    for li, (info, snr_grid) in enumerate(zip(loop_info, loop_snr_grids)):
        Z = snr_grid + li * z_gap
        hex_col = LOOP_COLORS_HEX[li % len(LOOP_COLORS_HEX)]

        hover = (f"<b>{info['date']}</b><br>"
                 "RA: %{x:.2f}h<br>"
                 f"{y_label}: " + "%{y:.3f}<br>"
                 "SNR: %{customdata:.2f} dB<extra></extra>")

        traces.append(go.Surface(
            x=X, y=Y, z=Z,
            colorscale=solid_colorscale(hex_col, alpha),
            showscale=False,
            name=info['date'],
            customdata=snr_grid,   # raw SNR without offset for hover
            hovertemplate=hover,
            opacity=alpha,
            lighting=dict(ambient=0.95, diffuse=0.9, specular=0.05, fresnel=0.0),
        ))

        print(f"  Loop {li+1}: {info['date']}  offset={li * z_gap:.1f} dB  color={hex_col}")

    fig = go.Figure(data=traces)

    # Add H-line reference plane (thin flat surface at SNR=0 across all loops)
    h_x_ref = F_REST if y_axis == 'freq' else 0.0
    fig.add_trace(go.Scatter3d(
        x=[ra_grid[0], ra_grid[-1]],
        y=[h_x_ref, h_x_ref],
        z=[0, 0],
        mode='lines',
        line=dict(color='white', dash='dash', width=2),
        name='H-line rest',
        hoverinfo='skip',
    ))

    date_legend = "  vs  ".join(
        f"<span style='color:{LOOP_COLORS_HEX[i]}'>{info['date']}</span>"
        for i, info in enumerate(loop_info)
    )

    fig.update_layout(
        title=dict(
            text=f'<b>H-Line Observatory</b> — Loop Comparison<br>'
                 f'DEC: {data[0]["dec"]:.1f}°  |  {date_legend}',
            x=0.5, xanchor='center'
        ),
        scene=dict(
            xaxis=dict(title='RA (hours)', gridcolor='rgba(128,128,128,0.3)',
                       showbackground=True, backgroundcolor='rgb(15,15,25)'),
            yaxis=dict(title=y_label, gridcolor='rgba(128,128,128,0.3)',
                       showbackground=True, backgroundcolor='rgb(15,15,25)'),
            zaxis=dict(title='SNR + offset (dB)', gridcolor='rgba(128,128,128,0.3)',
                       showbackground=True, backgroundcolor='rgb(15,15,25)'),
            bgcolor='rgb(10,10,20)',
            camera=dict(eye=dict(x=1.6, y=1.4, z=0.9)),
            aspectmode='manual',
            aspectratio=dict(x=1.2, y=1, z=0.8),
        ),
        paper_bgcolor='rgb(10,10,20)',
        font=dict(color='white'),
        legend=dict(x=0.01, y=0.99, font=dict(color='white')),
        width=1300,
        height=850,
        margin=dict(l=10, r=10, t=100, b=10),
    )

    if save_path is None:
        save_path = os.path.join(LOCAL_OUTPUT, "hline_compare_interactive.html")

    fig.write_html(save_path, include_plotlyjs=True, full_html=True)
    print(f"Saved: {save_path}")
    print("Open in browser — each loop is a separate colour, hover for SNR values!")

# ============== STACK / CO-ADD FUNCTIONS ==============

def stack_loops(processed_dirs, baseline_correct=True, use_vlsr=True,
                ra_bin_width=0.1, min_obs=1):
    """VLSR-correct every observation then co-add by RA bin.

    For a fixed-pointing dish the same sky drifts through the beam at the
    same RA each sidereal day.  Stacking multiple loops at the same RA
    improves SNR by sqrt(N).

    Steps:
      1. Load all loops
      2. VLSR-correct each spectrum onto a common LSR velocity grid
      3. Bin observations by RA (width = ra_bin_width hours)
      4. Median-average spectra within each bin across all loops
    Returns:
        stacked  : list of dicts {ra, dec, snr, n_stacked, dt_utc, dt_ist, hours}
        freq_out : frequency or LSR-velocity array
        meta     : summary dict
    """
    try:
        from scipy.interpolate import interp1d
    except ImportError:
        print("ERROR: scipy required for stacking.  pip install scipy")
        return [], None, {}

    all_data = []
    common_freq = None
    loop_summary = []

    for processed_dir in processed_dirs:
        data, freq = load_processed_data(processed_dir, baseline_correct=baseline_correct)
        if data:
            all_data.extend(data)
            if common_freq is None:
                common_freq = freq
            loop_summary.append({
                'dir': os.path.basename(processed_dir),
                'n': len(data),
                'ra_range': (min(d['ra'] for d in data), max(d['ra'] for d in data)),
            })

    if not all_data:
        print("ERROR: No data loaded!")
        return [], None, {}

    print(f"\nStack: {len(all_data)} total observations from {len(processed_dirs)} loop(s)")
    for li in loop_summary:
        print(f"  {li['dir']}: {li['n']} obs, RA {li['ra_range'][0]:.2f}-{li['ra_range'][1]:.2f}h")

    # VLSR correct onto common LSR velocity grid
    if use_vlsr and check_astropy():
        print("  Applying VLSR correction to each observation...")
        vel_grid, snr_corrected, vlsr_vals = build_vlsr_grid(all_data, common_freq, use_vlsr=True)
        for i, d in enumerate(all_data):
            d['snr_lsr'] = snr_corrected[i]
        freq_out = vel_grid
        axis_label = 'velocity_lsr'
        mean_vlsr = float(np.nanmean(vlsr_vals))
        print(f"  VLSR range: {vlsr_vals.min():.2f} to {vlsr_vals.max():.2f} km/s  "
              f"(mean={mean_vlsr:+.2f} km/s)")
    else:
        for d in all_data:
            d['snr_lsr'] = d['snr']
        freq_out = common_freq
        axis_label = 'freq'
        mean_vlsr = 0.0

    # Bin by RA and median stack
    all_ra = np.array([d['ra'] for d in all_data])
    ra_min, ra_max = all_ra.min(), all_ra.max()
    n_bins = max(1, int(round((ra_max - ra_min) / ra_bin_width)))
    bin_edges = np.linspace(ra_min, ra_max, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    stacked = []
    n_skipped = 0
    for bi, ra_center in enumerate(bin_centers):
        lo, hi = bin_edges[bi], bin_edges[bi + 1]
        # Include upper edge for last bin so ra_max obs aren't dropped
        is_last = (bi == len(bin_centers) - 1)
        members = [d for d in all_data if lo <= d['ra'] <= hi] if is_last \
                  else [d for d in all_data if lo <= d['ra'] < hi]
        if len(members) < min_obs:
            n_skipped += 1
            continue
        stack_snr = np.nanmedian(np.array([d['snr_lsr'] for d in members]), axis=0)
        members_sorted = sorted(members, key=lambda d: d['dt_utc'])
        rep = members_sorted[len(members_sorted) // 2]
        stacked.append({
            'ra':        ra_center,
            'dec':       float(np.mean([d['dec'] for d in members])),
            'snr':       stack_snr,
            'n_stacked': len(members),
            'dt_utc':    rep['dt_utc'],
            'dt_ist':    rep['dt_ist'],
            'hours':     float(bi),
            'loop_idx':  0,
        })

    if not stacked:
        print("ERROR: No bins had enough observations!")
        return [], None, {}

    counts = [s['n_stacked'] for s in stacked]
    print(f"  Stacked into {len(stacked)} RA bins (skipped {n_skipped} empty bins)")
    print(f"  Obs per bin: min={min(counts)}  max={max(counts)}  mean={np.mean(counts):.1f}")
    print(f"  Expected SNR improvement: ~{np.sqrt(np.mean(counts)):.1f}x")

    # Collect dates from actual loop data (not from stacked bins which may miss loops)
    loop_dates = []
    seen_dates = set()
    for processed_dir in processed_dirs:
        loop_data = [d for d in all_data if os.path.dirname(d['file']) == processed_dir]
        if loop_data:
            d = loop_data[0]['dt_ist'].strftime('%Y-%m-%d')
            if d not in seen_dates:
                loop_dates.append(d)
                seen_dates.add(d)

    meta = {
        'axis_label':   axis_label,
        'mean_vlsr':    mean_vlsr,
        'ra_bin_width': ra_bin_width,
        'n_loops':      len(processed_dirs),
        'n_total_obs':  len(all_data),
        'n_bins':       len(stacked),
        'loop_dates':   loop_dates,
    }
    return stacked, freq_out, meta


def plot_stack(processed_dirs, y_axis='velocity', baseline_correct=True,
               use_vlsr=True, ra_bin_width=0.1, save_path=None, format='show', view='3d'):
    """Plot the VLSR-corrected RA-stacked surface.

    X = RA (hours), Y = velocity or frequency, Z = median-stacked SNR.
    Hover (Plotly) shows how many observations were stacked per bin.
    """
    stacked, freq_out, meta = stack_loops(
        processed_dirs, baseline_correct=baseline_correct,
        use_vlsr=use_vlsr, ra_bin_width=ra_bin_width
    )
    if not stacked:
        return

    ra_arr  = np.array([s['ra']        for s in stacked])
    snr_arr = np.array([s['snr']       for s in stacked])
    cnt_arr = np.array([s['n_stacked'] for s in stacked])

    # Y axis
    if y_axis == 'velocity' and meta['axis_label'] == 'velocity_lsr':
        y_data  = freq_out
        y_label = f"Velocity LSR (km/s)  [mean VLSR={meta['mean_vlsr']:+.1f} km/s]"
    elif y_axis == 'velocity':
        y_data  = C_KMS * (F_REST - freq_out) / F_REST
        y_label = 'Velocity topocentric (km/s)'
    else:
        y_data  = freq_out
        y_label = 'Frequency (MHz)'

    # Title — use loop_dates from meta (guaranteed to have all loops)
    date_tags = meta.get('loop_dates', [])
    date_str = ' + '.join(date_tags[:4])
    if len(date_tags) > 4:
        date_str += f' (+{len(date_tags)-4} more)'
    title = (f'H-Line Observatory — Stacked ({meta["n_loops"]} loops, '
             f'{meta["n_total_obs"]} obs → {meta["n_bins"]} RA bins)\n'
             f'{date_str}  |  DEC: {stacked[0]["dec"]:.1f}°  |  '
             f'bin width={meta["ra_bin_width"]:.2f}h')

    # ── Plotly ─────────────────────────────────────────────────────────────
    if format == 'html':
        if not PLOTLY_AVAILABLE:
            print("ERROR: Plotly not installed.  pip install plotly"); return
        X, Y = np.meshgrid(ra_arr, y_data, indexing='ij')
        cnt_grid = np.broadcast_to(cnt_arr[:, np.newaxis], snr_arr.shape).copy()
        fig = go.Figure(data=[go.Surface(
            x=X, y=Y, z=snr_arr,
            colorscale='Viridis',
            colorbar=dict(title=dict(text='SNR (dB)', side='right')),
            customdata=cnt_grid,
            hovertemplate=(
                'RA: %{x:.2f}h<br>'
                'Y: %{y:.2f}<br>'
                'SNR: %{z:.3f} dB<br>'
                'Stacked: %{customdata} obs<extra></extra>'
            ),
            lighting=dict(ambient=0.9, diffuse=0.9, specular=0.1),
        )])
        fig.update_layout(
            title=dict(text=title.replace('\n', '<br>'), x=0.5, xanchor='center'),
            scene=dict(
                xaxis=dict(title='RA (hours)', gridcolor='rgba(128,128,128,0.3)',
                           showbackground=True, backgroundcolor='rgb(15,15,25)'),
                yaxis=dict(title=y_label, gridcolor='rgba(128,128,128,0.3)',
                           showbackground=True, backgroundcolor='rgb(15,15,25)'),
                zaxis=dict(title='SNR (dB)', gridcolor='rgba(128,128,128,0.3)',
                           showbackground=True, backgroundcolor='rgb(15,15,25)'),
                bgcolor='rgb(10,10,20)',
                camera=dict(eye=dict(x=1.6, y=1.4, z=0.9)),
                aspectmode='manual', aspectratio=dict(x=1.2, y=1, z=0.7),
            ),
            paper_bgcolor='rgb(10,10,20)', font=dict(color='white'),
            width=1300, height=850, margin=dict(l=10, r=10, t=110, b=10),
        )
        if save_path is None:
            save_path = os.path.join(LOCAL_OUTPUT, 'hline_stacked.html')
        fig.write_html(save_path, include_plotlyjs=True, full_html=True)
        print(f"Saved: {save_path}")
        print("Hover shows n_stacked per RA bin — coverage map built in!")
        return

    # ── Matplotlib 2D heatmap ───────────────────────────────────────────────
    if view == 'top':
        fig, ax = plt.subplots(figsize=(16, 8))
        X, Y = np.meshgrid(ra_arr, y_data, indexing='ij')
        mesh = ax.pcolormesh(ra_arr, y_data, snr_arr.T, cmap='viridis', shading='auto')
        cbar = fig.colorbar(mesh, ax=ax, label='SNR (dB)')
        ax.set_xlabel('RA (hours)')
        ax.set_ylabel(y_label)
        ax.set_title(title, fontsize=9)
        if y_axis == 'velocity':
            ax.axhline(y=0, color='white', linestyle='--', alpha=0.5, linewidth=1,
                       label='0 km/s (H-line rest)')
        else:
            ax.axhline(y=F_REST, color='white', linestyle='--', alpha=0.5, linewidth=1,
                       label=f'{F_REST:.3f} MHz')
        ax.legend(loc='upper right', fontsize=8)
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved: {save_path}")
            plt.close()
        else:
            plt.show()
        return

    # ── Matplotlib 3D surface ───────────────────────────────────────────────
    fig = plt.figure(figsize=(16, 10))
    ax  = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(ra_arr, y_data, indexing='ij')
    surf = ax.plot_surface(X, Y, snr_arr, cmap='viridis', alpha=0.9,
                           linewidth=0, antialiased=False,
                           rstride=max(1, len(ra_arr)//80),
                           cstride=max(1, len(y_data)//80))
    fig.colorbar(surf, ax=ax, shrink=0.5, label='SNR (dB)')
    ax.set_xlabel('RA (hours)')
    ax.set_ylabel(y_label)
    ax.set_zlabel('SNR (dB)')
    ax.set_title(title, fontsize=9)
    if y_axis == 'velocity':
        ax.plot([ra_arr[0], ra_arr[-1]], [0, 0],
                [snr_arr.min()]*2, 'w--', alpha=0.5, lw=1, label='0 km/s (H-line rest)')
    else:
        ax.plot([ra_arr[0], ra_arr[-1]], [F_REST, F_REST],
                [snr_arr.min()]*2, 'w--', alpha=0.5, lw=1, label=f'{F_REST:.3f} MHz')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {save_path}")
        plt.close()
    else:
        plt.show()




# ============== EZRA FUNCTIONS ==============
def run_ezra_pipeline(processed_dir, az=DISH_AZ, el=DISH_EL):
    """Convert to ezRA format (single directory)"""
    return run_ezra_pipeline_multi([processed_dir], az, el)

def run_ezra_pipeline_multi(processed_dirs, az=DISH_AZ, el=DISH_EL):
    """Convert multiple processed dirs to a single ezRA .txt file (virgo2ezra style).

    Bypasses ezColDSPIRA entirely — that tool ignores our AZ/EL from filenames
    (reads defaults instead) and then silently writes nothing to ezhr.
    Instead we build one combined .txt file and pass it directly to ezCon.py.

    Output format per line: YYYY-MM-DDTHH:MM:SS <p0> <p1> ... <pN>
    Header: lat, lon, amsl, freqMin, freqMax, freqBinQty, azDeg, elDeg
    """
    ezra_home     = os.path.expanduser("~/ezRA/ezRA")
    ezra_data_dir = os.path.join(ezra_home, "data")
    os.makedirs(ezra_data_dir, exist_ok=True)

    # Gather all processed CSVs
    all_files = []
    for processed_dir in processed_dirs:
        files = sorted(glob.glob(f"{processed_dir}/*_spectra_filtered.csv"))
        if files:
            all_files.extend(files)
            print(f"Found {len(files)} files in {os.path.basename(processed_dir)}")

    if not all_files:
        print("ERROR: No processed CSV files found!")
        return None

    print(f"Total: {len(all_files)} observations")

    # Read first file to build frequency grid (col 0 = freq, col 4 = filtered SNR)
    freqs = []
    with open(all_files[0], 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) >= 5:
                freq = float(row[0])
                if FREQ_MIN <= freq <= FREQ_MAX:
                    freqs.append(freq)
    if not freqs:
        print("ERROR: No frequency data found in range!")
        return None

    freq_min  = freqs[0]
    freq_max  = freqs[-1]
    freq_bins = len(freqs)

    # Wipe ALL files from data/ — old CSVs (from ezColDSPIRA runs) and old txt files
    # would confuse ezCon into mixing datasets
    old_files = [f for f in glob.glob(f"{ezra_data_dir}/*") if os.path.isfile(f)]
    for f in old_files:
        os.remove(f)
    if old_files:
        print(f"Cleared {len(old_files)} old file(s) from data/")

    # Build output filename from date range
    dt0  = ist_to_utc(datetime_from_filename(all_files[0]))
    dt1  = ist_to_utc(datetime_from_filename(all_files[-1]))
    name = f"HLine_{dt0.strftime('%Y%m%d')}_{dt1.strftime('%Y%m%d')}"
    out_path = os.path.join(ezra_data_dir, f"{name}.txt")

    written = 0
    with open(out_path, 'w') as out:
        out.write("from ezColHLinePipeline.py combined\n")
        out.write(f"lat {OBS_LAT} long {OBS_LON} amsl {OBS_ALT} name {name}\n")
        out.write(f"freqMin {freq_min:.4f} freqMax {freq_max:.4f} freqBinQty {freq_bins}\n")
        out.write(f"azDeg {az} elDeg {el}\n")
        out.write("# times are in UTC\n")

        for i, csv_file in enumerate(all_files):
            dt_utc = ist_to_utc(datetime_from_filename(csv_file))
            ts = dt_utc.strftime('%Y-%m-%dT%H:%M:%S')

            powers = []
            with open(csv_file, 'r') as f:
                reader = csv.reader(f)
                for row in reader:
                    if len(row) >= 5:
                        freq = float(row[0])
                        if FREQ_MIN <= freq <= FREQ_MAX:
                            powers.append(float(row[4]))  # col 5 = filtered SNR

            if powers:
                out.write(f"{ts} " + " ".join(f"{p:.5f}" for p in powers) + "\n")
                written += 1

            if (i + 1) % 100 == 0:
                print(f"  Processed {i + 1}/{len(all_files)}...")

    print(f"\nCreated: {out_path}")
    print(f"  {written} observations, {freq_bins} freq bins ({freq_min:.3f}–{freq_max:.3f} MHz)")

    # Run ezCon directly on the .txt file (no ezColDSPIRA needed)
    rel_path = os.path.relpath(out_path, ezra_home)
    print(f"\n=== Running ezCon ===")
    con_cmd = ["python3", "ezCon.py", rel_path]
    print("  " + " ".join(con_cmd))
    result = subprocess.run(con_cmd, cwd=ezra_home)
    if result.returncode != 0:
        print("ERROR: ezCon failed! Check output above.")
        return None

    print(f"\n=== ezRA complete! ===")
    return out_path

# ============== MAIN ==============
def main():
    global OBS_LAT, OBS_LON, OBS_ALT
    parser = argparse.ArgumentParser(description='H-Line Pipeline for H-Line Observatory')
    
    parser.add_argument('--data-dir', default=os.path.expanduser('~/hydrogen_obs'), help='Path to directory containing loop_* folders (default: ~/hydrogen_obs)')
    parser.add_argument('--list', action='store_true', help='List available local loops')
    parser.add_argument('--process', metavar='LOOP', nargs='?', const='latest', help='Reprocess fetched data')
    parser.add_argument('--lat', type=float, default=OBS_LAT, help='Observatory latitude (default: 0.0)')
    parser.add_argument('--lon', type=float, default=OBS_LON, help='Observatory longitude (default: 0.0)')
    parser.add_argument('--altitude', type=float, default=OBS_ALT, help='Observatory altitude in meters (default: 0)')
    parser.add_argument('--plot3d', metavar='LOOP', nargs='*', help='Generate 3D plot (multiple loops = combined)')
    parser.add_argument('--heatmap', metavar='LOOP', nargs='*', help='Generate 2D heatmap (cleaner than 3D)')
    parser.add_argument('--flipbook', metavar='LOOP', nargs='*', help='Generate flipbook (multiple loops = combined)')
    parser.add_argument('--single', metavar='HOUR', type=float, help='Plot single observation at hour N (use with --hours to see available)')
    parser.add_argument('--loop', metavar='LOOP', help='Specify loop for --single and --hours (default: latest)')
    parser.add_argument('--hours', metavar='LOOP', nargs='?', const='latest', help='List available hours')
    parser.add_argument('--ezra', metavar='LOOP', nargs='*', help='Run ezRA pipeline (multiple loops = combined)')
    parser.add_argument('--compare', metavar='LOOP', nargs='+', help='Overlay loops on same RA axis, each a different color (min 2 loops)')
    parser.add_argument('--z-gap', type=float, default=None, help='Vertical SNR offset between loops in compare mode (default: auto)')
    parser.add_argument('--all', metavar='LOOP', nargs='+', help='Fetch, process, and plot (accepts multiple loops)')
    
    # SETI options
    parser.add_argument('--seti', metavar='LOOP', nargs='?', const='latest', help='SETI analysis - waterfall plot and candidate search')
    parser.add_argument('--seti-data', choices=['raw', 'calibrated'], default='calibrated', 
                        help='SETI data type: raw (no processing) or calibrated (VIRGO processed, no median)')
    parser.add_argument('--seti-sigma', type=float, default=5, help='SETI detection threshold in sigma (default: 5)')
    parser.add_argument('--seti-candidates', metavar='LOOP', nargs='?', const='latest', help='Find SETI narrowband candidates')
    
    parser.add_argument('--x-axis', choices=['ra', 'time'], default='ra', help='X-axis for 3D plot')
    parser.add_argument('--y-axis', choices=['freq', 'velocity'], default='freq', help='Y-axis: frequency or velocity')
    parser.add_argument('--render', choices=['grid', 'solid'], default='grid', help='3D render style: grid (fast) or solid (smooth)')
    parser.add_argument('--quality', choices=['fast', 'medium', 'high'], default='fast', help='Matplotlib quality: fast (smooth rotation), medium, high (detailed but slow)')
    parser.add_argument('--view', choices=['3d', 'top', '2d'], default='3d', help='View angle: 3d (default), top/2d (flat heatmap)')
    parser.add_argument('--clean', action='store_true', help='Remove background panes and grids, keep only axes')
    parser.add_argument('--format', choices=['show', 'png', 'html'], default='show', help='Output: show (matplotlib), png (image), html (plotly smooth 3D)')
    parser.add_argument('--no-baseline', action='store_true', help='Disable baseline correction')
    parser.add_argument('--vlsr', action='store_true', help='Apply VLSR correction to velocity axis (requires astropy)')
    parser.add_argument('--stack', metavar='LOOP', nargs='+', help='VLSR-correct then co-add all loops by RA bin (best SNR, fixed dish)')
    parser.add_argument('--ra-bin', type=float, default=0.1, help='RA bin width in hours for --stack (default: 0.1h ~ 6 min)')
    parser.add_argument('--az', type=float, default=DISH_AZ, help=f'Dish azimuth (default: {DISH_AZ})')
    parser.add_argument('--el', type=float, default=DISH_EL, help=f'Dish elevation (default: {DISH_EL})')
    parser.add_argument('--cal', default=CAL_FILE, help='Calibration file path')
    parser.add_argument('--save', metavar='PATH', help='Save plot to file instead of showing')
    parser.add_argument('--export', metavar='LOOP', nargs='+',
                        help='Export loops to Zenodo-ready zip (CSV + plots).')
    parser.add_argument('--name', metavar='NAME', default=None,
                        help='Zip filename for --export (without .zip)')
    parser.add_argument('--dest', metavar='DIR', default=os.path.expanduser('~/Downloads'),
                        help='Destination folder for --export zip (default: ~/Downloads)')
    
    args = parser.parse_args()

    OBS_LAT = args.lat
    OBS_LON = args.lon
    OBS_ALT = args.altitude
    data_dir = os.path.expanduser(args.data_dir)

    # Create directories
    for d in [LOCAL_BASE, LOCAL_RAW, LOCAL_PROCESSED, LOCAL_OUTPUT, LOCAL_CAL]:
        os.makedirs(d, exist_ok=True)

    # Helper to find latest loop
    def find_loop(loop_arg, in_dir):
        if loop_arg == 'latest':
            loops = sorted(glob.glob(f"{in_dir}/loop_*"))
            if loops:
                return loops[-1]
            else:
                print(f"No loops found in {in_dir}")
                return None
        else:
            if not loop_arg.startswith("loop_"):
                # Try to find matching loop
                matches = glob.glob(f"{in_dir}/loop_{loop_arg}*")
                if matches:
                    return matches[0]
            path = os.path.join(in_dir, loop_arg)
            if os.path.exists(path):
                return path
        print(f"Loop not found: {loop_arg}")
        return None

    # Execute commands
    if args.list:
        list_local_loops(data_dir)

    elif args.process:
        loop_path = find_loop(args.process, LOCAL_RAW)
        if not loop_path:
            loop_path = find_loop(args.process, data_dir)
        if loop_path:
            reprocess_loop(loop_path, args.cal)
    
    elif args.plot3d is not None:
        baseline = not args.no_baseline
        
        # Determine save path based on format
        if args.format == 'png':
            save_path = args.save if args.save else os.path.join(LOCAL_OUTPUT, "hline_3d.png")
        elif args.format == 'html':
            save_path = args.save if args.save else os.path.join(LOCAL_OUTPUT, "hline_3d_interactive.html")
        else:
            save_path = args.save  # None = show interactive
        
        if len(args.plot3d) == 0:
            processed_path = find_loop('latest', LOCAL_PROCESSED)
        elif len(args.plot3d) == 1:
            processed_path = find_loop(args.plot3d[0], LOCAL_PROCESSED)
        else:
            # Multiple loops = combined
            processed_paths = []
            for loop in args.plot3d:
                p = find_loop(loop, LOCAL_PROCESSED)
                if p:
                    processed_paths.append(p)
            if processed_paths:
                if args.format == 'html':
                    plot_3d_plotly_multiday(processed_paths, y_axis=args.y_axis, baseline_correct=baseline, save_path=save_path, use_vlsr=args.vlsr)
                else:
                    plot_3d_multiday(processed_paths, y_axis=args.y_axis, baseline_correct=baseline, render=args.render, view=args.view, clean=args.clean, quality=args.quality, save_path=save_path, use_vlsr=args.vlsr)
            processed_path = None
        
        if processed_path:
            if args.format == 'html':
                plot_3d_plotly(processed_path, x_axis=args.x_axis, y_axis=args.y_axis, baseline_correct=baseline, save_path=save_path, use_vlsr=args.vlsr)
            else:
                plot_3d(processed_path, x_axis=args.x_axis, y_axis=args.y_axis, baseline_correct=baseline, render=args.render, view=args.view, clean=args.clean, quality=args.quality, save_path=save_path, use_vlsr=args.vlsr)
    
    elif args.heatmap is not None:
        baseline = not args.no_baseline
        if len(args.heatmap) == 0:
            # No args = latest
            processed_path = find_loop('latest', LOCAL_PROCESSED)
            if processed_path:
                plot_heatmap(processed_path, x_axis=args.x_axis, y_axis=args.y_axis, baseline_correct=baseline, save_path=args.save, use_vlsr=args.vlsr)
        elif len(args.heatmap) == 1:
            # Single loop
            processed_path = find_loop(args.heatmap[0], LOCAL_PROCESSED)
            if processed_path:
                plot_heatmap(processed_path, x_axis=args.x_axis, y_axis=args.y_axis, baseline_correct=baseline, save_path=args.save, use_vlsr=args.vlsr)
        else:
            # Multiple loops - combine data
            print("Heatmap multi-loop: combining data...")
            # For now just use first loop
            processed_path = find_loop(args.heatmap[0], LOCAL_PROCESSED)
            if processed_path:
                plot_heatmap(processed_path, x_axis=args.x_axis, y_axis=args.y_axis, baseline_correct=baseline, save_path=args.save, use_vlsr=args.vlsr)
    
    elif args.flipbook is not None:
        baseline = not args.no_baseline
        if len(args.flipbook) == 0:
            # No args = latest
            processed_path = find_loop('latest', LOCAL_PROCESSED)
            if processed_path:
                plot_flipbook(processed_path, y_axis=args.y_axis, baseline_correct=baseline, save_path=args.save)
        elif len(args.flipbook) == 1:
            # Single loop
            processed_path = find_loop(args.flipbook[0], LOCAL_PROCESSED)
            if processed_path:
                plot_flipbook(processed_path, y_axis=args.y_axis, baseline_correct=baseline, save_path=args.save)
        else:
            # Multiple loops = combined
            processed_paths = []
            for loop in args.flipbook:
                p = find_loop(loop, LOCAL_PROCESSED)
                if p:
                    processed_paths.append(p)
            if processed_paths:
                plot_flipbook_multiday(processed_paths, y_axis=args.y_axis, baseline_correct=baseline, save_path=args.save)
    
    elif args.single is not None:
        loop_arg = args.loop if args.loop else 'latest'
        processed_path = find_loop(loop_arg, LOCAL_PROCESSED)
        if processed_path:
            plot_single(processed_path, args.single)
        else:
            print("Use --hours to list available observations first")
    
    elif args.hours:
        processed_path = find_loop(args.hours, LOCAL_PROCESSED)
        if processed_path:
            list_hours(processed_path)

    elif args.compare:
        if len(args.compare) < 2:
            print("ERROR: --compare needs at least 2 loops.")
        else:
            processed_paths = []
            for loop in args.compare:
                p = find_loop(loop, LOCAL_PROCESSED)
                if p:
                    processed_paths.append(p)
            if len(processed_paths) < 2:
                print("ERROR: Could not find at least 2 valid loops.")
            else:
                baseline = not args.no_baseline
                if args.format == 'html':
                    save_path = args.save if args.save else os.path.join(LOCAL_OUTPUT, "hline_compare_interactive.html")
                    plot_compare_plotly(processed_paths, y_axis=args.y_axis,
                                        baseline_correct=baseline,
                                        z_gap=args.z_gap, save_path=save_path, use_vlsr=args.vlsr)
                else:
                    save_path = args.save if args.save else (
                        os.path.join(LOCAL_OUTPUT, "hline_compare.png") if args.format == 'png' else None
                    )
                    plot_compare(processed_paths, y_axis=args.y_axis,
                                 baseline_correct=baseline,
                                 z_gap=args.z_gap, quality=args.quality,
                                 save_path=save_path, use_vlsr=args.vlsr)

    elif args.stack:
        processed_paths = []
        for loop in args.stack:
            p = find_loop(loop, LOCAL_PROCESSED)
            if p:
                processed_paths.append(p)
        if not processed_paths:
            print("ERROR: No valid loops found.")
        else:
            baseline = not args.no_baseline
            y_axis = 'velocity' if args.vlsr else args.y_axis
            if args.format == 'html':
                save_path = args.save if args.save else os.path.join(LOCAL_OUTPUT, "hline_stacked.html")
                plot_stack(processed_paths, y_axis=y_axis,
                           baseline_correct=baseline, use_vlsr=args.vlsr,
                           ra_bin_width=args.ra_bin,
                           save_path=save_path, format='html')
            else:
                save_path = args.save if args.save else (
                    os.path.join(LOCAL_OUTPUT, "hline_stacked.png") if args.format == 'png' else None
                )
                plot_stack(processed_paths, y_axis=y_axis,
                           baseline_correct=baseline, use_vlsr=args.vlsr,
                           ra_bin_width=args.ra_bin,
                           save_path=save_path, format=args.format,
                           view=args.view)

    elif args.ezra is not None:
        if len(args.ezra) == 0:
            # No args = latest
            processed_path = find_loop('latest', LOCAL_PROCESSED)
            if processed_path:
                result = run_ezra_pipeline(processed_path, args.az, args.el)
                if not result:
                    print("ezRA conversion failed!")
        elif len(args.ezra) == 1:
            # Single loop
            processed_path = find_loop(args.ezra[0], LOCAL_PROCESSED)
            if processed_path:
                result = run_ezra_pipeline(processed_path, args.az, args.el)
                if not result:
                    print("ezRA conversion failed!")
        else:
            # Multiple loops = combined
            processed_paths = []
            for loop in args.ezra:
                p = find_loop(loop, LOCAL_PROCESSED)
                if p:
                    processed_paths.append(p)
            if processed_paths:
                result = run_ezra_pipeline_multi(processed_paths, args.az, args.el)
                if not result:
                    print("ezRA conversion failed!")
    
    elif args.all:
        print("=== FULL PIPELINE ===\n")

        all_processed = []
        baseline = not args.no_baseline

        for loop_name in args.all:
            print(f"\n--- Processing: {loop_name} ---")

            # Find loop locally
            print("Step 1: Locating loop data...")
            local_raw = find_loop(loop_name, data_dir)
            if not local_raw:
                local_raw = find_loop(loop_name, LOCAL_RAW)
            if not local_raw:
                print(f"Skipping {loop_name} - not found")
                continue

            # Process
            print("\nStep 2: Reprocessing with Virgo...")
            processed_path = reprocess_loop(local_raw, args.cal)
            if not processed_path:
                print(f"Skipping {loop_name} - reprocessing failed")
                continue

            all_processed.append(processed_path)

        if not all_processed:
            print("\nERROR: No loops processed successfully!")
            return

        # Plot combined data
        print(f"\n--- Generating combined plots for {len(all_processed)} loops ---")

        print("\nStep 3: Generating 3D plot...")
        if args.format == 'html':
            save_3d = os.path.join(LOCAL_OUTPUT, "hline_3d_combined.html")
            if len(all_processed) == 1:
                plot_3d_plotly(all_processed[0], x_axis=args.x_axis, y_axis=args.y_axis, baseline_correct=baseline, save_path=save_3d)
            else:
                plot_3d_plotly_multiday(all_processed, y_axis=args.y_axis, baseline_correct=baseline, save_path=save_3d)
        else:
            save_3d = os.path.join(LOCAL_OUTPUT, "hline_3d_combined.png")
            if len(all_processed) == 1:
                plot_3d(all_processed[0], x_axis=args.x_axis, y_axis=args.y_axis, baseline_correct=baseline, render=args.render, view=args.view, clean=args.clean, quality=args.quality, save_path=save_3d)
            else:
                plot_3d_multiday(all_processed, y_axis=args.y_axis, baseline_correct=baseline, render=args.render, view=args.view, clean=args.clean, quality=args.quality, save_path=save_3d)

        print("\nStep 4: Generating flipbook...")
        save_flip = os.path.join(LOCAL_OUTPUT, "hline_flipbook_combined.gif")
        if len(all_processed) == 1:
            plot_flipbook(all_processed[0], y_axis=args.y_axis, baseline_correct=baseline, save_path=save_flip)
        else:
            plot_flipbook_multiday(all_processed, y_axis=args.y_axis, baseline_correct=baseline, save_path=save_flip)

        print("\nStep 5: Converting to ezRA format...")
        run_ezra_pipeline_multi(all_processed, args.az, args.el)

        print("\n=== DONE ===")
        print(f"Processed loops: {all_processed}")
        print(f"Outputs: {LOCAL_OUTPUT}")
    
    elif args.seti is not None:
        # SETI waterfall plot
        if args.seti_data == 'raw':
            # Raw data needs the fetched loop (not processed)
            loop_path = find_loop(args.seti, LOCAL_RAW)
            if not loop_path:
                print(f"ERROR: Loop {args.seti} not found in {LOCAL_RAW}")
                print("Place loop data in your --data-dir first")
                return
        else:
            # Calibrated data needs processed loop
            loop_path = find_loop(args.seti, LOCAL_PROCESSED)
            if not loop_path:
                print(f"ERROR: Loop {args.seti} not found in {LOCAL_PROCESSED}")
                print("Process it first with: --process LOOP")
                return
        
        save_path = args.save  # None = show interactive
        plot_seti_waterfall(loop_path, data_type=args.seti_data, view=args.view, format=args.format, x_axis=args.x_axis, save_path=save_path)
    
    elif args.seti_candidates is not None:
        # SETI candidate search
        if args.seti_data == 'raw':
            loop_path = find_loop(args.seti_candidates, LOCAL_RAW)
            if not loop_path:
                print(f"ERROR: Loop {args.seti_candidates} not found in {LOCAL_RAW}")
                return
        else:
            loop_path = find_loop(args.seti_candidates, LOCAL_PROCESSED)
            if not loop_path:
                print(f"ERROR: Loop {args.seti_candidates} not found in {LOCAL_PROCESSED}")
                return
        
        candidates = seti_find_candidates(loop_path, data_type=args.seti_data, sigma=args.seti_sigma)
        seti_report(candidates)
    
    elif args.export:
        import zipfile
        from datetime import datetime

        processed_paths = []
        for loop in args.export:
            p = find_loop(loop, LOCAL_PROCESSED)
            if p:
                processed_paths.append(p)
            else:
                print(f"WARNING: Loop not found: {loop}")

        if not processed_paths:
            print("ERROR: No valid loops found for export.")
        else:
            # Auto-generate zip name from date range
            if args.name:
                zip_name = args.name
            else:
                dates = []
                for p in processed_paths:
                    d, _ = load_processed_data(p, baseline_correct=False)
                    if d:
                        dates.append(d[0]['dt_ist'].strftime('%Y%m%d'))
                if len(dates) >= 2:
                    zip_name = f"hline_21cm_{dates[0]}_to_{dates[-1]}_{len(processed_paths)}loops"
                elif dates:
                    zip_name = f"hline_21cm_{dates[0]}_{len(processed_paths)}loops"
                else:
                    zip_name = f"hline_21cm_export_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

            os.makedirs(args.dest, exist_ok=True)
            zip_path = os.path.join(args.dest, zip_name + '.zip')

            print(f"\nExporting {len(processed_paths)} loops -> {zip_path}\n")

            with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zf:

                # CSV files
                total_csv = 0
                for p in processed_paths:
                    loop_name = os.path.basename(p)
                    csvs = sorted(glob.glob(f"{p}/*_spectra_filtered.csv"))
                    for csv_path in csvs:
                        arcname = os.path.join('csv', loop_name, os.path.basename(csv_path))
                        zf.write(csv_path, arcname)
                        total_csv += 1
                print(f"  [1/3] Added {total_csv} CSV files")

                # Generate + add VLSR stacked HTML
                print("  [2/3] Generating VLSR stacked plot...")
                stack_html = os.path.join(LOCAL_OUTPUT, f"{zip_name}_stacked.html")
                plot_stack(processed_paths, y_axis='velocity',
                           baseline_correct=True, use_vlsr=True,
                           ra_bin_width=0.1, save_path=stack_html, format='html')
                if os.path.exists(stack_html):
                    zf.write(stack_html, f"plots/{zip_name}_stacked_vlsr.html")

                # Generate + add heatmap PNG
                print("  [3/3] Generating multiday heatmap...")
                heatmap_png = os.path.join(LOCAL_OUTPUT, f"{zip_name}_heatmap.png")
                plot_3d_multiday(processed_paths, y_axis='freq',
                                 baseline_correct=True, view='top',
                                 save_path=heatmap_png)
                if os.path.exists(heatmap_png):
                    zf.write(heatmap_png, f"plots/{zip_name}_heatmap.png")

                # README
                date_list = ', '.join(os.path.basename(p) for p in processed_paths)
                readme = (
                    "H-Line Observatory - Hydrogen Line Export\n"
                    "===============================================\n"
                    f"Generated : {datetime.now().strftime('%Y-%m-%d %H:%M')}\n"
                    f"Loops     : {len(processed_paths)}\n"
                    f"Datasets  : {date_list}\n\n"
                    "Contents:\n"
                    "  csv/    - Per-loop spectra_filtered.csv (frequency vs SNR)\n"
                    "  plots/  - VLSR stacked interactive HTML + multiday heatmap PNG\n\n"
                    "Observatory : H-Line Observatory\n"
                    f"Dish        : prime focus\n"
                    f"Frequency   : {FREQ_MIN}-{FREQ_MAX} MHz (H-line 1420.406 MHz)\n"
                )
                zf.writestr('README.txt', readme)

            size_mb = os.path.getsize(zip_path) / 1e6
            print(f"\nDone!  {zip_path}  ({size_mb:.1f} MB)")
            print("Ready to upload to Zenodo!")

    else:
        parser.print_help()

if __name__ == "__main__":

    main()
