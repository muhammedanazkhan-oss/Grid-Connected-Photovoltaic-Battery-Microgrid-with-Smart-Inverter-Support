#!/usr/bin/env python3
"""Reproduce the Python verification results for the PV-battery microgrid study.

Runs the analysis pipeline in order. All result files (JSON/NPZ) are written to
this directory. Requires only numpy and scipy (see ../../requirements.txt).
"""
import os, sys, subprocess, time

HERE = os.path.dirname(os.path.abspath(__file__))
STEPS = [
    ("microgrid_core.py", "single-diode PV + 24 h dispatch  -> results_core.json, ref_series.npz"),
    ("analyses.py",       "energy tables + sensitivity       -> results_ext.json, cloudy_series.npz"),
    ("econ_tariff.py",    "SAR techno-economics              -> updates results_ext.json"),
    ("thd_sim.py",        "time-domain LCL THD (FFT)         -> results_thd.json, thd_wave.npz"),
    ("month_analysis.py", "Riyadh 12-month analysis          -> results_riyadh.json, riyadh_*.npz"),
]

def main():
    os.chdir(HERE)
    for i, (script, desc) in enumerate(STEPS, 1):
        print(f"[{i}/{len(STEPS)}] {script:<18} {desc}", flush=True)
        t0 = time.time()
        r = subprocess.run([sys.executable, script])
        if r.returncode != 0:
            sys.exit(f"  ERROR: {script} failed")
        print(f"      done in {time.time()-t0:.1f} s", flush=True)
    print("\nAll Python results reproduced in:", HERE)

if __name__ == "__main__":
    main()
