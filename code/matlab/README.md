# PV–Battery Microgrid — MATLAB / Simulink model

Runnable MATLAB and Simulink model for the study *"Theoretical Formulation and
Simulation-Based Verification of a Grid-Connected Photovoltaic–Battery Microgrid
with Smart-Inverter Support for High-Irradiance Residential Applications in Saudi
Arabia"* (manuscript solar-4407776).

Two parts:

1. **MATLAB script suite** — reproduces the study numerically: single-diode PV,
   24 h EMS dispatch, time-domain LCL THD, twelve-month seasonal analysis,
   battery-size sensitivity and techno-economics. Runs in **base MATLAB R2024b**
   (and in **GNU Octave**). **No toolboxes required.**
2. **Simulink model** (`microgrid_lcl.slx`) — a switching-level single-phase
   inverter + LCL filter that yields the grid-current THD by FFT. Built from
   **core Simulink blocks only** (no Simscape / SimPowerSystems).

> Figure-plotting scripts are intentionally excluded from this archive.

## Requirements
- MATLAB R2024b (R2016b+ should also work) for the `.m` scripts — no toolboxes.
- Simulink (any recent version) **only** for `microgrid_lcl.slx`.

## Quick start
```matlab
% A) Reproduce the whole study (prints a summary of every headline result):
>> run_all

% B) Simulink inverter + LCL filter (grid-current THD):
>> build_microgrid_lcl_slx   % (re)builds microgrid_lcl.slx from core blocks
>> run_thd_simulink          % runs it and reports THD + spectrum
```
The prebuilt `microgrid_lcl.slx` is included; `build_microgrid_lcl_slx.m`
regenerates it from scratch.

## Files
| File | Purpose |
|------|---------|
| `run_all.m` | Master script — runs the full study and prints the summary |
| `pv_module_current.m` | Single-diode module current I(V,G,Tc) (Eq. 1–3) |
| `pv_module_mpp.m` | Module maximum-power point by I–V scan |
| `pv_array_power.m` | 6.25 kW array MPP power (fast persistent lookup) |
| `clear_sky_day.m` | Diurnal irradiance + NOCT cell-temperature profile |
| `residential_load.m` | AC-dominated residential load profile |
| `ems_dispatch.m` | Priority-dispatch EMS to cyclic steady state (Eq. 17–19) |
| `thd_lcl_timedomain.m` | Time-domain PWM+LCL, grid-current THD by FFT (Sec. 5.4) |
| `technoeconomics.m` | LCOE, payback, NPV, degradation, tariff sweep (Sec. 5.8) |
| `riyadh_day.m` | Build one representative day from measured monthly-mean data |
| `riyadh_monthly.m` | Twelve-month Riyadh analysis; edit GHI/Tamb vectors here |
| `build_microgrid_lcl_slx.m` | Builds the Simulink inverter+LCL model (core blocks) |
| `run_thd_simulink.m` | Runs the Simulink model and computes THD |
| `microgrid_lcl.slx` | Prebuilt Simulink model |
