"""
Grid-connected PV-battery microgrid: theoretical + verification model.
Core numeric engine: single-diode PV, 24h EMS dispatch, seasonal & sensitivity,
techno-economics. Produces machine-readable results consumed by figures + paper.
"""
import numpy as np
from scipy.optimize import brentq
import json

# ---------- physical constants ----------
Q = 1.602176634e-19
KB = 1.380649e-23

# ---------- module parameters (250 W, 60-cell mono) : Table 1 ----------
NS_CELL = 60
ISC_STC = 8.68        # A
VOC_STC = 37.51       # V
A_IDEAL = 1.18
RS = 0.1865           # ohm
RSH = 433.6           # ohm
KI = 0.0045           # A/degC  (+0.05%/degC)
KV = -0.117           # V/degC  (-0.31%/degC)
G_STC = 1000.0
T_STC = 25.0
NSER, NPAR = 5, 5
NMOD = NSER * NPAR    # 25 modules -> 6.25 kW array

# ---------- single-diode solver ----------
def _module_current(V, G, Tc):
    """Terminal current of one module at voltage V, irradiance G, cell temp Tc(C)."""
    if G <= 0:
        return 0.0
    T = Tc + 273.15
    Vt = NS_CELL * KB * T / Q
    Iph = (ISC_STC + KI * (Tc - T_STC)) * G / G_STC
    Voc_T = VOC_STC + KV * (Tc - T_STC)
    I0 = (ISC_STC + KI * (Tc - T_STC)) / (np.exp(Voc_T / (A_IDEAL * Vt)) - 1.0)

    def f(I):
        return Iph - I0 * (np.exp((V + I * RS) / (A_IDEAL * Vt)) - 1.0) - (V + I * RS) / RSH - I
    B = abs(Iph) + 50.0            # f is monotonic in I; this bracket always brackets the root
    try:
        return brentq(f, -B, B, xtol=1e-9, maxiter=200)
    except ValueError:
        return 0.0

def module_mpp(G, Tc, n=400):
    """Return (Vmpp, Impp, Pmpp, Voc, Isc) of one module by scanning the I-V curve."""
    if G <= 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0
    Voc_T = VOC_STC + KV * (Tc - T_STC)
    Vs = np.linspace(0.0, max(Voc_T, 1.0), n)
    Is = np.array([_module_current(v, G, Tc) for v in Vs])
    Is = np.clip(Is, 0.0, None)
    Ps = Vs * Is
    j = int(np.argmax(Ps))
    Isc = _module_current(0.0, G, Tc)
    return Vs[j], Is[j], Ps[j], Voc_T, Isc

def array_mpp_power(G, Tc):
    return NMOD * module_mpp(G, Tc)[2]

# ---------- verify module / array at STC and characteristic table ----------
def pv_characteristic_table():
    rows = []
    for G in (200, 400, 600, 800, 1000):
        for Tc in (25, 45):
            Vm, Im, Pm, Voc, Isc = module_mpp(G, Tc)
            rows.append(dict(G=G, Tc=Tc, Voc=round(Voc, 2), Isc=round(Isc, 2),
                             Vmpp=round(Vm, 2), Impp=round(Im, 2),
                             Pmod=round(Pm, 1), Parr_kW=round(NMOD * Pm / 1000.0, 2)))
    return rows

# ---------- diurnal input profiles ----------
def clear_sky_day(dt_min=1, peak_G=1000.0, sunrise=6.5, sunset=17.5, shape_exp=1.5,
                  Tamb_min=18.0, Tamb_max=30.0, NOCT=45.0):
    """Return t(h), G(W/m2), Tcell(C) for a clear-sky day.
    shape_exp>1 narrows the bell so the daily insolation matches a realistic
    ~5.4 peak-sun-hours while retaining a 1000 W/m2 solar-noon peak."""
    n = int(24 * 60 / dt_min)
    t = np.arange(n) * dt_min / 60.0
    G = np.zeros(n)
    day = (t > sunrise) & (t < sunset)
    G[day] = peak_G * np.sin(np.pi * (t[day] - sunrise) / (sunset - sunrise)) ** shape_exp
    G = np.clip(G, 0.0, None)
    # ambient: diurnal sinusoid peaking ~15:00
    Tamb = Tamb_min + (Tamb_max - Tamb_min) * np.clip(
        np.sin(np.pi * (t - (sunrise - 1.0)) / ((sunset + 2.0) - (sunrise - 1.0))), 0, None)
    Tcell = Tamb + (NOCT - 20.0) / 800.0 * G
    return t, G, Tcell

def residential_load(t, total_kwh=54.4, base=1.5, evening_peak=3.6):
    """AC-dominated residential profile (kW), scaled to a daily energy."""
    # shape: base + morning bump + strong afternoon/evening cooling extending past sunset
    shape = base * np.ones_like(t)
    shape += 0.9 * np.exp(-((t - 7.5) ** 2) / (2 * 1.1 ** 2))          # morning
    shape += 2.1 * np.exp(-((t - 15.5) ** 2) / (2 * 2.2 ** 2))         # afternoon AC
    shape += 1.9 * np.exp(-((t - 20.5) ** 2) / (2 * 1.8 ** 2))         # evening AC
    shape = np.clip(shape, base * 0.8, None)
    dt = (t[1] - t[0])
    e = np.sum(shape) * dt
    shape *= total_kwh / e
    return shape

# ---------- PV generation over a day (ideal-MPP assumption) ----------
def pv_power_series(G, Tcell):
    return np.array([array_mpp_power(g, tc) for g, tc in zip(G, Tcell)]) / 1000.0  # kW

# ---------- EMS dispatch to cyclic steady state ----------
def dispatch(t, Ppv, Pload, Ebat_nom=13.5, soc_min=0.40, soc_max=0.80,
             eta_c=0.9747, eta_d=0.9747, Pmax=6.75, slew=None, soc0=None,
             max_days=60, tol=1e-4):
    """Priority EMS dispatch; iterate day until end SOC == start SOC (cyclic steady state)."""
    dt = t[1] - t[0]                      # hours
    n = len(t)
    if soc0 is None:
        soc0 = 0.5 * (soc_min + soc_max)
    soc_start = soc0
    for _ in range(max_days):
        soc = soc_start
        SOC = np.zeros(n); Pb = np.zeros(n); Pgrid = np.zeros(n)
        pv2load = np.zeros(n); pv2bat = np.zeros(n); bat2load = np.zeros(n)
        grid2load = np.zeros(n); pv2grid = np.zeros(n)
        prev_pb = 0.0
        for i in range(n):
            SOC[i] = soc
            net = Ppv[i] - Pload[i]
            headroom = max(0.0, (soc_max - soc)) * Ebat_nom / (eta_c * dt)   # kW absorbable
            avail = max(0.0, (soc - soc_min)) * Ebat_nom * eta_d / dt        # kW deliverable
            if net >= 0:
                pch = min(net, Pmax, headroom)
                if slew is not None:  # W/s -> kW/step
                    pch = min(pch, prev_pb + slew * 3600.0 * dt / 1000.0) if prev_pb >= 0 else pch
                pdis = 0.0
                pv2load[i] = Pload[i]
                pv2bat[i] = pch
                pv2grid[i] = net - pch          # export
                Pgrid[i] = -(net - pch)
            else:
                pdis = min(-net, Pmax, avail)
                if slew is not None:
                    step = slew * 3600.0 * dt / 1000.0
                    pdis = min(pdis, abs(prev_pb) + step) if prev_pb <= 0 else pdis
                pch = 0.0
                pv2load[i] = Ppv[i]
                bat2load[i] = pdis
                grid2load[i] = (-net) - pdis    # import
                Pgrid[i] = (-net) - pdis
            Pb[i] = pch - pdis                  # +charge / -discharge
            prev_pb = Pb[i]
            soc += (eta_c * pch - pdis / eta_d) / Ebat_nom * dt
            soc = min(max(soc, soc_min), soc_max)
        if abs(soc - soc_start) < tol:
            break
        soc_start = soc                          # converge to periodic SOC
    E = lambda x: float(np.sum(x) * dt)
    res = dict(
        Epv=E(Ppv), Eload=E(Pload),
        Epv2load=E(pv2load), Epv2bat=E(pv2bat), Ebat2load=E(bat2load),
        Egrid_imp=E(grid2load), Egrid_exp=E(pv2grid),
        Ppv_peak=float(np.max(Ppv)), Pimp_peak=float(np.max(np.clip(Pgrid, 0, None))),
        soc_min_reached=float(np.min(SOC)), soc_max_reached=float(np.max(SOC)),
    )
    res['self_suff'] = 100.0 * (1 - res['Egrid_imp'] / res['Eload'])
    res['self_cons'] = 100.0 * (res['Epv2load'] + res['Epv2bat']) / res['Epv']
    res['load_share_pv'] = 100.0 * res['Epv2load'] / res['Eload']
    res['load_share_bat'] = 100.0 * res['Ebat2load'] / res['Eload']
    res['load_share_grid'] = 100.0 * res['Egrid_imp'] / res['Eload']
    series = dict(t=t, Ppv=Ppv, Pload=Pload, Pb=Pb, Pgrid=Pgrid, SOC=SOC)
    return res, series

# =====================================================================
if __name__ == "__main__":
    out = {}

    # ---- 1. STC verification ----
    Vm, Im, Pm, Voc, Isc = module_mpp(G_STC, T_STC)
    out['stc'] = dict(Vmpp=round(Vm, 2), Impp=round(Im, 2), Pmod=round(Pm, 2),
                      Parr_kW=round(NMOD * Pm / 1000.0, 3),
                      nameplate_err_pct=round(abs(NMOD * Pm - 6250) / 6250 * 100, 3))
    # 30 C reference peak
    out['p_30C_kW'] = round(array_mpp_power(1000.0, 30.0) / 1000.0, 3)
    out['p_45C_kW'] = round(array_mpp_power(1000.0, 45.0) / 1000.0, 3)
    out['derate_25_45_pct'] = round((array_mpp_power(1000, 25) - array_mpp_power(1000, 45))
                                     / array_mpp_power(1000, 25) * 100, 2)
    out['pv_table'] = pv_characteristic_table()

    # ---- 2. reference-day dispatch (baseline & optimised) ----
    t, G, Tc = clear_sky_day()
    Ppv = pv_power_series(G, Tc)
    Pload = residential_load(t, total_kwh=54.4)
    out['Epv_day_kWh'] = round(float(np.sum(Ppv) * (t[1]-t[0])), 2)
    out['Ppv_peak_kW'] = round(float(np.max(Ppv)), 2)
    out['Tcell_noon'] = round(float(Tc[720]), 1)

    base_res, base_series = dispatch(t, Ppv, Pload, soc_min=0.50, soc_max=0.80, slew=None,
                                     eta_c=0.9747, eta_d=0.9747)
    opt_res, opt_series = dispatch(t, Ppv, Pload, soc_min=0.40, soc_max=0.80, slew=500.0,
                                   eta_c=0.9747, eta_d=0.9747)
    out['baseline'] = {k: round(v, 3) for k, v in base_res.items()}
    out['optimised'] = {k: round(v, 3) for k, v in opt_res.items()}

    # save reference series for figures
    np.savez('ref_series.npz', t=t, G=G, Tc=Tc, Ppv=Ppv, Pload=Pload,
             base_Pb=base_series['Pb'], base_Pgrid=base_series['Pgrid'], base_SOC=base_series['SOC'],
             opt_Pb=opt_series['Pb'], opt_Pgrid=opt_series['Pgrid'], opt_SOC=opt_series['SOC'])

    print(json.dumps(out, indent=2))
    with open('results_core.json', 'w') as f:
        json.dump(out, f, indent=2)
