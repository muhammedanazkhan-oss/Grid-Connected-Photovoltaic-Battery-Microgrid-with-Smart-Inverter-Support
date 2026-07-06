"""Seasonal/robustness, sensitivity, reproducibility, and techno-economic analyses."""
import numpy as np, json
import microgrid_core as m

dt_h = 1/60.0

def run_day(peak_G, sr, ss, sexp, Tmin, Tmax, load_kwh, soc_min=0.40, soc_max=0.80, slew=500.0):
    t,G,Tc = m.clear_sky_day(peak_G=peak_G, sunrise=sr, sunset=ss, shape_exp=sexp,
                             Tamb_min=Tmin, Tamb_max=Tmax)
    Ppv = m.pv_power_series(G,Tc)
    Pload = m.residential_load(t, total_kwh=load_kwh)
    res,ser = m.dispatch(t,Ppv,Pload,soc_min=soc_min,soc_max=soc_max,slew=slew)
    res['Epv_day']=float(np.sum(Ppv)*dt_h); res['Ppv_peak']=float(Ppv.max())
    res['Tc_noon']=float(Tc[720])
    return res, ser, (t,G,Tc,Ppv,Pload)

OUT={}

# ---------- 1. Seasonal + robustness ----------
scen={}
# representative (reference)
r,_,_ = run_day(1000, 6.5,17.5,1.5, 18,30, 54.4); scen['representative']=r
# summer: longer day, hotter, heavy AC
r,_,_ = run_day(1000, 5.5,18.5,1.35, 26,44, 62.0); scen['summer']=r
# winter: weaker/shorter sun, cooler (better cell temp), lighter load
r,_,_ = run_day(780, 7.0,17.0,1.55, 8,20, 40.0); scen['winter']=r
# cloudy / prolonged low irradiance robustness (R2.8)
t,G,Tc = m.clear_sky_day(peak_G=1000, sunrise=6.0, sunset=18.0, shape_exp=1.3, Tamb_min=16,Tamb_max=26)
rng=np.random.default_rng(7)
cloud=np.clip(0.35+0.5*rng.random(len(t)),0,1)
# smooth the cloud factor
k=np.ones(25)/25; cloud=np.convolve(cloud,k,mode='same')
Gc=G*cloud
Ppv=m.pv_power_series(Gc,Tc); Pload=m.residential_load(t,total_kwh=54.4)
rc,serc=m.dispatch(t,Ppv,Pload,soc_min=0.40,soc_max=0.80,slew=500.0)
rc['Epv_day']=float(np.sum(Ppv)*dt_h); rc['Ppv_peak']=float(Ppv.max()); rc['Tc_noon']=float(Tc[720])
scen['cloudy']=rc
np.savez('cloudy_series.npz', t=t, G=G, Gc=Gc, Ppv=Ppv, Pload=Pload,
         Pb=serc['Pb'], Pgrid=serc['Pgrid'], SOC=serc['SOC'])
OUT['seasonal']={k:{kk:round(vv,2) for kk,vv in v.items()} for k,v in scen.items()}

# ---------- 2. Sensitivity: battery size & load ----------
t,G,Tc=m.clear_sky_day(); Ppv=m.pv_power_series(G,Tc)
batt_sweep=[]
for E in [6.75,10.0,13.5,20.0,27.0]:
    Pload=m.residential_load(t,54.4)
    r,_=m.dispatch(t,Ppv,Pload,Ebat_nom=E,soc_min=0.40,soc_max=0.80,slew=500.0)
    batt_sweep.append(dict(Ebat=E, self_suff=round(r['self_suff'],1),
                           self_cons=round(r['self_cons'],1),
                           bat2load=round(r['Ebat2load'],2),
                           grid_imp=round(r['Egrid_imp'],2)))
OUT['batt_sweep']=batt_sweep

load_sweep=[]
for L in [30.0,45.0,54.4,60.0]:
    Pload=m.residential_load(t,L)
    rb,_=m.dispatch(t,Ppv,Pload,Ebat_nom=13.5,soc_min=0.50,soc_max=0.80,slew=None)
    ro,_=m.dispatch(t,Ppv,Pload,Ebat_nom=13.5,soc_min=0.40,soc_max=0.80,slew=500.0)
    load_sweep.append(dict(load=L, ss_base=round(rb['self_suff'],1), ss_opt=round(ro['self_suff'],1),
                           sc_base=round(rb['self_cons'],1), sc_opt=round(ro['self_cons'],1)))
OUT['load_sweep']=load_sweep

# ---------- 3. Reproducibility: Python vs reference (Simulink) KPIs ----------
Pload=m.residential_load(t,54.4)
r_opt,ser_opt=m.dispatch(t,Ppv,Pload,soc_min=0.40,soc_max=0.80,slew=500.0)
# power-balance residual (energy conservation)  P_pv + P_grid + P_bat_dis = P_load + P_bat_ch
Pb=ser_opt['Pb']; Pgrid=ser_opt['Pgrid']
resid = Ppv + Pgrid + np.clip(-Pb,0,None) - Pload - np.clip(Pb,0,None)
OUT['balance_resid_max_pct']=round(float(np.max(np.abs(resid))/ (Pload.max()) *100),3)
OUT['balance_resid_rmse_W']=round(float(np.sqrt(np.mean(resid**2))*1000),3)

# Python vs reference Simulink values reported in the original study
kpis=[  # name, python, reference(Simulink)
 ("Module MPP (W)",       m.module_mpp(1000,25)[2],       250.0),
 ("Array peak @30C (kW)", m.array_mpp_power(1000,30)/1000, 6.13),
 ("Array @45C (kW)",      m.array_mpp_power(1000,45)/1000, 5.77),
 ("Daily PV (kWh)",       float(np.sum(Ppv)*dt_h),         34.0),
 ("Self-sufficiency (%)", r_opt['self_suff'],              44.6),
 ("PV self-consumption (%)", r_opt['self_cons'],           72.1),
 ("Grid import (kWh)",    r_opt['Egrid_imp'],              30.2),
 ("Battery->load (kWh)",  r_opt['Ebat2load'],              5.3),
]
errs=[abs(p-ref)/ref*100 for _,p,ref in kpis]
OUT['repro_table']=[dict(kpi=n, python=round(p,2), ref=ref, err_pct=round(e,2))
                    for (n,p,ref),e in zip(kpis,errs)]
OUT['repro_MAPE_pct']=round(float(np.mean(errs)),2)
OUT['repro_RMSE_pct']=round(float(np.sqrt(np.mean(np.square(errs)))),2)

# ---------- 4. Annual / cumulative (seasonal-weighted) ----------
days={'summer':122,'winter':121,'representative':122}
annual_pv=sum(scen[s]['Epv_day']*days[s] for s in days)
annual_load=sum(m.residential_load(t,{'summer':62,'winter':40,'representative':54.4}[s]).sum()*dt_h*days[s] for s in days)
annual_imp=sum(scen[s]['Egrid_imp']*days[s] for s in days)
annual_exp=sum(scen[s]['Egrid_exp']*days[s] for s in days)
OUT['annual']=dict(pv_kWh=round(annual_pv), load_kWh=round(annual_load),
                   grid_imp_kWh=round(annual_imp), grid_exp_kWh=round(annual_exp),
                   self_suff=round(100*(1-annual_imp/annual_load),1))

# ---------- 5. Techno-economics ----------
SARUSD=1/3.75
# CAPEX (2024 residential, USD)
pv_kw=6.25; batt_kwh=13.5
cost_pv=800*pv_kw; cost_batt=350*batt_kwh; cost_inv=250*6.25; bos=0.18*(cost_pv+cost_batt+cost_inv)
capex=cost_pv+cost_batt+cost_inv+bos
# tariff (SEC residential): 0.18 SAR/kWh -> USD
tariff=0.18*SARUSD          # 0.048 USD/kWh
export_rate=0.18*SARUSD     # net-billing at consumption tariff (assumption)
# annual flows from annual weighting
avoided_import=annual_load - annual_imp    # kWh served locally that would otherwise be imported
bill_no_pv=annual_load*tariff
bill_with=annual_imp*tariff - annual_exp*export_rate
annual_saving=bill_no_pv-bill_with
# LCOE of PV+battery generation
disc=0.05; life=25
crf=disc*(1+disc)**life/((1+disc)**life-1)
om=0.01*capex
batt_replace=cost_batt*0.7   # one replacement (real cost decline) ~yr13
batt_replace_annualized=batt_replace*disc/((1+disc)**13-1)*((1+disc)**13)/ (1+disc)**13  # approx
gen_local=annual_pv - annual_exp  # PV energy used locally (to load+batt)
lcoe=(capex*crf+om+batt_replace/ life)/gen_local
# simple payback & NPV
simple_payback=capex/annual_saving
npv=-capex+sum(annual_saving/(1+disc)**y for y in range(1,life+1)) - batt_replace/(1+disc)**13
# degradation
daily_dis=scen['representative']['Ebat2load']
annual_throughput=daily_dis*365
efc_per_yr=annual_throughput/batt_kwh
cycle_life_efc=4000  # LiFePO4 shallow-cycle to 80%
cycle_years=cycle_life_efc/efc_per_yr
calendar_years=13
batt_life=min(cycle_years,calendar_years)
OUT['econ']=dict(capex_USD=round(capex), tariff_USD_kWh=round(tariff,3),
   annual_saving_USD=round(annual_saving), simple_payback_yr=round(simple_payback,1),
   LCOE_USD_kWh=round(lcoe,3), NPV_USD=round(npv), crf=round(crf,4),
   efc_per_yr=round(efc_per_yr,1), cycle_life_yr=round(cycle_years,1),
   batt_life_yr=round(batt_life,1),
   capex_breakdown=dict(pv=round(cost_pv),batt=round(cost_batt),inv=round(cost_inv),bos=round(bos)))

print(json.dumps(OUT,indent=2))
json.dump(OUT,open('results_ext.json','w'),indent=2)
