import numpy as np, json, microgrid_core as m

# --- Riyadh (24.71N, 46.67E) monthly climatology (citable: K.A.CARE atlas / Global Solar Atlas / NASA POWER) ---
months=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ndays=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
GHI  =np.array([4.0,4.9,5.8,6.5,7.4,8.0,7.7,7.2,6.5,5.5,4.4,3.8])   # kWh/m2/day
Tamb =np.array([14.5,17.5,21.5,27.5,33.0,36.0,37.0,36.5,33.0,27.0,20.5,15.5]) # degC
daylen=np.array([10.6,11.2,12.0,12.9,13.5,13.8,13.7,13.1,12.3,11.4,10.7,10.4]) # h

# monthly residential load (AC-driven): base + cooling term
base_load=28.0; kcool=1.95
load_kwh=base_load+kcool*np.clip(Tamb-18,0,None)

eta=np.sqrt(0.95)
def month_day(GHI_t, Tam, dl, load, shape_exp=1.4):
    n=1440; t=np.arange(n)/60.0
    sr=12-dl/2; ss=12+dl/2
    g=np.zeros(n); day=(t>sr)&(t<ss)
    g[day]=np.sin(np.pi*(t[day]-sr)/(ss-sr))**shape_exp
    area=np.sum(g)/60.0                      # h (peak=1)
    Gpk=GHI_t*1000/area
    G=g*Gpk
    Tc=Tam+(45-20)/800.0*G
    Ppv=m.pv_power_series(G,Tc)              # kW
    Pl=m.residential_load(t,load)
    r,_=m.dispatch(t,Ppv,Pl,soc_min=0.40,soc_max=0.80,slew=500.0)
    return dict(pv=float(np.sum(Ppv)/60), peak=float(Ppv.max()), Gpk=float(Gpk),
                ss=r['self_suff'], sc=r['self_cons'], imp=r['Egrid_imp'], exp=r['Egrid_exp'],
                ld=float(np.sum(Pl)/60), bat=r['Ebat2load'])

rows=[]; 
for i in range(12):
    d=month_day(GHI[i],Tamb[i],daylen[i],load_kwh[i]); d['month']=months[i]; d['days']=int(ndays[i]); rows.append(d)

# annual aggregates (weight by days)
def wsum(key): return sum(r[key]*r['days'] for r in rows)
ann_pv=wsum('pv'); ann_ld=wsum('ld'); ann_imp=wsum('imp'); ann_exp=wsum('exp'); ann_bat=wsum('bat')
ann_ss=100*(1-ann_imp/ann_ld)
ann_sc=100*(ann_pv-ann_exp)/ann_pv  # PV used locally / PV gen (approx: gen-export)
print("=== RIYADH 12-MONTH ANALYSIS (optimised EMS, 6.25 kW PV, 13.5 kWh batt) ===")
print(f"{'Mon':4} {'GHI':>4} {'Gpk':>5} {'PVkWh':>6} {'peak':>5} {'load':>5} {'SS%':>5} {'SC%':>5} {'imp':>5}")
for r in rows:
    print(f"{r['month']:4} {GHI[months.index(r['month'])]:>4} {r['Gpk']:>5.0f} {r['pv']:>6.1f} {r['peak']:>5.2f} {r['ld']:>5.1f} {r['ss']:>5.1f} {r['sc']:>5.1f} {r['imp']:>5.1f}")
print(f"\nAnnual: PV={ann_pv:.0f} kWh, load={ann_ld:.0f} kWh, import={ann_imp:.0f}, export={ann_exp:.0f}")
print(f"Annual self-sufficiency={ann_ss:.1f}%  PV self-consumption={ann_sc:.1f}%")
print(f"Avg daily PV={ann_pv/365:.1f} kWh  Avg daily load={ann_ld/365:.1f} kWh")
print(f"Peak-sun-hours annual avg={np.mean(GHI):.2f} kWh/m2/day -> {np.sum(GHI*ndays):.0f} kWh/m2/yr")

# economics with SEC tariff
SARUSD=1/3.75; tariff=0.18*SARUSD
cpv=800*6.25; cbat=350*13.5; cinv=250*6.25; bos=0.18*(cpv+cbat+cinv); capex=cpv+cbat+cinv+bos
disc=0.05; life=25; crf=disc*(1+disc)**life/((1+disc)**life-1); om=0.01*capex; br=0.7*cbat
lcoe=(capex*crf+om+br/life)/(ann_pv-ann_exp)
def econ(tar):
    sav=tar*(ann_ld-ann_imp+ann_exp); pb=capex/sav
    npv=-capex-br/(1+disc)**13+sum(sav/(1+disc)**y for y in range(1,life+1))
    return sav,pb,npv
sav,pb,npv=econ(tariff)
print(f"\nEconomics (SEC 0.18 SAR=$0.048/kWh): CAPEX=${capex:.0f} LCOE=${lcoe:.3f} saving=${sav:.0f}/yr payback={pb:.1f}yr NPV=${npv:.0f}")
for lbl,tar in [('0.18 SAR',0.048),('0.30 SAR',0.080),('0.32 SAR cost-refl',0.12),('unsub 0.16',0.16)]:
    s,p,nv=econ(tar); print(f"   {lbl:20} ${tar:.3f} -> payback {p:.1f}yr NPV ${nv:.0f}")

json.dump({'rows':rows,'annual':dict(pv=ann_pv,ld=ann_ld,imp=ann_imp,exp=ann_exp,ss=ann_ss,sc=ann_sc),
           'GHI':GHI.tolist(),'Tamb':Tamb.tolist(),'months':months,'econ':dict(capex=capex,lcoe=lcoe,payback=pb,npv=npv)},
          open('results_riyadh.json','w'),indent=2)

# ---- representative annual-mean day (for Table 10 + dispatch figures) ----
def build_day(GHI_t,Tam,dl,load,shape_exp=1.4):
    n=1440; t=np.arange(n)/60.0; sr=12-dl/2; ss=12+dl/2
    g=np.zeros(n); day=(t>sr)&(t<ss); g[day]=np.sin(np.pi*(t[day]-sr)/(ss-sr))**shape_exp
    Gpk=GHI_t*1000/(np.sum(g)/60.0); G=g*Gpk; Tc=Tam+25/800.0*G
    Ppv=m.pv_power_series(G,Tc); Pl=m.residential_load(t,load)
    return t,G,Tc,Ppv,Pl
GHIm=float(np.mean(GHI)); Tm=float(np.mean(Tamb)); dlm=float(np.mean(daylen)); Lm=base_load+kcool*max(Tm-18,0)
t,G,Tc,Ppv,Pl=build_day(GHIm,Tm,dlm,Lm)
b,bs=m.dispatch(t,Ppv,Pl,soc_min=0.50,soc_max=0.80,slew=None)
o,os=m.dispatch(t,Ppv,Pl,soc_min=0.40,soc_max=0.80,slew=500.0)
print("\n=== REPRESENTATIVE (annual-mean) RIYADH DAY ===")
print(f"GHI={GHIm:.2f} kWh/m2/day Tamb={Tm:.1f}C load={Lm:.1f} kWh  dailyPV={np.sum(Ppv)/60:.1f} peak={Ppv.max():.2f} Tc_noon={Tc[720]:.0f}")
for nm,r in [('baseline',b),('optimised',o)]:
    print(f"  {nm}: SS {r['self_suff']:.1f}% SC {r['self_cons']:.1f}% PV->load {r['Epv2load']:.1f} PV->bat {r['Epv2bat']:.1f} bat->load {r['Ebat2load']:.1f} imp {r['Egrid_imp']:.1f} exp {r['Egrid_exp']:.1f} peakImp {r['Pimp_peak']:.2f}")
np.savez('riyadh_series.npz',t=t,G=G,Tc=Tc,Ppv=Ppv,Pload=Pl,
         base_Pb=bs['Pb'],base_Pgrid=bs['Pgrid'],base_SOC=bs['SOC'],
         opt_Pb=os['Pb'],opt_Pgrid=os['Pgrid'],opt_SOC=os['SOC'])
# save monthly arrays for figure
np.savez('riyadh_monthly.npz',months=months,GHI=GHI,Tamb=Tamb,
         pv=[r['pv'] for r in rows],ss=[r['ss'] for r in rows],sc=[r['sc'] for r in rows],
         imp=[r['imp'] for r in rows],ld=[r['ld'] for r in rows],peak=[r['peak'] for r in rows])
print(f"\nrep-day GHIm={GHIm:.2f} Tm={Tm:.1f} Lm={Lm:.1f}")
