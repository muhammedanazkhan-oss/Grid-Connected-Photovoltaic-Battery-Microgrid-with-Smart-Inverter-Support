import json
d=json.load(open('results_ext.json'))
capex=d['econ']['capex_USD']; a=d['annual']
al, ai, ae = a['load_kWh'], a['grid_imp_kWh'], a['grid_exp_kWh']
disc=0.05; life=25; batt_replace=0.7*d['econ']['capex_breakdown']['batt']
rows=[]
for name,tar in [('Current subsidised (0.18 SAR)',0.048),('Partial reform (0.30 SAR)',0.080),
                 ('Cost-reflective (~0.45 SAR)',0.120),('Unsubsidised (~0.60 SAR)',0.160)]:
    saving=tar*(al-ai+ae)
    payback=capex/saving
    npv=-capex+sum(saving/(1+disc)**y for y in range(1,life+1))-batt_replace/(1+disc)**13
    rows.append(dict(scenario=name, tariff_USD=tar, saving_USD=round(saving),
                     payback_yr=round(payback,1), NPV_USD=round(npv)))
d['econ']['tariff_sensitivity']=rows
json.dump(d,open('results_ext.json','w'),indent=2)
for r in rows: print(f"{r['scenario']:32s} tariff ${r['tariff_USD']:.3f}  saving ${r['saving_USD']:4d}/yr  payback {r['payback_yr']:5.1f} yr  NPV ${r['NPV_USD']}")
