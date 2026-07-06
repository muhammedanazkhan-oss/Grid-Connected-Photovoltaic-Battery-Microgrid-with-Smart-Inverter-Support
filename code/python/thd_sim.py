"""Time-domain single-phase full-bridge unipolar-SPWM inverter + LCL filter.
Computes grid-current THD via FFT and compares LCL vs L-only. Addresses R1.9/R2.3/R3.4."""
import numpy as np
import json

# ---- grid / inverter ----
Vg_rms = 230.0
Vg_pk  = Vg_rms*np.sqrt(2)
f0     = 60.0
w0     = 2*np.pi*f0
Vdc    = 400.0
fsw    = 10_000.0
Pn     = 6250.0
In_rms = Pn/Vg_rms          # 27.2 A
In_pk  = In_rms*np.sqrt(2)

# ---- LCL (Table 6) ----
L1 = 0.87e-3
L2 = 0.43e-3
Cf = 15.67e-6
Rd = 1.43
Lg = 30e-6                   # small stiff-grid inductance
L2t = L2 + Lg

fs   = 2_000_000.0          # 2 MHz integration
dt   = 1/fs
Tend = 6/f0                 # 6 fundamental cycles
N    = int(Tend/dt)
tt   = np.arange(N)*dt

# fundamental inverter voltage phasor to push rated current at unity PF at PCC:
# Vinv1 = Vgrid + jw(L1+L2t) I   ->  amplitude/phase of modulating ref
Ltot = L1 + L2t
Vinv1 = Vg_pk + 1j*w0*Ltot*In_pk           # peak phasor (grid as reference 0 deg)
ma    = np.abs(Vinv1)/Vdc
phi   = np.angle(Vinv1)
carrier = np.abs(2*(tt*fsw - np.floor(tt*fsw+0.5)))     # 0..1 triangle
mod = ma*np.sin(w0*tt + phi)

def run_lcl(vinv):
    iL1=np.zeros(N); iL2=np.zeros(N); vcf=np.zeros(N)
    vg = Vg_pk*np.sin(w0*tt)
    for n in range(N-1):
        inode = iL1[n]-iL2[n]
        vnode = vcf[n] + Rd*inode
        diL1 = (vinv[n]-vnode)/L1
        diL2 = (vnode-vg[n])/L2t
        dvcf = inode/Cf
        iL1[n+1]=iL1[n]+dt*diL1
        iL2[n+1]=iL2[n]+dt*diL2
        vcf[n+1]=vcf[n]+dt*dvcf
    return iL2, vg

def run_lonly(vinv):
    iL=np.zeros(N); vg=Vg_pk*np.sin(w0*tt)
    for n in range(N-1):
        iL[n+1]=iL[n]+dt*(vinv[n]-vg[n])/Ltot
    return iL, vg

# unipolar SPWM: two legs, references +mod and -mod
def unipolar(mod):
    va = np.where(mod >  carrier*1.0, Vdc/2, -Vdc/2)
    vb = np.where(-mod > carrier*1.0, Vdc/2, -Vdc/2)
    return va - vb    # bridge output -Vdc..+Vdc

vinv = unipolar(mod)

iL2, vg = run_lcl(vinv)
iL_l, _ = run_lonly(vinv)

# analyse last 3 cycles (integer cycles for clean FFT)
def thd(sig):
    ncyc=3
    Ns=int(ncyc/f0/dt)
    s=sig[-Ns:]
    S=np.fft.rfft(s*np.hanning(len(s)))
    freqs=np.fft.rfftfreq(len(s),dt)
    # fundamental bin near 60 Hz
    k1=np.argmin(np.abs(freqs-f0))
    # sum window around each bin to capture leakage
    def amp(k):
        return np.sqrt(np.sum(np.abs(S[max(k-2,0):k+3])**2))
    A1=amp(k1)
    # total distortion = everything except fundamental cluster and DC
    mask=np.ones(len(S),bool); mask[:3]=False
    for kk in range(k1-2,k1+3):
        if 0<=kk<len(mask): mask[kk]=False
    Ah=np.sqrt(np.sum(np.abs(S[mask])**2))
    I1_rms=A1/len(s)*2/np.sqrt(2)*np.sqrt(8/3)  # hanning correction approx
    return 100*Ah/A1, I1_rms

thd_lcl, i1 = thd(iL2)
thd_l, i1l  = thd(iL_l)

# fundamental grid current amplitude (robust): fit 60 Hz sine on last 3 cycles
Ns=int(3/f0/dt); s=iL2[-Ns:]; ts=tt[-Ns:]
A=np.c_[np.sin(w0*ts),np.cos(w0*ts),np.ones_like(ts)]
c,_,_,_=np.linalg.lstsq(A,s,rcond=None)
I1_pk=np.hypot(c[0],c[1]); I1_rms=I1_pk/np.sqrt(2)

# also L-only fundamental
sl=iL_l[-Ns:]; cl,_,_,_=np.linalg.lstsq(A,sl,rcond=None)
I1l_rms=np.hypot(cl[0],cl[1])/np.sqrt(2)

out=dict(ma=round(ma,3), phi_deg=round(np.degrees(phi),2),
         I1_rms_LCL=round(I1_rms,2), I1_rms_L=round(I1l_rms,2),
         THD_LCL_pct=round(thd_lcl,2), THD_L_pct=round(thd_l,1),
         fres_Hz=round(1/(2*np.pi)*np.sqrt((L1+L2t)/(L1*L2t*Cf)),0),
         IEEE519_limit_pct=5.0, pass_IEEE519=bool(thd_lcl<5.0))
print(json.dumps(out,indent=2))
json.dump(out,open('results_thd.json','w'),indent=2)
# save waveforms for figure (decimate)
d=20
np.savez('thd_wave.npz', t=tt[::d], vinv=vinv[::d], iL2=iL2[::d], iL_l=iL_l[::d],
         vg=vg[::d], fs=fs, f0=f0)
