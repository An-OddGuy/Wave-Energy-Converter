import numpy as np, pandas as pd
Hs, Tp, gam, N = 0.04, 0.8, 3.3, 20
g, wp = 9.81, 2*np.pi/0.8
d = pd.read_csv("results/mat_sweep.csv")
w0, w1 = d.omega_rad_s.min(), d.omega_rad_s.max()   # Fe table range
dw = (w1 - w0)/N
w = w0 + dw*(np.arange(N) + 0.5)
sig = np.where(w <= wp, 0.07, 0.09)
S = w**-5*np.exp(-1.25*(wp/w)**4)*gam**np.exp(-(w-wp)**2/(2*sig**2*wp**2))
a = np.sqrt(2*S*dw)
a *= (Hs/4)/np.sqrt(np.sum(a**2)/2)      # band-limited: total Hs = 0.04 m
F = np.interp(w, d.omega_rad_s, d.Fe_no_mat)
ph = np.random.default_rng(1).uniform(0, 2*np.pi, N)
f = lambda x: "{" + ", ".join(f"{v:.5f}" for v in x) + "}"
print(f"  parameter Integer N = {N};")
print(f"  parameter Real wi[N] = {f(w)};")
print(f"  parameter Real ai[N] = {f(a)};")
print(f"  parameter Real Fi[N] = {f(F)};")
print(f"  parameter Real ph[N] = {f(ph)};")