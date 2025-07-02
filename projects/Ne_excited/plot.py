import numpy as np
import matplotlib.pyplot as plt

dir = 'xsec'
initial = 1
final = 1
title = 'He, doublet A'
xlabel = 'Energy (eV)'
ylabel = 'Cross section (m²)'
e_unit = 1.0

# a.u. -> m² 変換係数
bohr_to_meter = 0.529177210903e-10
au_to_m2 = bohr_to_meter ** 2

# cross section 単位変換
x_unit = au_to_m2

filename = f"{dir}/xsec.doublet.A.from_initial_state_{initial}.geom1"

data = np.loadtxt(filename)

energy = data[:, 0] * e_unit
cross_section = data[:, final + 1] * x_unit

plt.figure()
plt.plot(energy, cross_section, label='He atom')
plt.xscale('log')
plt.yscale('log')
plt.title(title)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend()
plt.grid(True)
plt.savefig("xsec")

