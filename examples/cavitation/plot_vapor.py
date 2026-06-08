import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('monitor/simulation', skiprows=3)
time      = data[:,1]
VFint     = data[:,17]
vap_mass  = data[:,20]
Yv_int_r  = data[:,21]

fig, axes = plt.subplots(3, 1, figsize=(7,10), sharex=True)

axes[0].plot(time, vap_mass, '-o', ms=2)
axes[0].set_ylabel('Vapor mass [kg]')
axes[0].set_title('Vapor mass vs time')
axes[0].grid(True)

axes[1].plot(time, VFint, '-o', ms=2, color='C1')
axes[1].set_ylabel('VFint [m$^3$]')
axes[1].set_title('VFint (liquid volume integral) vs time')
axes[1].grid(True)

axes[2].plot(time, Yv_int_r, '-o', ms=2, color='C2')
axes[2].set_ylabel(r'$\int Y_v\,dV$ ($r\leq0.01$) [m$^3$]')
axes[2].set_xlabel('Time [s]')
axes[2].set_title('Yv integral within r=0.01 vs time')
axes[2].grid(True)

fig.tight_layout()
fig.savefig('vapor_mass_VFint.pdf')
fig.savefig('vapor_mass_VFint.png', dpi=150)
print('wrote vapor_mass_VFint.pdf and .png')
