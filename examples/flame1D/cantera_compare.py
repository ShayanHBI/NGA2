import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

# gas = ct.Solution('cti/gri30.cti')
# gas.set_equivalence_ratio(1.0, "CH4", "O2:1.0,N2:3.76")
# gas.TP = 300, 101325

gas = ct.Solution('cti/reducedS152R621_0.cti')
gas.set_equivalence_ratio(1.0, "XC12H26:0.8649, HMN:0.1351", "O2:1.0,N2:3.76")
gas.TP = 700.0, 3.4e6

# width = 0.2
width = 3e-4
f = ct.FreeFlame(gas, width=width)
f.transport_model = 'Mix'

f.solve(refine_grid='refine')

print(f"{'Flame speed':<20} : {f.u[0]:.6f}")
print(f"{'Flame location':<20} : {f.grid[np.argmax(f.X[gas.species_index('OH'), :])]:.6f}")

label = "Unburned"
print(f"{label} {'temperature':<10} : {f.T[0]:.3f}")
for i, name in enumerate(f.gas.species_names):
    Y = f.Y[i, 0]
    if Y > 1e-6:
        print(f"{label} {name:<11} : {Y:.6f}")

label = "Burned"
print(f"{label} {'temperature':<13} : {f.T[-1]:.3f}")
for i, name in enumerate(f.gas.species_names):
    Y = f.Y[i, -1]
    if Y > 1e-6:
        print(f"{label} {name:<13} : {Y:.6f}")

# plt.figure()
# plt.plot(f.grid, f.u)
# plt.xlabel('x (m)')
# plt.ylabel('u (m/s)')
# plt.show()

# plt.figure()
# plt.plot(f.grid, f.T)
# plt.xlabel('x (m)')
# plt.ylabel('T (K)')
# plt.show()

# plt.figure()
# plt.plot(f.grid, f.Y[gas.species_index('OH'), :])
# plt.xlabel('x (m)')
# plt.ylabel('Y_OH')
# plt.show()

# Load data
data = np.loadtxt('./NGA2.dat', skiprows=1)
x       = data[:, 0] * 1000
XC12H26 = data[:, 1]
HMN     = data[:, 2]
N2      = data[:, 3]
OH      = data[:, 4]
CO      = data[:, 5]
T       = data[:, 6]
u       = data[:, 7]

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Plot

# XC12H26
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(f.grid*1000, f.Y[gas.species_index('XC12H26'), :], ls='-', lw=2, color='k')
plt.plot(x, XC12H26, ls='-', lw=2, color='b')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$x~(mm)$', fontsize=12)
plt.ylabel(r'$Y_{XC12H26}$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$CANTERA$', r'$NGA2$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./flame_profiles/XC12H26.pdf')

# HMN
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(f.grid*1000, f.Y[gas.species_index('HMN'), :], ls='-', lw=2, color='k')
plt.plot(x, HMN, ls='-', lw=2, color='b')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$x~(mm)$', fontsize=12)
plt.ylabel(r'$Y_{HMN}$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$CANTERA$', r'$NGA2$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./flame_profiles/HMN.pdf')

# N2
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(f.grid*1000, f.Y[gas.species_index('N2'), :], ls='-', lw=2, color='k')
plt.plot(x, N2, ls='-', lw=2, color='b')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$x~(mm)$', fontsize=12)
plt.ylabel(r'$Y_{N2}$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$CANTERA$', r'$NGA2$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./flame_profiles/N2.pdf')

# OH
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(f.grid*1000, f.Y[gas.species_index('OH'), :], ls='-', lw=2, color='k')
plt.plot(x, OH, ls='-', lw=2, color='b')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$x~(mm)$', fontsize=12)
plt.ylabel(r'$Y_{OH}$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$CANTERA$', r'$NGA2$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./flame_profiles/OH.pdf')

# CO
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(f.grid*1000, f.Y[gas.species_index('CO'), :], ls='-', lw=2, color='k')
plt.plot(x, CO, ls='-', lw=2, color='b')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$x~(mm)$', fontsize=12)
plt.ylabel(r'$Y_{CO}$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$CANTERA$', r'$NGA2$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./flame_profiles/CO.pdf')

# T
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(f.grid*1000, f.T, ls='-', lw=2, color='k')
plt.plot(x, T, ls='-', lw=2, color='b')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$x~(mm)$', fontsize=12)
plt.ylabel(r'$T~(K)$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$CANTERA$', r'$NGA2$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./flame_profiles/T.pdf')

# u
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(f.grid*1000, f.u, ls='-', lw=2, color='k')
plt.plot(x, u, ls='-', lw=2, color='b')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$x~(mm)$', fontsize=12)
plt.ylabel(r'$u~(m/s)$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$CANTERA$', r'$NGA2$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./flame_profiles/u.pdf')