import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

# Load data
data = np.loadtxt('./monitor/ignition', skiprows=1)
t_nga       = data[:, 0]
T_nga       = data[:, 1]
XC12H26_nga = data[:, 2]
HMN_nga     = data[:, 3]
N2_nga      = data[:, 4]
OH_nga      = data[:, 5]
CO_nga      = data[:, 6]
P_nga       = data[:, 7]

gas = ct.Solution('cti/reducedS152R621_0.cti')
gas.set_equivalence_ratio(1.0, "XC12H26:0.8649, HMN:0.1351", "O2:1.0,N2:3.76")
gas.TP = 700.0, 3.4e6

label = "Initial"
print(f"{label} {'temperature':<10} : {gas.T:.3f}")
for i, name in enumerate(gas.species_names):
    Y = gas.Y[i]
    if Y > 0:
        print(f"{label} {name:<11} : {Y:.6f}")

r = ct.IdealGasReactor(gas) # Constant volume reactor
sim = ct.ReactorNet([r])

nt = int(7e5)
dt = 1e-8
t_can       = np.zeros(nt)
T_can       = np.zeros(nt)
XC12H26_can = np.zeros(nt)
HMN_can     = np.zeros(nt)
N2_can      = np.zeros(nt)
OH_can      = np.zeros(nt)
CO_can      = np.zeros(nt)
P_can       = np.zeros(nt)
time = 0
for n in range(nt):
    time += dt
    sim.advance(time)
    t_can[n] = time
    T_can[n] = r.T
    XC12H26_can[n] = r.thermo.Y[gas.species_index('XC12H26')]
    HMN_can[n]     = r.thermo.Y[gas.species_index('HMN')]
    N2_can [n]     = r.thermo.Y[gas.species_index('N2')]
    OH_can [n]     = r.thermo.Y[gas.species_index('OH')]
    CO_can [n]     = r.thermo.Y[gas.species_index('CO')]
    P_can  [n]     = r.thermo.P

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Temperature
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(t_can, T_can, ls='-' , lw=2, color='k')
plt.plot(t_nga, T_nga, ls='--', lw=2, color='b')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$t~(s)$', fontsize=12)
plt.ylabel(r'$T~(K)$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$CANTERA$', r'$NGA2$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./ignition_profiles/T.pdf')

# XC12H26
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(t_can, XC12H26_can, ls='-' , lw=2, color='k')
plt.plot(t_nga, XC12H26_nga, ls='--', lw=2, color='b')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$t~(s)$', fontsize=12)
plt.ylabel(r'$Y_{XC12H26}$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$CANTERA$', r'$NGA2$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./ignition_profiles/XC12H26.pdf')

# HMN
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(t_can, HMN_can, ls='-' , lw=2, color='k')
plt.plot(t_nga, HMN_nga, ls='--', lw=2, color='b')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$t~(s)$', fontsize=12)
plt.ylabel(r'$Y_{HMN}$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$CANTERA$', r'$NGA2$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./ignition_profiles/HMN.pdf')

# N2
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(t_can, N2_can, ls='-' , lw=2, color='k')
plt.plot(t_nga, N2_nga, ls='--', lw=2, color='b')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$t~(s)$', fontsize=12)
plt.ylabel(r'$Y_{N2}$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$CANTERA$', r'$NGA2$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./ignition_profiles/N2.pdf')

# OH
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(t_can, OH_can, ls='-' , lw=2, color='k')
plt.plot(t_nga, OH_nga, ls='--', lw=2, color='b')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$t~(s)$', fontsize=12)
plt.ylabel(r'$Y_{OH}$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$CANTERA$', r'$NGA2$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./ignition_profiles/OH.pdf')

# CO
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(t_can, CO_can, ls='-' , lw=2, color='k')
plt.plot(t_nga, CO_nga, ls='--', lw=2, color='b')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$t~(s)$', fontsize=12)
plt.ylabel(r'$Y_{CO}$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$CANTERA$', r'$NGA2$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./ignition_profiles/CO.pdf')

# P
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plt.plot(t_can, P_can, ls='-' , lw=2, color='k')
plt.plot(t_nga, P_nga, ls='--', lw=2, color='b')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$t~(s)$', fontsize=12)
plt.ylabel(r'$p~(Pa)$', fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid(which='major', axis='both', color='k', linestyle='--', linewidth=0.6, alpha=0.25)
legends = ax.legend([r'$CANTERA$', r'$NGA2$'], frameon=False, loc='lower right', fontsize=10)
ax.add_artist(legends)
plt.tight_layout()
plt.savefig('./ignition_profiles/P.pdf')