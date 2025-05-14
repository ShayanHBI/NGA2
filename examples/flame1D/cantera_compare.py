import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

# gas = ct.Solution('cti/gri30.cti')
# gas.set_equivalence_ratio(1.0, "CH4", "O2:1.0,N2:3.76")
# gas.TP = 300, 101325

gas = ct.Solution('cti/reducedS152R621_0.cti')
gas.set_equivalence_ratio(1.0, "XC12H26:0.8649, HMN:0.1351", "O2:1.0,N2:3.76")
gas.TP = 700.0, 3.4e6
width = 3e-4

# width = 0.2
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

plt.figure()
plt.plot(f.grid, f.u)
plt.xlabel('x (m)')
plt.ylabel('u (m/s)')
plt.show()

plt.figure()
plt.plot(f.grid, f.T)
plt.xlabel('x (m)')
plt.ylabel('T (K)')
plt.show()