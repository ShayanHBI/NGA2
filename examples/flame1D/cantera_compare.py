import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

gas = ct.Solution('cti/gri30.cti')
gas.set_equivalence_ratio(1.0, "CH4", "O2:1.0,N2:3.76")
gas.TP = 300, 101325

width = 0.2
f = ct.FreeFlame(gas, width=width)
f.transport_model = 'Mix'

f.solve(refine_grid='refine')

print('Flame speed = {0:7f} m/s'.format(f.u[0]))
print('Flame location = {0:7f} m/s'.format(f.grid[np.argmax(f.X[gas.species_index('OH'), :])]))

label = "Unburned"
for i, name in enumerate(f.gas.species_names):
    Y = f.Y[i, 0]
    if Y > 1e-6:
        print(f"{label} {name:<10} : {Y:.6f}")

label = "Burned"
for i, name in enumerate(f.gas.species_names):
    Y = f.Y[i, -1]
    if Y > 1e-6:
        print(f"{label} {name:<10} : {Y:.6f}")

plt.figure()
plt.plot(f.grid, f.u)
plt.xlabel('x (m)')
plt.ylabel('u (m/s)')
plt.show()