import cantera as ct
import numpy as np


import matplotlib.pyplot as plt
import plot_params
from plot_params import lw, ls, ms, mk, c
from matplotlib.lines import Line2D



file = "monitor/simulation"
with open(file, 'r') as f:
    lines = f.readlines()

names = lines[0].split()
lines.pop(0)
lines.pop(0)

data = []
for line in lines:
    split_line = line.split()
    data.append([float(val) for val in split_line])
data = np.array(data)

time_array = data[:,names.index('Time')]
T_array = data[:,names.index('Temperature')]
Y_O2_array = data[:,names.index('Y_O2')]
Y_N2_array = data[:,names.index('Y_N2')]
Y_HMN_array = data[:,names.index('Y_HMN')]

gas = ct.Solution('cti/reducedS152R621_0.cti')

gas.set_equivalence_ratio(1.0, "XC12H26:0.8649, HMN:0.1351", "O2:1.0,N2:3.76")
gas.TP = 750.0, 3.4e6

print()
print("Initial {:10}    :    {:.7e}".format("enthalpy", gas.enthalpy_mass))
for i,spec in enumerate(gas.species_names):
    if gas.Y[i] > 1.0e-10: 
        print("Initial {:10}    :    {:.7e}".format(spec, gas.Y[i]))
print()


r = ct.IdealGasConstPressureReactor(gas)


sim = ct.ReactorNet([r])

states = ct.SolutionArray(gas, extra='t')
states.append(r.thermo.state, t=0.0)

time = 0.0
dt = 1.0e-6

for i in range(1,10000):
    sim.advance(time + dt)
    time += dt

    states.append(r.thermo.state, t=time)

fig, ax = plt.subplots(1, 2, figsize=(10, 5))

ax[0].plot(
    states.t * 1000,
    states.T,
    lw=lw,
    label="Cantera",
    color='blue',
)

ax[0].plot(
    time_array * 1000,
    T_array,
    lw=lw,
    ls="--",
    label="NGA2",
    color='red'
)

ax[0].set(
    xlabel="Time [ms]",
    ylabel="Temperature [K]",
    xlim=[0, 10],
    ylim=[750, 2600]
)

ax[0].grid()

ax[0].legend(frameon=False, loc="upper left")

ax[1].plot(
    states.t * 1000,
    states.Y[:, gas.species_index('O2')],
    lw=lw,
    label="Cantera",
    color='blue',
)

ax[1].plot(
    time_array * 1000,
    Y_O2_array,
    lw=lw,
    ls="--",
    label="NGA2",
    color='red'
)

ax[1].set(
    xlabel="Time [ms]",
    ylabel="O2 mass fraction [1]",
    xlim=[0, 10]
)

ax[1].grid()

# ax[1, 0].plot(
#     states.t * 1000,
#     states.Y[:, gas.species_index('N2')],
#     lw=lw,
#     label="Cantera",
#     color='blue',
# )

# ax[1, 0].plot(
#     time_array * 1000,
#     Y_N2_array,
#     lw=lw,
#     ls="--",
#     label="NGA2",
#     color='red'
# )

# ax[1, 0].set(
#     xlabel="Time [ms]",
#     ylabel="N2 mass fraction [1]",
#     xlim=[0, 10]
# )

# ax[1, 0].ticklabel_format(axis='y', style='sci', scilimits=(0,0))

# ax[1, 1].plot(
#     states.t * 1000,
#     states.Y[:, gas.species_index('HMN')],
#     lw=lw,
#     label="Cantera",
#     color='blue',
# )

# ax[1, 1].plot(
#     time_array * 1000,
#     Y_HMN_array,
#     lw=lw,
#     ls="--",
#     label="NGA2",
#     color='red'
# )

# ax[1, 1].set(
#     xlabel="Time [ms]",
#     ylabel="HMN mass fraction [1]",
#     xlim=[0, 10]
# )

# ax[1, 1].ticklabel_format(axis='y', style='sci', scilimits=(0,0))

# ax[1, 1].grid()

plt.tight_layout()

plt.show()
