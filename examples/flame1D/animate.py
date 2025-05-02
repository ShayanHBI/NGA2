import cantera as ct
import numpy as np
from pathlib import Path

import os
import matplotlib.pyplot as plt
import plot_params
from plot_params import lw, ls, ms, mk, c
from matplotlib.lines import Line2D

import pandas as pd


gas = ct.Solution('cti/YAO_reduced.cti')

gas.set_equivalence_ratio(1.0, "NXC12H26", "O2:1.0,N2:3.76")
gas.TP = 300.0, 3.4e6
width = 3e-4  # m
loglevel = 1  # amount of diagnostic output (0 to 8)

restore = True

# Set up flame object
f = ct.FreeFlame(gas, width=width)
f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
f.transport_model = 'multicomponent'

output = "./laminar_flame_data/flame.yaml"
# f.show()

if restore:
    try:
        f.restore(output, name="multi")
    except:
        print("Couldn't restore...")
else:
    # Solve with multi-component transport properties
    f.solve(loglevel)  # don't use 'auto' on subsequent solves
    # f.show()

    f.save(output, name="multi", description="solution with multicomponent transport")

print(f"multicomponent flamespeed = {f.velocity[0]:7f} m/s")

Yin = f.Y[:,0]
Yb = f.Y[:,-1]

for i,spec in enumerate(gas.species_names):
    if Yin[i] > 1.0e-10: 
        print("Inflow {:10}    :    {:.7e}".format(spec, Yin[i]))
print()
for i,spec in enumerate(gas.species_names):
    if Yb[i] > 1.0e-10: 
        print("Burnt {:10}    :    {:.7e}".format(spec, Yb[i]))


#compute flame thickness
z= f.flame.grid
T = f.T
size = len(z)-1
grad = np.zeros(size)
for i in range(size):
  grad[i] = (T[i+1]-T[i])/(z[i+1]-z[i])
thickness = (max(T) -min(T)) / max(grad)
print('laminar flame thickness = ', thickness)

files = os.listdir('ensight/vdjet/extract/')

files = [val for val in files if "extract" in val]
files.sort()


for i,file in enumerate(files):
    df = pd.read_csv('ensight/vdjet/extract/{}'.format(file))
    df = df[df['T'].notnull()]

    fig, axs = plt.subplots(3, 1, figsize=(10, 12))

    ax = axs[0]
    ax.plot(
        f.grid * 1000,
        f.T,
        lw=lw,
        label="Cantera",
        color=c[0],
    )

    ax.plot(
        df["Points:0"] * 1000,
        df["T"],
        lw=lw/2,
        ls='--',
        label="NGA2",
        color=c[0],
    )

    ax.set(
        ylabel="Temperature [K]",
        ylim=[200, 2600],
    )
    ax.legend(frameon=False, loc="upper left")
    ax.grid()

    ax = axs[1]

    ax.plot(
        f.grid * 1000,
        f.Y[gas.species_index('NXC12H26'), :],
        lw=lw,
        label="NC12H26",
        color=c[1],
    )
    ax.plot(
        f.grid * 1000,
        f.Y[gas.species_index('H2O'), :],
        lw=lw,
        label="H2O",
        color=c[2],
    )
    ax.plot(
        f.grid * 1000,
        f.Y[gas.species_index('CO'), :],
        lw=lw,
        label="CO",
        color=c[3],
    )

    ax.plot(
        f.grid * 1000,
        f.Y[gas.species_index('OH'), :] * 20.0,
        lw=lw,
        label="OH x20",
        color=c[4],
    )


    ax.plot(
        df["Points:0"] * 1000,
        df["YNC12H26"],
        lw=lw/2,
        ls='--',
        color=c[1],
    )
    ax.plot(
        df["Points:0"] * 1000,
        df["YH2O"],
        lw=lw/2,
        ls='--',
        # label="H2O",
        color=c[2],
    )
    ax.plot(
        df["Points:0"] * 1000,
        df["YCO"],
        lw=lw/2,
        ls='--',
        # label="CO",
        color=c[3],
    )

    ax.plot(
        df["Points:0"] * 1000,
        df["YOH"] * 20.0,
        lw=lw/2,
        ls='--',
        color=c[4],
    )

    ax.set(
        ylabel="Mass Fractions [-]",
        ylim=[-0.02, 0.14],
    )
    ax.grid()
    ax.legend(frameon=False, loc="upper center", ncols=5)

    ax = axs[2]
    ax.plot(
        f.grid * 1000,
        f.velocity,
        lw=lw,
        label="Velocity",
        color=c[0],
    )

    ax.plot(
        df["Points:0"] * 1000,
        df["velocity:0"],
        lw=lw/2,
        ls='--',
        # label="Temperature",
        color=c[0],
    )

    ax.set(
        ylabel="Velocity [m/s]",
        ylim=[-0.01, 2.0],
    )
    ax.legend(frameon=False, loc="upper left")
    ax.grid()


    plt.tight_layout()
    plt.savefig("figs/test.{:03}.png".format(i))
    plt.close()