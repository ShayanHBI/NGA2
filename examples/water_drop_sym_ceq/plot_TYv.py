# Import libraries
import os
import re
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Date file reader
def file_to_data(file, data_names):

    # Read lines
    with open(file, 'r') as sf:
        lines = sf.readlines()

    # Store header
    names = lines[0].split()

    # Remove header line
    lines.pop(0)

    fdata = []
    for line in lines:
        split_line = line.split()
        fdata.append([float(val) for val in split_line])
    fdata = np.array(fdata)

    data = {}
    for name in data_names:
        data[name] = fdata[:, names.index(name)]
    
    return data

# Use latex font
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Path where the files are located
data_dir = './monitor'

# Get and sort data
arr_names = np.array(['r', 'T', 'Yv'])
data = []
time = []

# Matches a float, with or without scientific notation
pattern = r'([-+]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)'

for filename in os.listdir(data_dir):
    if filename.startswith('T_Yv_'):
        match = re.search(pattern, filename)
        if match:
            time.append(float(filename[5:]))
            data.append(file_to_data(data_dir + '/' + filename, arr_names))

sorted_lists = sorted(zip(time, data))
time, data = zip(*sorted_lists)

# Make colors
colors = plt.cm.plasma(np.linspace(0, 1, len(time)))

fig, ax = plt.subplots(1, 1, figsize=(8,6))

# Visualize T
for i, t in enumerate(time):
    ax.plot(data[i]['r'], data[i]['T'], '-',  color=colors[i], label = r'$t = {:.3e}\ \mathrm{{s}}$'.format(t))

# Custom legend
# custom_lines = [
#     Line2D([0], [0], color='k', lw=1.5, ls='-',  label=r'$Exact$'),
#     Line2D([0], [0], color='k', lw=1.5, ls='--', label=r'$Numerical$'),
# ]
# first_legend = ax.legend(custom_lines, [r'$Exact$', r'$Numerical$'], frameon=False, loc='lower right', bbox_to_anchor=(1, 0.3), fontsize=14)
# ax.add_artist(first_legend)
ax.legend(frameon=False, loc='lower right', fontsize=14)
plt.grid(which='major', axis='both', color='gray', linestyle='-', linewidth=0.5, alpha=0.25)
plt.xlabel(r'$r$', fontsize=14)
plt.ylabel(r'$T$', fontsize=14)
for spine in plt.gca().spines.values():
    spine.set_linewidth(1.2)
plt.tight_layout()
plt.savefig('./T.pdf')

fig, ax = plt.subplots(1, 1, figsize=(8,6))

# Visualize Yv
for i, t in enumerate(time):
    ax.plot(data[i]['r'], data[i]['Yv'], '-',  color=colors[i], label = r'$t = {:.3e}\ \mathrm{{s}}$'.format(t))

# Custom legend
# custom_lines = [
#     Line2D([0], [0], color='k', lw=1.5, ls='-',  label=r'$Exact$'),
#     Line2D([0], [0], color='k', lw=1.5, ls='--', label=r'$Numerical$'),
# ]
# first_legend = ax.legend(custom_lines, [r'$Exact$', r'$Numerical$'], frameon=False, loc='lower right', bbox_to_anchor=(1, 0.3), fontsize=14)
# ax.add_artist(first_legend)
ax.legend(frameon=False, loc='lower right', fontsize=14)
plt.grid(which='major', axis='both', color='gray', linestyle='-', linewidth=0.5, alpha=0.25)
plt.xlabel(r'$r$', fontsize=14)
plt.ylabel(r'$Y_{H_2O(G)}$', fontsize=14)
for spine in plt.gca().spines.values():
    spine.set_linewidth(1.2)
plt.tight_layout()
plt.savefig('./Yv.pdf')
