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
data_dir = "/Users/shayanhbi/Repositories/NGA2/examples/stephan/temperature"

# Get and sort data
arr_names = np.array(['x', 'Tg'])
data_ext = []
data_num = []
time_num = []
time_ext = []

pattern = r'[-+]?\d*\.\d+|\d+'
for filename in os.listdir(data_dir):
    if filename.endswith('.dat'):
        match = re.search(pattern, filename)
        time = float(match.group())
        if 'num' in filename:
            time_num.append(time)
            data_num.append(file_to_data(data_dir +'/' + filename, arr_names))
        elif 'ext' in filename:
            time_ext.append(time)
            data_ext.append(file_to_data(data_dir + '/' + filename, arr_names))
        else:
            pass

sorted_lists = sorted(zip(time_num, data_num))
time_num, data_num = zip(*sorted_lists)

sorted_lists = sorted(zip(time_ext, data_ext))
time_ext, data_ext = zip(*sorted_lists)


# Make colors
colors = plt.cm.plasma(np.linspace(0, 1, len(time_ext)))

fig, ax = plt.subplots(1, 1, figsize=(8,6))

# Visualize
for i, t in enumerate(time_ext):
    ax.plot(1000*data_ext[i]['x'], data_ext[i]['Tg'], '-',  linewidth=4, color=colors[i], label = r'$t = {:.3f}~(s)$'.format(t))
    ax.plot(1000*data_num[i]['x'], data_num[i]['Tg'], '--', linewidth=4, color=colors[i])

# Custom legend
custom_lines = [
    Line2D([0], [0], color='k', lw=4, ls='-',  label=r'$Analytical$'),
    Line2D([0], [0], color='k', lw=4, ls='--', label=r'$Numerical$'),
]
first_legend = ax.legend(custom_lines, [r'$Analytical$', r'$Numerical$'], frameon=False, bbox_to_anchor=(0.58, 1.0), loc='upper right', fontsize=20)
ax.add_artist(first_legend)
ax.legend(frameon=False, loc='lower right', bbox_to_anchor=(1, 0.42), fontsize=20)
plt.grid(which='major', axis='both', color='gray', linestyle='-', linewidth=0.7, alpha=0.25)
plt.xlabel(r'$x~(mm)$',  fontsize=26)
plt.ylabel(r'$T_g~(K)$', fontsize=26)
ax.tick_params(axis='both', which='major', labelsize=22)
for spine in plt.gca().spines.values():
    spine.set_linewidth(2)
plt.tight_layout()
plt.savefig('./T_g.pdf')
