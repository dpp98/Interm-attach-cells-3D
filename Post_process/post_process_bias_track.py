#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Jan 27 20:34:17 2026

@author: devi
"""

import math
import numpy as np
import os
import matplotlib.pyplot as plt


# --------------------------------------------------
# Input files
# --------------------------------------------------

filename1 = "particles.xyz"
filehandle1 = open(filename1, 'r')


# --------------------------------------------------
# Simulation parameters
# --------------------------------------------------

Ncells = 1470
Ntotal = 7000 + Ncells
N = 1001
Lx = 1000.0
kd = 0.1

xp = np.zeros(Ncells, dtype=float)
zp = np.zeros(Ncells, dtype=float)
bias = np.zeros(Ncells, dtype=float)

recoc = np.zeros(N, dtype=float)

c0 = 100.0

xgrid = np.linspace(0.0, Lx, N)


# --------------------------------------------------
# Skip initial frames
# --------------------------------------------------

skip = 3000

for i in range(skip):

    filehandle1.readline()
    filehandle1.readline()

    for j in range(Ntotal):

        s = filehandle1.readline()
        split = s.split()

        typ = int(split[0])
        idx = int(split[1])

        if typ == 1:
            xp[idx] = float(split[2])
            zp[idx] = float(split[4])
            bias[idx] = float(split[8])


# --------------------------------------------------
# Number of frames and time step
# --------------------------------------------------

npar = 1000

L = np.zeros(npar, dtype=float)
time = np.zeros(npar, dtype=float)

# Bias histories for two selected cells
bias_cell1 = np.zeros(npar, dtype=float)
bias_cell2 = np.zeros(npar, dtype=float)
bias_min = np.zeros(npar, dtype=float)
bias_max = np.zeros(npar, dtype=float)

# --------------------------------------------------
# Initial aggregate size
# --------------------------------------------------

xp_pos = xp[zp > 0.5]
L0 = xp_pos.max() - xp_pos.min()


# Time between frames
dt = 100 / 3600       # hours


# --------------------------------------------------
# Randomly select two cells
# --------------------------------------------------

rng = np.random.default_rng()

selected_cells = rng.choice(Ncells, size=2, replace=False)

cell1 = selected_cells[0]
cell2 = selected_cells[1]

print("Selected cell IDs:")
print("Cell 1:", cell1)
print("Cell 2:", cell2)


# --------------------------------------------------
# Main loop
# --------------------------------------------------

for i in range(npar):

    filehandle1.readline()
    filehandle1.readline()

    for j in range(Ntotal):

        s = filehandle1.readline()
        split = s.split()

        typ = int(split[0])
        idx = int(split[1])

        if typ == 1:
            xp[idx] = float(split[2])
            zp[idx] = float(split[4])
            bias[idx] = float(split[8])


    # ----------------------------------------------
    # Aggregate size
    # ----------------------------------------------

    xp_pos = xp[zp > 0.5]

    L[i] = (
        (xp_pos.max() - xp_pos.min()) / L0
        if xp_pos.size
        else 0.0
    )


    aggregate_mask = zp > 0.5

    if np.any(aggregate_mask):
        bias_min[i] = np.min(bias[aggregate_mask])
        bias_max[i] = np.max(bias[aggregate_mask])

    # ----------------------------------------------
    # Store bias of the two selected cells
    # ----------------------------------------------

    bias_cell1[i] = bias[cell1]
    bias_cell2[i] = bias[cell2]


    # ----------------------------------------------
    # Time
    # ----------------------------------------------

    # We skipped 'skip' frames before starting
    time[i] = dt * (i)


# --------------------------------------------------
# Close file
# --------------------------------------------------

filehandle1.close()


# --------------------------------------------------
# Save aggregate size
# --------------------------------------------------

folder_name = os.path.basename(os.getcwd())

np.save(folder_name + ".npy", L)


global_bias_min = np.min(bias_min)
global_bias_max = np.max(bias_max)

print("Minimum bias in aggregate:", global_bias_min)
print("Maximum bias in aggregate:", global_bias_max)

# --------------------------------------------------
# Plot bias of the two randomly selected cells
# --------------------------------------------------

# --------------------------------------------------
# Plot bias of the two randomly selected cells
# --------------------------------------------------

plt.figure(figsize=(6.0, 4.0))

# Cell 1 -- dark navy
plt.plot(
    time,
    bias_cell1,
    color='#1B2A49',
    linewidth=1.2,
    label=f"Cell {cell1}"
)

# Cell 2 -- brown
plt.plot(
    time,
    bias_cell2,
    color='#994C00',
    linewidth=1.2,
    label=f"Cell {cell2}"
)

# Minimum and maximum bias in the aggregate
plt.axhline(
    global_bias_min,
    color='black',
    linestyle=':',
    linewidth=1.2
)

plt.axhline(
    global_bias_max,
    color='black',
    linestyle=':',
    linewidth=1.2
)


# Axis labels
plt.xlabel("Time (hours)", fontsize=13)
plt.ylabel("Bias in cell polarity", fontsize=13)

# Axis limits
#plt.xlim(0, 11)
#plt.ylim(0.49, 0.67)


# Tick formatting
plt.tick_params(
    axis='both',
    which='major',
    labelsize=11
)


plt.tight_layout()
plt.savefig('bias_track.svg')
plt.show()
