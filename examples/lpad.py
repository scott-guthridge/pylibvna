#!/usr/bin/python3
import numpy as np
from libvna.conv import ztos, stozi

#
# Find the Z parameters for the 75-ohm 50-ohm impedance maching L-pad
#
z1 = 75
z2 = 50
r1 = np.sqrt(z1 * (z1 - z2))
r2 = z2 * np.sqrt(z1 / (z1 - z2))
z = [[r1+r2, r2], [r2, r2]]
print("z:\n", np.asarray(z))

#
# Find the S parameters
#
s = ztos(z, [z1, z2])
print("\ns:\n", s)

#
# Find the impedances looking into each port
#
zin = stozi(s, [z1, z2])
print("\nzin:\n", zin)
