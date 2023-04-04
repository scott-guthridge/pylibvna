#!/usr/bin/python3
import numpy as np
import vna.conv as vc

#
# Find the Z parameters for the 75-ohm 50-ohm impedance maching L-pad
#
z1 = 75
z2 = 50
r1 = np.sqrt(z1 * (z1 - z2))
r2 = z2 * np.sqrt(z1 / (z1 - z2))
z = [[r1+r2, r2], [r2, r2]]

#
# Find the S parameters
#
s = vc.ztos(z, [z1, z2])
print(s)

#
# Find the impedances looking into each port
#
zin = vc.stozi(s, [z1, z2])
print(zin)
