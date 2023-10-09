#!/usr/bin/python3
# Vector Network Analyzer Library
# Copyright © 2020-2023 D Scott Guthridge <scott_guthridge@rompromity.net>
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
Example of "through", "reflect", "line" (TRL) calibration in 10-term T
and E parameters, where the reflection parameter and line parameter are
only partially known.
"""

import sim_VNA as sim
import cmath as c
from math import sqrt
import numpy as np
from libvna.cal import Calset, CalType, Parameter, UnknownParameter, Solver
from libvna.data import NPData, PType

# Misc constants:
C = 2.9979246e+08                       # speed of light in vacuum (m/s)
Z0 = 50.0                               # reference impedance (ohms)

# Calibration frequency range and number of points
C_FMIN = 1.0e+9
C_FMAX = 8.0e+9
C_FREQUENCIES = 50

# Actual DUT parameter file
ACTUAL_FILE = "JA4220-AL.s2p"

# Calculate length of the line standard in meters:
FC = ((C_FMIN + C_FMAX) / 2.0)          # center frequency
ER_EFF = 8.25                           # effective permittivity of line
VF = (1.0 / sqrt(ER_EFF))               # velocity factor
LINE_LENGTH = (0.25 * C / FC * VF)      # line length in meters


def ideal_gamma(f):
    """
    Return the propagation constant for the ideal line.
    """
    return 2.0j * np.pi * f / (C * VF)


def make_calibration():
    """
    Solve for VNA error terms based on measurements of calibration
    standards and save the calibration to a file.
    """
    # Create the simulated VNA
    vna = sim.VNA(C_FMIN, C_FMAX, C_FREQUENCIES)

    # Create a new container to hold the calibration.
    calset = Calset()

    # Make the through measurement.  With the returned frequency vector,
    # create the error term solver, then add the measurement.
    (f_vector, a, b) = vna.measure(sim.Measurement.THROUGH)
    solver = Solver(calset, CalType.TE10, 2, 2, f_vector)
    solver.add_through(a, b, 1, 2)

    # Create the unknown reflect parameter with initial value -1 (short).
    unknown_reflect = UnknownParameter(calset, -1)

    # Make the reflect measurement.
    (_, a, b) = vna.measure(sim.Measurement.IMPERFECT_REFLECT)
    solver.add_double_reflect(a, b, unknown_reflect, unknown_reflect, 1, 2)

    # Find the ideal transmission coefficients of the line and form them
    # into a vector parameter.  Make the unknown line parameter using
    # the ideal as initial value.
    l_ideal = np.empty((C_FREQUENCIES,), dtype=np.complex128)
    for findex, f in enumerate(f_vector):
        gl = ideal_gamma(f) * LINE_LENGTH
        l_ideal[findex] = c.exp(-gl)
    unknown_line = UnknownParameter(calset, (f_vector, l_ideal))

    # Make the line measurement and add.
    (_, a, b) = vna.measure(sim.Measurement.IMPERFECT_LINE)
    s = [[0.0, unknown_line], [unknown_line, 0.0]]
    solver.add_line(a, b, s, 1, 2)

    # Solve for the error terms.  Add them back to the Calset and save
    # the calibration to a file.
    solver.solve()
    solver.add_to_calset("TE10")
    calset.save("TRL.vnacal")

    # Plot the solved unknown parameters.
    sim.plot_RL((C_FMIN, C_FMAX), f_vector, unknown_reflect, l_ideal,
                unknown_line)


def apply_calibration():
    """
    Apply the calibration we made above to measured data to find the
    true S parameters of the device under test.
    """
    # Set up the simulated VNA with the device under test
    vna = sim.VNA(dut_file="JA4220-AL.s2p")

    # Load the calibration set and find the calibration named "TE10"
    calset = Calset("TRL.vnacal")
    calibration = calset.calibrations["TE10"]

    # Measure the device under test (DUT) with the flawed VNA
    (f_vector, a, b) = vna.measure(sim.Measurement.DUT)

    # Apply the calibration correction
    result = calibration.apply(f_vector, a, b)

    # Plot results
    sim.plot_correction((1.0e+9, 8.0e+9), f_vector, np.swapaxes(b, 0, 2),
                        vna.dut.data_array, result.data_array)


make_calibration()
apply_calibration()
