#!/usr/bin/python3
#
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
#
"""
Example of "through", "reflect", "line" (TRL) calibration in 10-term T
and E parameters, where the reflection parameter and line parameter are
only partially known.

Most of the complexity in this example is in simulating the flawed VNA
and standards.  It's not until we get to the make_calibration() function
that we start to see what a user of the library needs to do.
"""

from cmath import exp, cosh, sinh
from enum import Enum
from math import sqrt
from matplotlib import ticker
from matplotlib.pyplot import show, subplots
import numpy as np
from libvna.cal import CalSet, CalType, Parameter, Solver
from libvna.conv import atos, stoa, ztos
from libvna.data import Data, PType
from warnings import warn

# Misc constants:
C = 2.9979246e+08		        # speed of light in vacuum (m/s)
Z0 = 50.0                       # reference impedance (ohms)
ER_EFF = 8.25                   # effective permittivity of line dielectric
NP_PER_DB = 0.11512925          # neper per dB
MM_PER_M = 1000.0               # mm per meter

# Calibration frequency range and number of points
C_FMIN = 1.0e+9
C_FMAX = 8.0e+9
C_FREQUENCIES = 50

# Actual DUT parameter file
ACTUAL_FILE = "JA4220-AL.s2p"

# VNA Port 1 parasitic elements
#   To insert flaws in VNA port 1, from the directional coupler facing
#   toward the DUT, we add L1 and R1 are in series and C1 shunted across
#   the port.

R1 = 10.0                               # ohms
L1 = 3.979e-9                           # henries
C1 = 1.592e-12                          # farads

# VNA Port 2 parasitic elements
#   To insert flaws in VNA port 2, from the directional coupler facing
#   toward the DUT, we add L2 and C2 in series with R2 shunted across
#   the port.
R2 = 100.0                              # ohms
L2 = 1.326e-9                           # henries
C2 = 530.5e-15                          # farads

# Errors in the Reflect Standard
#    The standard should be a short, but our imperfect standard
#    has resistor RR in series with inductor, RL.
RR = 5.0                                # ohms
RL = 707.4e-12                          # henries

# Errors in the Line Standard
#    The line standard should be a perfect transmission line, but
#    ours has both loss and a phase error.
LINE_LOSS = 0.5                         # dB/mm
PHASE_ERROR = 10.                       # degrees

# Calculate length of the line standard in meters:
FC = ((C_FMIN + C_FMAX) / 2.0)          # center frequency
VF = (1.0 / sqrt(ER_EFF))               # velocity factor
LINE_LENGTH = (0.25 * C / FC * VF)      # line length in meters


def ideal_gamma(f):
    """
    Return the propagation constant for the ideal line.
    """
    return 2.0j * np.pi * (f) / (C * VF)


def actual_gamma(f):
    """
    Return the propagation constant of the actual line.
    """
    return ideal_gamma(f) * exp(1.0j * np.pi * PHASE_ERROR / 180.0) \
        + NP_PER_DB * MM_PER_M * LINE_LOSS


class Measurement(Enum):
    """
    Select which measurement the simulated VNA should return.
    """
    THROUGH = 0
    REFLECT = 1
    LINE = 2
    DUT = 3


class VNA:
    """
    Simulated VNA
    """
    def __init__(self, dut_file=None):
        if dut_file is not None:
            self.dut = Data()
            self.dut.load(dut_file)
            self.dut.convert(PType.S)
        else:
            self.dut = None

    def measure(self, measurement: Measurement):
        """
        Simulate a measurement of a standard or the DUT.

        Parameters
        ----------
        measurement:
            Select what to measure: THROUGH, REFLECT, LINE or DUT

        f_vector:
            Vector of frequency points.  Use for THROUGH, REFLECT and
            LINE only.

        Returns
        -------
        f_vector:
            vector of frequencies where measurements were made

        a_matrix:
            (2 x 2 x frequencies) measured reference matrix

        b_matrix:
            (2 x 2 x frequencies) measured value matrix
        """
        # Get frequency vector.
        if measurement == Measurement.DUT:
            if self.dut is None:
                raise ValueError("measure: for DUT measurement, dut_file "
                                 "must be given to constructor")
            f_vector = self.dut.frequency_vector
        else:
            f_vector = np.linspace(C_FMIN, C_FMAX, C_FREQUENCIES)

        # Allocate result matrices
        a_matrix = np.empty((2, 2, len(f_vector)), dtype=np.complex128)
        b_matrix = np.empty((2, 2, len(f_vector)), dtype=np.complex128)

        # For all measurements, fill a to simulate leakage in the VNA
        # switch.  Send 2/3 of the signal to the intended port and 1/3
        # to the other.
        a = np.asarray([[2.0 / 3.0, 1.0 / 3.0], [1.0 / 3.0, 2.0 / 3.0]],
                       dtype=np.complex128)

        # For each frequency...
        for findex, f in enumerate(f_vector):
            s = 2.0j * np.pi * f

            # Find ABCD parameters for the errors at VNA port 1 using
            # temporary matrices u, v and w, and storing the result into
            # port1_abcd.  Detector is on the left; DUT is on the right.

            # series inductor L1
            u = np.asarray([[1.0, L1 * s], [0.0, 1.0]], dtype=np.complex128)

            # series resistor R1
            v = np.asarray([[1.0, R1], [0.0, 1.0]], dtype=np.complex128)
            w = np.matmul(u, v)

            # shunt capacitor C1
            u = np.asarray([[1.0, 0.0], [C1 * s, 1.0]], dtype=np.complex128)
            port1_abcd = np.matmul(w, u)

            # Find ABCD parameters for the errors at VNA port 2, storing
            # the result into port2_abcd.  DUT is on the left; detector
            # is on the right.

            # shunt resistor R2
            u = np.asarray([[1.0, 0.0], [1.0 / R2, 1.0]], dtype=np.complex128)

            # series inductor L2
            v = np.asarray([[1.0, L2 * s], [0.0, 1.0]], dtype=np.complex128)
            w = np.matmul(u, v)

            # series capacitor C2
            u = np.asarray([[1.0, 1.0 / (C2 * s)], [0.0, 1.0]],
                           dtype=np.complex128)
            port2_abcd = np.matmul(w, u)

            # Calculate the b matrix for the requested measurement
            if measurement == Measurement.THROUGH:
                # Multiply the ABCD parameters of the two error boxes,
                # convert to s-parameters and find b = s a.
                u = np.matmul(port1_abcd, port2_abcd)
                v = atos(u, Z0)
                b = np.matmul(v, a)

            elif measurement == Measurement.REFLECT:
                # We have VNA port 1 on the left side of the port1_abcd
                # error box with reflect standard on the right.  And we
                # have VNA port 2 on the right side of the port2_abcd
                # error box with the reflected standard on the left.
                # For port 1, start with the definition of ABCD
                # parameters:
                #
                #     [ v1 ]   [ A11 A12 ] [  v2 ]
                #     [    ] = [         ] [     ]
                #     [ i1 ]   [ A21 A22 ] [ -i2 ]
                #
                # with zr on the right, set v2 == -i2 zr.  The minus sign
                # is needed because i2 is defined as current from the
                # reflect standard into the error box.  The impedance
                # the VNA sees looking into the left side of the error
                # box is: i1 / i1.  This simplifies to:
                #
                #                 P1_A12 + P1_A11 zr
                #    Zin_left =  --------------------
                #                 P1_A22 + P1_A21 zr
                #
                # For port 2, we can multiply each side of the ABCD
                # equation above on the left with A^-1, set v1 = -i1 zr,
                # and find v2 / i2.  This simplifies to:
                #
                #                 P2_A12 + P2_A22 zr
                #    Zin_right = --------------------
                #                 P2_A11 + P2_A21 zr
                #
                # From these, we can construct the Z parameters of the
                # two error boxes with double reflect in the middle.
                # No signal passes through the standard, so the
                # off-diagonal entries are zero.
                #
                #        [ Zin_left  0         ]
                #    Z = [                     ]
                #        [ 0         Zin_right ]
                zr = RR + RL * s
                z = np.asarray([[(port1_abcd[0, 1] + port1_abcd[0, 0] * zr) /
                                 (port1_abcd[1, 1] + port1_abcd[1, 0] * zr),
                                 0.0],
                                [0.0,
                                 (port2_abcd[0, 1] + port2_abcd[1, 1] * zr) /
                                 (port2_abcd[0, 0] + port2_abcd[1, 0] * zr)]],
                               dtype=np.complex128)
                v = ztos(z, Z0)
                b = np.matmul(v, a)

            elif measurement == Measurement.LINE:
                # Multiply the ABCD parameters of the first error box,
                # the line and the second error box. Next, convert to
                # s-parameters and find b = s a.
                gl = actual_gamma(f) * LINE_LENGTH
                u = np.asarray([[cosh(gl), sinh(gl) * Z0],
                                [sinh(gl) / Z0, cosh(gl)]],
                               dtype=np.complex128)
                v = np.matmul(port1_abcd, u)
                u = np.matmul(v, port2_abcd)
                v = atos(u, Z0)
                b = np.matmul(v, a)

            elif measurement == Measurement.DUT:
                # Convert the actual s-parameters of the DUT to ABCD
                # parameters.  Multiply the ABCD parameters of the first
                # error box, the DUT and the second error box.  Finally,
                # convert to s-parameters and find b = s a.
                u = self.dut.data_array[findex, ...]
                v = stoa(u, Z0)
                w = np.matmul(port1_abcd, v)
                u = np.matmul(w, port2_abcd)
                v = atos(u, Z0)
                b = np.matmul(v, a)

            else:
                raise AssertionError("VNA.measure: invalid measurement")

            # Copy a and b to the return matrices.
            a_matrix[..., findex] = a
            b_matrix[..., findex] = b

        return (f_vector, a_matrix, b_matrix)


def make_calibration():
    """
    Solve for VNA error terms based on measurements of calibration
    standards and save the calibration to a file.
    """

    # Create a new container to hold the calibration.
    calset = CalSet.create()

    # Create the simulated VNA
    vna = VNA()

    # Make the through measurement.  With the returned frequency vector,
    # create the error term solver then add the measurement.
    (f_vector, a, b) = vna.measure(Measurement.THROUGH)
    solver = calset.make_solver(CalType.TE10, f_vector, 2, 2)
    solver.add_through(a, b, 1, 2)

    # Create the unknown reflect parameter with initial value -1 (short).
    unknown_reflect = calset.make_unknown(-1)

    # Make the reflect measurement.
    (_, a, b) = vna.measure(Measurement.REFLECT)
    solver.add_double_reflect(a, b, unknown_reflect, unknown_reflect, 1, 2)

    # Find the ideal transmission coefficients of the line and form
    # them into a vector parameter.  Make the unknown line parameter
    # with ideal as the initial value.
    l_ideal = np.empty((C_FREQUENCIES,), dtype=np.complex128)
    for findex, f in enumerate(f_vector):
        gl = ideal_gamma(f) * LINE_LENGTH
        l_ideal[findex] = exp(-gl)
    unknown_line = calset.make_unknown((f_vector, l_ideal))

    # Make the line measurement and add.
    (_, a, b) = vna.measure(Measurement.LINE)
    s = [[0.0, unknown_line], [unknown_line, 0.0]]
    solver.add_line(a, b, s, 1, 2)

    # Solve for the error terms.  Add back to the CalSet and save
    # the calibration to a file.
    solver.solve()
    calset.add(solver, "TE10")
    calset.save("TRL.vnacal")

    # Plot the solved unknown parameters.
    make_RL_plots(f_vector, unknown_reflect, l_ideal, unknown_line)


def apply_calibration():
    """
    Apply the calibration we made above to measured data to find the
    true S parameters of the device under test.
    """
    # Load the calibration set and find the calibration named "TE10"
    calset = CalSet.load("TRL.vnacal")
    calibration = calset.calibrations["TE10"]

    # Set up the simulated VNA with a device under test
    vna = VNA(dut_file="JA4220-AL.s2p")

    # Measure the device with the flawed VNA
    (f_vector, a, b) = vna.measure(Measurement.DUT)

    # Apply the calibration correction
    result = calibration.apply(f_vector, a, b)

    # Plot results
    plot_DUT(f_vector, np.swapaxes(b, 0, 2),
             vna.dut.data_array, result.data_array)


def make_RL_plots(f_vector, unknown_reflect, l_ideal, unknown_line):
    """
    Plot ideal, actual and solved unknown parameter values.
    """
    # Find the actual and corrected reflection coefficient
    r_actual = np.empty((len(f_vector),), dtype=np.complex128)
    r_solved = np.empty((len(f_vector),), dtype=np.complex128)
    for findex, f in enumerate(f_vector):
        s = 2.0j * np.pi * f
        zr = RR + RL * s
        r_actual[findex] = (zr - Z0) / (zr + Z0)
        r_solved[findex] = unknown_reflect.get_value(f)

    # Reflect
    fig, axs = subplots(2, 1)
    axs[0].plot(f_vector, -np.ones((C_FREQUENCIES,)),
                label="$reflect_r$ ideal",
                color="C0", linestyle=":")
    axs[0].plot(f_vector, np.zeros((C_FREQUENCIES,)),
                label="$reflect_i$ ideal",
                color="C1", linestyle=":")
    axs[0].plot(f_vector, r_actual.real, 'o',
                label="$reflect_r$ actual",
                color="C0", markersize=2)
    axs[0].plot(f_vector, r_actual.imag, 's',
                label="$reflect_i$ actual",
                color="C1", markersize=2)
    axs[0].plot(f_vector, r_solved.real,
                label="$reflect_r$ solved",
                color="C0")
    axs[0].plot(f_vector, r_solved.imag,
                label="$reflect_i$ solved",
                color="C1")
    axs[0].set_title('Unknown Reflection Coefficient')
    axs[0].set_xlabel('frequency (Hz)')
    axs[0].set_ylabel('reflection coefficient')
    axs[0].set_xlim([C_FMIN, C_FMAX])
    axs[0].set_ylim([-1.0, +1.0])
    axs[0].grid()
    axs[0].legend()

    # Find the actual and corrected transmission coefficient
    l_actual = np.empty((len(f_vector),), dtype=np.complex128)
    l_solved = np.empty((len(f_vector),), dtype=np.complex128)
    for findex, f in enumerate(f_vector):
        l_actual[findex] = exp(-actual_gamma(f) * LINE_LENGTH)
        l_solved[findex] = unknown_line.get_value(f)

    # Line
    axs[1].plot(f_vector, l_ideal.real,
                label="$line_r$ ideal",
                color="C2", linestyle=":")
    axs[1].plot(f_vector, l_ideal.imag,
                label="$line_i$ ideal",
                color="C3", linestyle=":")
    axs[1].plot(f_vector, l_actual.real, 'o',
                label="$line_r$ actual",
                color="C2", markersize=2)
    axs[1].plot(f_vector, l_actual.imag, 's',
                label="$line_i$ actual",
                color="C3", markersize=2)
    axs[1].plot(f_vector, l_solved.real,
                label="$line_r$ solved",
                color="C2")
    axs[1].plot(f_vector, l_solved.imag,
                label="$line_i$ solved",
                color="C3")
    axs[1].set_title('Unknown Transmission Coefficient')
    axs[1].set_xlabel('frequency (Hz)')
    axs[1].set_ylabel('transmission coefficient')
    axs[1].set_xlim([C_FMIN, C_FMAX])
    axs[1].set_ylim([-1.0, +1.0])
    axs[1].grid()
    axs[1].legend()

    show()


def plot_DUT(f_vector, measured, actual, corrected):
    #
    # S11
    #
    fig, axs = subplots(2, 2)
    axs[0, 0].plot(f_vector, measured[:, 0, 0].real,
                   label="$s11_r$ measured",
                   color="C0", linestyle=":")
    axs[0, 0].plot(f_vector, measured[:, 0, 0].imag,
                   label="$s11_i$ measured",
                   color="C1", linestyle=":")
    axs[0, 0].plot(f_vector, actual[:, 0, 0].real, 'o',
                   label="$s11_r$ actual",
                   color="C0", markersize=2)
    axs[0, 0].plot(f_vector, actual[:, 0, 0].imag, 's',
                   label="$s11_i$ actual",
                   color="C1", markersize=2)
    axs[0, 0].plot(f_vector, corrected[:, 0, 0].real,
                   label="$s11_r$ corrected",
                   color="C0")
    axs[0, 0].plot(f_vector, corrected[:, 0, 0].imag,
                   label="$s11_i$ corrected",
                   color="C1")
    axs[0, 0].set_title('S11')
    axs[0, 0].set_xlabel('frequency (Hz)')
    axs[0, 0].set_ylabel('reflection')
    axs[0, 0].set_xlim([f_vector[0], f_vector[-1]])
    axs[0, 0].set_ylim([-1.0, +1.0])
    axs[0, 0].grid()
    axs[0, 0].legend()

    #
    # S12
    #
    axs[0, 1].plot(f_vector, measured[:, 0, 1].real,
                   label="$s12_r$ measured",
                   color="C2", linestyle=":")
    axs[0, 1].plot(f_vector, measured[:, 0, 1].imag,
                   label="$s12_i$ measured",
                   color="C3", linestyle=":")
    axs[0, 1].plot(f_vector, actual[:, 0, 1].real, 'o',
                   label="$s12_r$ actual",
                   color="C2", markersize=2)
    axs[0, 1].plot(f_vector, actual[:, 0, 1].imag, 's',
                   label="$s12_i$ actual",
                   color="C3", markersize=2)
    axs[0, 1].plot(f_vector, corrected[:, 0, 1].real,
                   label="$s12_r$ corrected",
                   color="C2")
    axs[0, 1].plot(f_vector, corrected[:, 0, 1].imag,
                   label="$s12_i$ corrected",
                   color="C3")
    axs[0, 1].set_title('S12')
    axs[0, 1].set_xlabel('frequency (Hz)')
    axs[0, 1].set_ylabel('transmission')
    axs[0, 1].set_xlim([f_vector[0], f_vector[-1]])
    axs[0, 1].set_ylim([-1.0, +1.0])
    axs[0, 1].grid()
    axs[0, 1].legend()

    #
    # S21
    #
    axs[1, 0].plot(f_vector, measured[:, 1, 0].real,
                   label="$s21_r$ measured",
                   color="C4", linestyle=":")
    axs[1, 0].plot(f_vector, measured[:, 1, 0].imag,
                   label="$s21_i$ measured",
                   color="C5", linestyle=":")
    axs[1, 0].plot(f_vector, actual[:, 1, 0].real, 'o',
                   label="$s21_r$ actual",
                   color="C4", markersize=2)
    axs[1, 0].plot(f_vector, actual[:, 1, 0].imag, 's',
                   label="$s21_i$ actual",
                   color="C5", markersize=2)
    axs[1, 0].plot(f_vector, corrected[:, 1, 0].real,
                   label="$s21_r$ corrected",
                   color="C4")
    axs[1, 0].plot(f_vector, corrected[:, 1, 0].imag,
                   label="$s21_i$ corrected",
                   color="C5")
    axs[1, 0].set_title('S21')
    axs[1, 0].set_xlabel('frequency (Hz)')
    axs[1, 0].set_ylabel('transmission')
    axs[1, 0].set_xlim([f_vector[0], f_vector[-1]])
    axs[1, 0].set_ylim([-1.0, +1.0])
    axs[1, 0].grid()
    axs[1, 0].legend()

    #
    # S22
    #
    axs[1, 1].plot(f_vector, measured[:, 1, 1].real,
                   label="$s22_r$ measured",
                   color="C6", linestyle=":")
    axs[1, 1].plot(f_vector, measured[:, 1, 1].imag,
                   label="$s22_i$ measured",
                   color="C7", linestyle=":")
    axs[1, 1].plot(f_vector, actual[:, 1, 1].real, 'o',
                   label="$s22_r$ actual",
                   color="C6", markersize=2)
    axs[1, 1].plot(f_vector, actual[:, 1, 1].imag, 's',
                   label="$s22_i$ actual",
                   color="C7", markersize=2)
    axs[1, 1].plot(f_vector, corrected[:, 1, 1].real,
                   label="$s22_r$ corrected",
                   color="C6")
    axs[1, 1].plot(f_vector, corrected[:, 1, 1].imag,
                   label="$s22_i$ corrected",
                   color="C7")
    axs[1, 1].set_title('S22')
    axs[1, 1].set_xlabel('frequency (Hz)')
    axs[1, 1].set_ylabel('reflection')
    axs[1, 1].set_xlim([f_vector[0], f_vector[-1]])
    axs[1, 1].set_ylim([-1.0, +1.0])
    axs[1, 1].grid()
    axs[1, 1].legend()

    show()


make_calibration()
apply_calibration()
