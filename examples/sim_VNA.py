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
Simulate an imperfect VNA and imperfect calibration standards.
"""

import cmath as c
from enum import auto, Enum
from math import sqrt
from matplotlib import ticker
from matplotlib.pyplot import show, subplots
import numpy as np
from libvna.conv import ytos, ztos
from libvna.data import NPData, PType

# Misc constants
C = 2.9979246e+08                       # speed of light in vacuum (m/s)
Z0 = 50.0                               # reference impedance (ohms)
NP_PER_DB = 0.11512925                  # neper per dB
MM_PER_M = 1000.0                       # mm per meter
GHz = 1.0e+9                            # GHz per Hz

# Parasitic Elements: we represent the errors in our simulated VNA as
# a four-port passive network with its ports connected to VNA ports 1,
# 2 and DUT ports 1, 2, respectively.  Each port has an admittance to
# each other port and each has an admittance to ground.

# The parallel combination of R1, L1 and C1 connects VNA port1 to
# DUT port1.
R1 = 300.0              # ohms
L1 = 3.979e-9           # henries
C1 = 1.592e-12          # farads

# Parallel combination of R2 with the series combination of L2 and C2
# connect VNA port 2 to DUT port 2.
R2 = 200.0              # ohms
L2 = 1.137e-9           # henries
C2 = 454.73e-15         # farads

# C3 shunts VNA port 1 to ground
C3 = 2.52693e-12        # farads

# L3 shunts VNA port 2 to ground
L3 = 3.07487e-9         # henries

# C4 shunts DUT port 1 to ground
C4 = 423.208e-15        # farads

# L4 shunts DUT port 2 to ground
L4 = 1.93492e-9         # henries

# L5 connects VNA port 1 to VNA port 2 (16 error terms only)
L5 = 1.27154e-9         # henries

# L6 connects VNA port 1 to DUT port 2 (16 error terms only)
L6 = 1.15560e-9         # henries

# C5 connects VNA port 2 to DUT port 1 (16 error terms only)
C5 = 1.15119e-12        # farads

# C6 connects DUT port 1 to DUT port 2 (16 error terms only)
C6 = 455.068e-15        # farads

# Internal VNA Leakage Terms: these arbitrary leakage terms are used
# only for the 10 error terms calibrations.
E12 = 0.123541+0.286863j
E21 = 0.275458-0.376982j

# Arbitrary reflection coefficient for terminator used on ports
# not being directly measured.
TERMINATOR_GAMMA = 0.5+0.3j

# Errors in the Unknown Reflect Standard
#    The standard should be a short, but our imperfect standard
#    has resistor RR in series with inductor, RL.
RR = 5.0                                # ohms
RL = 707.4e-12                          # henries

# These constants are used to calculate length of the line standard
# in meters:
ER_EFF = 8.25                           # effective permittivity of line
VF = (1.0 / sqrt(ER_EFF))               # velocity factor

# Errors in the Unknown Line Standard
#    The line standard should be a perfect transmission line, but
#    ours has both loss and a phase error from the ideal values.
LINE_LOSS = 0.5                         # dB/mm
PHASE_ERROR = 10.                       # degrees


def calc_ideal_line(f):
    """
    Return the propagation constant for the ideal line.
    """
    return 2.0j * np.pi * f / (C * VF)


def calc_actual_line(f):
    """
    Return the propagation constant of the actual line.
    """
    return calc_ideal_line(f) * c.exp(1.0j * np.pi * PHASE_ERROR / 180.0) \
        + NP_PER_DB * MM_PER_M * LINE_LOSS


def calc_line_length(f_min, f_max):
    """
    Return the length of the line standard in meters.
    """
    fc = (f_min + f_max) / 2.0
    return 0.25 * C / fc * VF


class Measurement(Enum):
    """
    Select which measurement the simulated VNA should return.
    """
    THROUGH = auto()            # ideal through
    SHORT1 = auto()             # short on port 1, 100+100j on port 2
    SHORT2 = auto()             # short on port 2, 100+100j on port 1
    OPEN1 = auto()              # open  on port 1, 100+100j on port 2
    OPEN2 = auto()              # open  on port 2, 100+100j on port 1
    LOAD1 = auto()              # load  on port 1, 100+100j on port 2
    LOAD2 = auto()              # load  on port 2, 100+100j on port 1
    MATCH_MATCH = auto()        # match on port 1, match on port 2
    MATCH_OPEN = auto()         # match on port 1, open  on port 2
    MATCH_SHORT = auto()        # match on port 1, short on port 2
    OPEN_MATCH = auto()         # open  on port 1, match on port 2
    OPEN_OPEN = auto()          # open  on port 1, open  on port 2
    OPEN_SHORT = auto()         # open  on port 1, short on port 2
    SHORT_MATCH = auto()        # short on port 1, match on port 2
    SHORT_OPEN = auto()         # short on port 1, open  on port 2
    SHORT_SHORT = auto()        # short on port 1, short on port 2
    IDEAL_LINE = auto()         # ideal line
    IMPERFECT_REFLECT = auto()  # symmetrical imperfect short on both ports
    IMPERFECT_LINE = auto()     # imperfect line
    DUT = auto()                # DUT loaded from file


class VNA:
    """
    Simulated VNA

    Args:
        f_min: (float, optional) start of calibration frequency range
        f_max: (float, optional) end of calibration frequency range
        frequencies: (int, optional) number of calibration frequencies
        n_errors (int, optional): number of error terms:
            8, 10 or 16 (default 10)
        return_ab (bool, optional): if True, return the frequency vector,
            reference matrix (a) and measurement matrix (b); if False,
            return the frequency_vector and measurement matrix only.
        dut_file (str, optional): filename containing the correct network
            parameters of the DUT.  Needed for Measurement.DUT.
    """
    def __init__(self, f_min=None, f_max=None, frequencies=None, *,
                 n_errors=10, return_ab=True, dut_file=None):
        if n_errors != 8 and n_errors != 10 and n_errors != 16:
            raise ValueError("n_errors must be 8, 10 or 16")
        dut = None
        if dut_file is not None:
            dut = NPData()
            dut.load(dut_file)
            dut.convert(PType.S)

        self.f_min = f_min
        self.f_max = f_max
        self.frequencies = frequencies
        self.n_errors = n_errors
        self.return_ab = return_ab
        self.dut = dut

    def measure(self, measurement: Measurement):
        """
        Simulate a measurement of a calibration standard or device under
        test with an imperfect VNA and test fixture.

        Args:
            measurement (Measurement): select what to measure:
                THROUGH: ideal through
                SHORT1: short on port 1, 100+100j on port 2
                SHORT2: short on port 2, 100+100j on port 1
                OPEN1: open  on port 1, 100+100j on port 2
                OPEN2: open  on port 2, 100+100j on port 1
                LOAD1: load  on port 1, 100+100j on port 2
                LOAD2: load  on port 2, 100+100j on port 1
                MATCH_MATCH: match on port 1, match on port 2
                MATCH_OPEN: match on port 1, open  on port 2
                MATCH_SHORT: match on port 1, short on port 2
                OPEN_MATCH: open  on port 1, match on port 2
                OPEN_OPEN: open  on port 1, open  on port 2
                OPEN_SHORT: open  on port 1, short on port 2
                SHORT_MATCH: short on port 1, match on port 2
                SHORT_OPEN: short on port 1, open  on port 2
                SHORT_SHORT: short on port 1, short on port 2
                IDEAL_LINE: ideal line
                IMPERFECT_REFLECT: symmetrical imperfect short on both
                IMPERFECT_LINE: imperfect line
                DUT: device under test s-parameters loaded from file

        Returns:
            tuple (f_vector, [a_matrix, ] b_matrix):
                f_vector: frequency_vector
                a_matrix: reference matrix (if return_ab is True)
                b_matrix: measurement matrix
        """
        # Get frequency vector.
        if measurement == Measurement.DUT:
            if self.dut is None:
                raise ValueError("measure: for DUT measurement, dut_file "
                                 "must be given to constructor")
            f_vector = self.dut.frequency_vector
        else:
            f_vector = np.linspace(self.f_min, self.f_max, self.frequencies)

        # Allocate result matrices
        if self.return_ab:
            a_matrix = np.empty((2, 2, len(f_vector)), dtype=np.complex128)
        b_matrix = np.empty((2, 2, len(f_vector)), dtype=np.complex128)

        # If return_ab is True, create an "a" matrix for all measurements.
        # This simple examples sends 2/3 of the signal to the intended
        # port and 1/3 to the other at all frequencies.
        if self.return_ab:
            a = np.asarray([[2.0 / 3.0, 1.0 / 3.0], [1.0 / 3.0, 2.0 / 3.0]],
                           dtype=np.complex128)

        # If frequencies were given, calculate the length of the line
        # standard.
        if self.f_min is not None and self.f_max is not None:
            line_length = calc_line_length(self.f_min, self.f_max)

        # For each frequency...
        for findex, f in enumerate(f_vector):
            jω = 2.0j * np.pi * f

            # We simulate the flaws in our VNA and test fixture using
            # a four-port error box with ports 1 and 2 connected to
            # VNA ports 1 and 2, respectively; and error box ports
            # 3 and 4 connected to DUT ports 1 and 2, respectively.
            # For simplicity, we use a passive extended pi network
            # where each port of the errox box has an admittance to
            # ground and each has an admittance to each other port.
            # A consequence of this approach is that the S-parameters of
            # the error box are symmetrical, thus we really have fewer
            # error terms than we could have.

            # Conductance between VNA port 1 and DUT port 1 is a
            # parallel combination of R1, L1 and C1.
            y_v1d1 = 1 / (L1 * jω) + 1 / R1 + C1 * jω

            # Conductance between VNA port 2 and DUT port 2 is R2
            # in parallel with the series combination of L2 and C2.
            y_v2d2 = 1 / R2 + 1 / (L2 * jω + 1 / (C2 * jω))

            # VNA port 1 is shunted to ground by C3
            y_v1g = C3 * jω

            # VNA port 2 is shunted to ground by L3
            y_v2g = 1 / (L3 * jω)

            # DUT port 1 is shunted to ground by L4
            y_d1g = 1 / (L4 * jω)

            # DUT port 2 is shunted to ground by C4
            y_d2g = C4 * jω

            # If we're building errors for a 16-term error model,
            # add the cross port leakage terms.
            if self.n_errors == 16:
                # L5 bridges VNA port 1 to VNA port 2
                y_v1v2 = 1 / (L5 * jω)

                # L6 bridges VNA port 1 to DUT port 2
                y_v1d2 = 1 / (L6 * jω)

                # C5 bridges VNA port 2 to DUT port 1
                y_v2d1 = C5 * jω

                # C6 bridges DUT port 1 to DUT port 2
                y_d1d2 = C6 * jω

            else:
                y_v1v2 = 0
                y_v1d2 = 0
                y_v2d1 = 0
                y_d1d2 = 0

            # Form the Y parameters of the error box and convert to S
            s_terms = ytos([
                [y_v1g + y_v1v2 + y_v1d1 + y_v1d2, -y_v1v2, -y_v1d1, -y_v1d2],
                [-y_v1v2, y_v2g + y_v1v2 + y_v2d1 + y_v2d2, -y_v2d1, -y_v2d2],
                [-y_v1d1, -y_v2d1, y_d1g + y_v1d1 + y_v2d1 + y_d1d2, -y_d1d2],
                [-y_v1d2, -y_v2d2, -y_d1d2, y_d2g + y_v1d2 + y_v2d2 + y_d1d2]])

            # Partition s_terms into a 2x2 matrix of 2x2 matrices.
            S11 = s_terms[0:2, 0:2] # directivity terms
            S12 = s_terms[0:2, 2:4] # reflection tracking terms
            S21 = s_terms[2:4, 0:2] # transmission tracking terms
            S22 = s_terms[2:4, 2:4] # port match terms

            # Convert from S parameters to T parameters (with matrix elements)
            # to create the T error terms.
            T22 = np.linalg.inv(S21)
            T21 = -np.matmul(T22, S22)
            T12 = np.matmul(S11, T22)
            T11 = S12 - np.matmul(T12, S22)

            # Find the s matrix for each requested measurement
            if measurement == Measurement.SHORT1:
                s = [[-1, 0], [0, TERMINATOR_GAMMA]]

            elif measurement == Measurement.SHORT2:
                s = [[TERMINATOR_GAMMA, 0], [0, -1]]

            elif measurement == Measurement.OPEN1:
                s = [[1, 0], [0, TERMINATOR_GAMMA]]

            elif measurement == Measurement.OPEN2:
                s = [[TERMINATOR_GAMMA, 0], [0, 1]]

            elif measurement == Measurement.LOAD1:
                s = [[0, 0], [0, TERMINATOR_GAMMA]]

            elif measurement == Measurement.LOAD2:
                s = [[TERMINATOR_GAMMA, 0], [0, 1]]

            elif measurement == Measurement.MATCH_MATCH:
                s = [[0, 0], [0, 0]]

            elif measurement == Measurement.MATCH_OPEN:
                s = [[0, 0], [0, 1]]

            elif measurement == Measurement.MATCH_SHORT:
                s = [[0, 0], [0, -1]]

            elif measurement == Measurement.OPEN_MATCH:
                s = [[1, 0], [0, 0]]

            elif measurement == Measurement.OPEN_OPEN:
                s = [[1, 0], [0, 1]]

            elif measurement == Measurement.OPEN_SHORT:
                s = [[1, 0], [0, -1]]

            elif measurement == Measurement.SHORT_MATCH:
                s = [[-1, 0], [0, 0]]

            elif measurement == Measurement.SHORT_OPEN:
                s = [[-1, 0], [0, 1]]

            elif measurement == Measurement.SHORT_SHORT:
                s = [[-1, 0], [0, -1]]

            elif measurement == Measurement.IMPERFECT_REFLECT:
                gamma = ztos([[RR + RL * jω]], Z0)[0, 0]
                s = [[gamma, 0], [0, gamma]]

            elif measurement == Measurement.THROUGH:
                s = [[0, 1], [1, 0]]

            elif measurement == Measurement.IDEAL_LINE:
                gl = calc_ideal_line(f) * line_length
                temp = c.exp(-gl)
                s = [[0, temp], [temp, 0]]

            elif measurement == Measurement.IMPERFECT_LINE:
                gl = calc_actual_line(f) * line_length
                temp = c.exp(-gl)
                s = [[0, temp], [temp, 0]]

            elif measurement == Measurement.DUT:
                s = self.dut.data_array[findex, ...]

            # Embed the error terms into the measurement.
            m = np.matmul(np.matmul(T11, s) + T12,
                       np.linalg.inv(np.matmul(T21, s) + T22))

            # If we have leakage terms outside of the linear system,
            # add them.
            if self.n_errors == 10:
                m += [[0, E12], [E21, 0]]

            # Copy the result to the result matrix.
            if self.return_ab:
                b = np.matmul(m, a)
                a_matrix[..., findex] = a
                b_matrix[..., findex] = b
            else:
                b_matrix[..., findex] = m

        # Return results
        if self.return_ab:
            return (f_vector, a_matrix, b_matrix)
        else:
            return (f_vector, b_matrix)


def plot_RL(f_range, f_vector, unknown_reflect, l_ideal, unknown_line):
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

    # Find the length of the line standard.
    line_length = calc_line_length(*f_range)

    # Scale f_range and  f_vector to GHz.
    f_range = np.asarray(f_range) / GHz
    f_vector_GHz = f_vector / GHz

    # Reflect
    fig, axs = subplots(2, 1)
    axs[0].plot(f_vector_GHz, -np.ones((len(f_vector),)),
                label="$reflect_r$ ideal",
                color="C0", linestyle=":")
    axs[0].plot(f_vector_GHz, np.zeros((len(f_vector),)),
                label="$reflect_i$ ideal",
                color="C1", linestyle=":")
    axs[0].plot(f_vector_GHz, r_actual.real, 'o',
                label="$reflect_r$ actual",
                color="C0", markersize=2)
    axs[0].plot(f_vector_GHz, r_actual.imag, 's',
                label="$reflect_i$ actual",
                color="C1", markersize=2)
    axs[0].plot(f_vector_GHz, r_solved.real,
                label="$reflect_r$ solved",
                color="C0")
    axs[0].plot(f_vector_GHz, r_solved.imag,
                label="$reflect_i$ solved",
                color="C1")
    axs[0].set_title('Unknown Reflection Coefficient')
    axs[0].set_xlabel('frequency (GHz)')
    axs[0].set_ylabel('reflection coefficient')
    axs[0].set_xlim(f_range)
    axs[0].set_ylim([-1.0, +1.0])
    axs[0].grid()
    axs[0].legend()

    # Find the actual and corrected transmission coefficient
    l_actual = np.empty((len(f_vector),), dtype=np.complex128)
    l_solved = np.empty((len(f_vector),), dtype=np.complex128)
    for findex, f in enumerate(f_vector):
        l_actual[findex] = c.exp(-calc_actual_line(f) * line_length)
        l_solved[findex] = unknown_line.get_value(f)

    # Line
    axs[1].plot(f_vector_GHz, l_ideal.real,
                label="$line_r$ ideal",
                color="C2", linestyle=":")
    axs[1].plot(f_vector_GHz, l_ideal.imag,
                label="$line_i$ ideal",
                color="C3", linestyle=":")
    axs[1].plot(f_vector_GHz, l_actual.real, 'o',
                label="$line_r$ actual",
                color="C2", markersize=2)
    axs[1].plot(f_vector_GHz, l_actual.imag, 's',
                label="$line_i$ actual",
                color="C3", markersize=2)
    axs[1].plot(f_vector_GHz, l_solved.real,
                label="$line_r$ solved",
                color="C2")
    axs[1].plot(f_vector_GHz, l_solved.imag,
                label="$line_i$ solved",
                color="C3")
    axs[1].set_title('Unknown Transmission Coefficient')
    axs[1].set_xlabel('frequency (GHz)')
    axs[1].set_ylabel('transmission coefficient')
    axs[0].set_xlim(f_range)
    axs[1].set_ylim([-1.0, +1.0])
    axs[1].grid()
    axs[1].legend()

    show()


def plot_correction(f_range, f_vector, measured, actual, corrected):
    """
    Plot the measured, expected and corrected values.
    """
    # scale frequency vector to GHz
    f_vector_GHz = np.asarray(f_vector) / GHz

    # scale f_range to GHz
    f_range = np.asarray(f_range) / GHz

    # S11
    fig, axs = subplots(2, 2)
    axs[0, 0].plot(f_vector_GHz, measured[:, 0, 0].real,
                   label="$s11_r$ measured",
                   color="C0", linestyle=":")
    axs[0, 0].plot(f_vector_GHz, measured[:, 0, 0].imag,
                   label="$s11_i$ measured",
                   color="C1", linestyle=":")
    axs[0, 0].plot(f_vector_GHz, actual[:, 0, 0].real, 'o',
                   label="$s11_r$ actual",
                   color="C0", markersize=2)
    axs[0, 0].plot(f_vector_GHz, actual[:, 0, 0].imag, 's',
                   label="$s11_i$ actual",
                   color="C1", markersize=2)
    axs[0, 0].plot(f_vector_GHz, corrected[:, 0, 0].real,
                   label="$s11_r$ corrected",
                   color="C0")
    axs[0, 0].plot(f_vector_GHz, corrected[:, 0, 0].imag,
                   label="$s11_i$ corrected",
                   color="C1")
    axs[0, 0].set_title('S11')
    axs[0, 0].set_xlabel('frequency (GHz)')
    axs[0, 0].set_ylabel('reflection')
    axs[0, 0].set_xlim(f_range)
    axs[0, 0].set_ylim([-1.0, +1.0])
    axs[0, 0].grid()
    axs[0, 0].legend()

    #
    # S12
    #
    axs[0, 1].plot(f_vector_GHz, measured[:, 0, 1].real,
                   label="$s12_r$ measured",
                   color="C2", linestyle=":")
    axs[0, 1].plot(f_vector_GHz, measured[:, 0, 1].imag,
                   label="$s12_i$ measured",
                   color="C3", linestyle=":")
    axs[0, 1].plot(f_vector_GHz, actual[:, 0, 1].real, 'o',
                   label="$s12_r$ actual",
                   color="C2", markersize=2)
    axs[0, 1].plot(f_vector_GHz, actual[:, 0, 1].imag, 's',
                   label="$s12_i$ actual",
                   color="C3", markersize=2)
    axs[0, 1].plot(f_vector_GHz, corrected[:, 0, 1].real,
                   label="$s12_r$ corrected",
                   color="C2")
    axs[0, 1].plot(f_vector_GHz, corrected[:, 0, 1].imag,
                   label="$s12_i$ corrected",
                   color="C3")
    axs[0, 1].set_title('S12')
    axs[0, 1].set_xlabel('frequency (GHz)')
    axs[0, 1].set_ylabel('transmission')
    axs[0, 1].set_xlim(f_range)
    axs[0, 1].set_ylim([-1.0, +1.0])
    axs[0, 1].grid()
    axs[0, 1].legend()

    #
    # S21
    #
    axs[1, 0].plot(f_vector_GHz, measured[:, 1, 0].real,
                   label="$s21_r$ measured",
                   color="C4", linestyle=":")
    axs[1, 0].plot(f_vector_GHz, measured[:, 1, 0].imag,
                   label="$s21_i$ measured",
                   color="C5", linestyle=":")
    axs[1, 0].plot(f_vector_GHz, actual[:, 1, 0].real, 'o',
                   label="$s21_r$ actual",
                   color="C4", markersize=2)
    axs[1, 0].plot(f_vector_GHz, actual[:, 1, 0].imag, 's',
                   label="$s21_i$ actual",
                   color="C5", markersize=2)
    axs[1, 0].plot(f_vector_GHz, corrected[:, 1, 0].real,
                   label="$s21_r$ corrected",
                   color="C4")
    axs[1, 0].plot(f_vector_GHz, corrected[:, 1, 0].imag,
                   label="$s21_i$ corrected",
                   color="C5")
    axs[1, 0].set_title('S21')
    axs[1, 0].set_xlabel('frequency (GHz)')
    axs[1, 0].set_ylabel('transmission')
    axs[1, 0].set_xlim(f_range)
    axs[1, 0].set_ylim([-1.0, +1.0])
    axs[1, 0].grid()
    axs[1, 0].legend()

    #
    # S22
    #
    axs[1, 1].plot(f_vector_GHz, measured[:, 1, 1].real,
                   label="$s22_r$ measured",
                   color="C6", linestyle=":")
    axs[1, 1].plot(f_vector_GHz, measured[:, 1, 1].imag,
                   label="$s22_i$ measured",
                   color="C7", linestyle=":")
    axs[1, 1].plot(f_vector_GHz, actual[:, 1, 1].real, 'o',
                   label="$s22_r$ actual",
                   color="C6", markersize=2)
    axs[1, 1].plot(f_vector_GHz, actual[:, 1, 1].imag, 's',
                   label="$s22_i$ actual",
                   color="C7", markersize=2)
    axs[1, 1].plot(f_vector_GHz, corrected[:, 1, 1].real,
                   label="$s22_r$ corrected",
                   color="C6")
    axs[1, 1].plot(f_vector_GHz, corrected[:, 1, 1].imag,
                   label="$s22_i$ corrected",
                   color="C7")
    axs[1, 1].set_title('S22')
    axs[1, 1].set_xlabel('frequency (GHz)')
    axs[1, 1].set_ylabel('reflection')
    axs[1, 1].set_xlim(f_range)
    axs[1, 1].set_ylim([-1.0, +1.0])
    axs[1, 1].grid()
    axs[1, 1].legend()

    show()
