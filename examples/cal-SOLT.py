#!/usr/bin/python3
#
# This example shows how to create a calibration file for a flawed
# VNA and how to apply the calibration to measurements of a device
# to return the correct S parameters.  This example closely follows
# vnacal-SOLT-example.c from libvna.
#
from enum import Enum
from math import log10
from numpy import complex128, empty, flipud, hstack, logspace, matmul, ones, pi
from matplotlib import ticker
from matplotlib.pyplot import show, subplots
from libvna.cal import Calset, CalType, Solver
from libvna.data import NPData

#
# Set frequency range.
#
F_MIN = 10.0e+3
F_MAX = 100e+6

#
# Set parameters used for calibration.  Our simple VNA measures only
# S11 and S21 directly, hence C_COLUMNS = 1.
#
C_ROWS = 2
C_COLUMNS = 1
C_FREQUENCIES = 79

#
# Set parameters used for applying the calibration to measurements of
# a device under test (DUT).
#
M_ROWS = 2
M_COLUMNS = 2
M_FREQUENCIES = 100

#
# Set constants for our flawed VNA.
#
Z0 = 50.0
W1 = (2 * pi * 10.0e+6)
W2 = (2 * pi * 1.00e+6)


class Measurement(Enum):
    """
    Which measurement should our simulated VNA return?
    """
    SHORT_CALIBRATION = 0
    OPEN_CALIBRATION = 1
    LOAD_CALIBRATION = 2
    THROUGH_CALIBRATION = 3
    FORWARD_MEASUREMENT = 4
    REVERSE_MEASUREMENT = 5


def measure(measurement: Measurement, f_vector):
    """
    Simulate the requested VNA measurement.
        Our simulated VNA has two flaws: first, there is a
        stray capacitance of 1 / (Z0 * W1) [318pF] between
        port 1 and ground; second, there is an inductance of
        Z0 / W1 [796 nH] in series with port 2.

        The simulated device under test (DUT) is a second
        order LC divider low pass filter with L = Z0 / W2
        [7.96 μH] and C = 1 / (Z0 * W2) [3.18nF].

    Parameters:
        measurement: what to measure
        f_vector: vector of frequency points

    Return:
        (2 x 1 x frequencies) complex array of measurements
          - First index is which detector
          - Second index is which port generates signal (1st only for this VNA)
          - Third is frequency.  Frequency index is innermost rather
            than outermost here to match that many VNAs return single
            vectors of measurements.
    """
    el21 = 0.1j     # leakage error into port 2
    frequencies = len(f_vector)
    if measurement == Measurement.SHORT_CALIBRATION:
        #
        # The shorted calibration standard on port 1 shunts out the
        # stray capacitance, giving a perfect gamma value of -1.  Port 2
        # receives some phase shifted internal leakage signal.
        #
        result = ones(shape=(1, frequencies), dtype=complex128, order="C")
        result = matmul([[-1.0], [el21]], result)
        result = result.reshape((2, 1, frequencies))
        return result

    elif measurement == Measurement.OPEN_CALIBRATION:
        #
        # The open calibration standard exposes the stray capacitance
        # on port 1.  Port 2 continues to pick up internal leakage.
        #
        result = empty(shape=(2, 1, frequencies), dtype=complex128, order="C")
        for findex, f in enumerate(f_vector):
            s = 2j * pi * f
            result[:, 0, findex] = [(1.0 - s / W1) / (1.0 + s / W1), el21]
        return result

    elif measurement == Measurement.LOAD_CALIBRATION:
        #
        # The load calibration is in parallel with the stray capacitance
        # on port 1.  As always, port 2 continues to pick up internal
        # leakage.
        #
        result = empty(shape=(2, 1, frequencies), dtype=complex128, order="C")
        for findex, f in enumerate(f_vector):
            s = 2j * pi * f
            result[:, 0, findex] = [-s / (s + 2 * W1), el21]
        return result

    elif measurement == Measurement.THROUGH_CALIBRATION:
        #
        # In the through configuration, the stray capacitance on
        # port 1 and the stray inductance on port 2 form a resonant
        # circuit with a high-pass reflected signal and low-pass
        # transmitted signal.
        #
        result = empty(shape=(2, 1, frequencies), dtype=complex128, order="C")
        for findex, f in enumerate(f_vector):
            s = 2j * pi * f
            d = s**2 + 2 * W1 * s + 2 * W1**2
            result[:, 0, findex] = [-s**2 / d, 2 * W1**2 / d + el21]
        return result

    elif measurement == Measurement.FORWARD_MEASUREMENT:
        #
        # In the forward configuration, the DUT forms a fourth order
        # resonant circuit with the stray impedances of the VNA.
        #
        result = empty(shape=(2, 1, frequencies), dtype=complex128, order="C")
        for findex, f in enumerate(f_vector):
            s = 2j * pi * f
            d = s**4 + 2 * W1 * s**3 + ((W1 + W2) * s)**2 \
                + 2 * W1 * W2 * (W1 + W2) * s + 2 * (W1 * W2)**2
            result[:, 0, findex] = [
                -(s**4 - (W1**2 - 2 * W1 * W2 - W2**2) * s**2) / d,
                2 * (W1 * W2)**2 / d + el21
            ]
        return result

    elif measurement == Measurement.REVERSE_MEASUREMENT:
        #
        # In the reverse configuration, the stray capacitance on
        # port 1 is in parallel with the DUT capacitor and the stray
        # inductance on port 2 is in series with the DUT inductor,
        # forming a second order resonant circuit.
        #
        result = empty(shape=(2, 1, frequencies), dtype=complex128, order="C")
        for findex, f in enumerate(f_vector):
            s = 2j * pi * f
            d = s**2 + 2 * W1 * W2 / (W1 + W2) * s \
                + 2 * (W1 * W2)**2 / (W1 + W2)**2
            result[:, 0, findex] = [
                -s**2 / d,
                2 * (W1 * W2)**2 / (W1 + W2)**2 / d + el21
            ]
        return result

    else:
        assert(False)


def make_calibration():
    """
    Solve VNA calibration error terms from measurements of short, open,
    match and through standards.
    """

    # Create a new container to hold the calibration.
    calset = Calset()

    # Create the calibration frequency vector.
    f_vector = logspace(log10(F_MIN), log10(F_MAX), C_FREQUENCIES)

    # Create the error term solver with E12 error terms.
    solver = Solver(calset, CalType.E12, C_ROWS, C_COLUMNS, f_vector)

    #
    # Perform each of short, open and load calibration on port 1, assuming
    # a terminator attached to port 2, then perform through calibration
    # between ports 1 and 2.  Normally, we'd have to interact with the
    # user to connect each standard, but in our simulated environment, we
    # can simply call measure to get the measurements for each standard.
    #

    # Short calibration
    m = measure(Measurement.SHORT_CALIBRATION, f_vector)
    solver.add_single_reflect(None, m, -1.0, 1)

    # Open calibration
    m = measure(Measurement.OPEN_CALIBRATION, f_vector)
    solver.add_single_reflect(None, m, 1.0, 1)

    # Load calibration
    m = measure(Measurement.LOAD_CALIBRATION, f_vector)
    solver.add_single_reflect(None, m, 0.0, 1)

    # Through calibration
    m = measure(Measurement.THROUGH_CALIBRATION, f_vector)
    solver.add_through(None, m, 1, 2)

    # Solve for the error terms, add them to the calset and save to a file.
    solver.solve()
    solver.add_to_calset("cal_2x1")
    calset.save("SOLT.vnacal")


def apply_calibration():
    """
    Apply the calibration to measured data.
    """
    #
    # Load the calibration set from a file and select the cal_2x1
    # calibration we created above.
    #
    calset = Calset("SOLT.vnacal")
    calibration = calset.calibrations["cal_2x1"]

    #
    # Make a vector of frequencies for the measurements.  Note that
    # except for the endpoints, the frequency points we're using here
    # don't coincide with those used to make the calibration.  The apply
    # function uses rational function interpolation to interpolate the
    # error terms as needed.
    #
    f_vector = logspace(log10(F_MIN), log10(F_MAX), M_FREQUENCIES)

    #
    # Measure the device under test (DUT) with VNA ports 1 and 2
    # connected to DUT ports 1 and 2, respectively.
    #
    m1 = measure(Measurement.FORWARD_MEASUREMENT, f_vector)

    #
    # Measure the device under test with connections reversed.  Our
    # simulated VNA measures only S11 and S21 directly; by reversing
    # the connections, we get S22 and S12.  Note, that even though we
    # make the raw measurements separately, we need all four S parameters
    # together to apply the correction because reflections caused by match
    # errors on the VNA ports cause S11 and S21 to not be independent
    # S12 and S22.
    #
    # Normally, we would need to interact with the user or control a
    # relay to change the connections.  If the VNA has a relay, though,
    # then we would have created a 2x2 calibration above to calibrate
    # for the relay also.
    #
    m2 = measure(Measurement.REVERSE_MEASUREMENT, f_vector)

    #
    # Concatenate the two columns we measured above to form a
    # (2 x 2 x frequencies) array.  We have to flip m2 up-down
    # because, with ports reversed, VNA port 1 measures S22 and
    # port 2 measures S12.
    #
    measured = hstack((m1, flipud(m2)))

    #
    # Apply the calibration.  We're not using a reference matrix,
    # so we pass "a" as None.
    #
    result = calibration.apply(f_vector, None, measured)

    #
    # Plot the results
    #
    corrected = result.data_array
    make_plots(f_vector, measured, corrected)


def make_plots(f_vector, measured, corrected):
    """
    Plot measured, expected and corrected.
    """

    #
    # Calculate the expected result.
    #
    expected = empty(shape=(M_FREQUENCIES, M_ROWS, M_COLUMNS),
                     dtype=complex128)
    for findex, f in enumerate(f_vector):
        s = 2j * pi * f
        d = s**2 + 2*W2*s + 2*W2**2
        expected[findex, 0, 0] = +s**2 / d
        expected[findex, 0, 1] = +2*W2**2 / d
        expected[findex, 1, 0] = +2*W2**2 / d
        expected[findex, 1, 1] = -s**2 / d

    #
    # S11
    #
    fig, axs = subplots(2, 2)
    axs[0, 0].semilogx(f_vector, measured[0, 0, :].real,
                       label="$s11_r$ measured",
                       color="C0", linestyle=":")
    axs[0, 0].semilogx(f_vector, measured[0, 0, :].imag,
                       label="$s11_i$ measured",
                       color="C1", linestyle=":")
    axs[0, 0].semilogx(f_vector, expected[:, 0, 0].real, 'o',
                       label="$s11_r$ expected",
                       color="C0", markersize=2)
    axs[0, 0].semilogx(f_vector, expected[:, 0, 0].imag, 's',
                       label="$s11_i$ expected",
                       color="C1", markersize=2)
    axs[0, 0].semilogx(f_vector, corrected[:, 0, 0].real,
                       label="$s11_r$ corrected",
                       color="C0")
    axs[0, 0].semilogx(f_vector, corrected[:, 0, 0].imag,
                       label="$s11_i$ corrected",
                       color="C1")
    axs[0, 0].set_title('S11')
    axs[0, 0].set_xlabel('frequency (Hz)')
    axs[0, 0].set_ylabel('reflection')
    axs[0, 0].set_xlim([F_MIN, F_MAX])
    axs[0, 0].set_ylim([-1.0, +1.0])
    axs[0, 0].grid()
    axs[0, 0].legend()

    #
    # S12
    #
    axs[0, 1].semilogx(f_vector, measured[0, 1, :].real,
                       label="$s12_r$ measured",
                       color="C2", linestyle=":")
    axs[0, 1].semilogx(f_vector, measured[0, 1, :].imag,
                       label="$s12_i$ measured",
                       color="C3", linestyle=":")
    axs[0, 1].semilogx(f_vector, expected[:, 0, 1].real, 'o',
                       label="$s12_r$ expected",
                       color="C2", markersize=2)
    axs[0, 1].semilogx(f_vector, expected[:, 0, 1].imag, 's',
                       label="$s12_i$ expected",
                       color="C3", markersize=2)
    axs[0, 1].semilogx(f_vector, corrected[:, 0, 1].real,
                       label="$s12_r$ corrected",
                       color="C2")
    axs[0, 1].semilogx(f_vector, corrected[:, 0, 1].imag,
                       label="$s12_i$ corrected",
                       color="C3")
    axs[0, 1].set_title('S12')
    axs[0, 1].set_xlabel('frequency (Hz)')
    axs[0, 1].set_ylabel('transmission')
    axs[0, 1].set_xlim([F_MIN, F_MAX])
    axs[0, 1].set_ylim([-1.0, +1.0])
    axs[0, 1].grid()
    axs[0, 1].legend()

    #
    # S21
    #
    axs[1, 0].semilogx(f_vector, measured[1, 0, :].real,
                       label="$s21_r$ measured",
                       color="C4", linestyle=":")
    axs[1, 0].semilogx(f_vector, measured[1, 0, :].imag,
                       label="$s21_i$ measured",
                       color="C5", linestyle=":")
    axs[1, 0].semilogx(f_vector, expected[:, 1, 0].real, 'o',
                       label="$s21_r$ expected",
                       color="C4", markersize=2)
    axs[1, 0].semilogx(f_vector, expected[:, 1, 0].imag, 's',
                       label="$s21_i$ expected",
                       color="C5", markersize=2)
    axs[1, 0].semilogx(f_vector, corrected[:, 1, 0].real,
                       label="$s21_r$ corrected",
                       color="C4")
    axs[1, 0].semilogx(f_vector, corrected[:, 1, 0].imag,
                       label="$s21_i$ corrected",
                       color="C5")
    axs[1, 0].set_title('S21')
    axs[1, 0].set_xlabel('frequency (Hz)')
    axs[1, 0].set_ylabel('transmission')
    axs[1, 0].set_xlim([F_MIN, F_MAX])
    axs[1, 0].set_ylim([-1.0, +1.0])
    axs[1, 0].grid()
    axs[1, 0].legend()

    #
    # S22
    #
    axs[1, 1].semilogx(f_vector, measured[1, 1, :].real,
                       label="$s22_r$ measured",
                       color="C6", linestyle=":")
    axs[1, 1].semilogx(f_vector, measured[1, 1, :].imag,
                       label="$s22_i$ measured",
                       color="C7", linestyle=":")
    axs[1, 1].semilogx(f_vector, expected[:, 1, 1].real, 'o',
                       label="$s22_r$ expected",
                       color="C6", markersize=2)
    axs[1, 1].semilogx(f_vector, expected[:, 1, 1].imag, 's',
                       label="$s22_i$ expected",
                       color="C7", markersize=2)
    axs[1, 1].semilogx(f_vector, corrected[:, 1, 1].real,
                       label="$s22_r$ corrected",
                       color="C6")
    axs[1, 1].semilogx(f_vector, corrected[:, 1, 1].imag,
                       label="$s22_i$ corrected",
                       color="C7")
    axs[1, 1].set_title('S22')
    axs[1, 1].set_xlabel('frequency (Hz)')
    axs[1, 1].set_ylabel('reflection')
    axs[1, 1].set_xlim([F_MIN, F_MAX])
    axs[1, 1].set_ylim([-1.0, +1.0])
    axs[1, 1].grid()
    axs[1, 1].legend()

    show()


make_calibration()
apply_calibration()
