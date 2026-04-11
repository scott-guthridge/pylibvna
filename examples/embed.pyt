%# Imports for the error term and test data generator.
%[
from libvna.cal import CalType, Calset
from libvna.data import NPData, PType
import math
import numpy as np
import random_error_terms as et


%# Set up random error parameters for the VNA.
# Set up random generator.
rng = np.random.default_rng(seed=1)

# Set frequency range and number of points
fmin = 100.0e+6
fmax = 8.5e+9
points = 10
f_vector = np.logspace(np.log10(fmin), np.log10(fmax), num=points)
eterms = et.RandomErrorTerms(rng, CalType.E12, 1, 1, fmin, fmax)
calset = Calset()

%# Test harness needs the standards, too.
short_standard = calset.short_standard(
    offset_z0=50.0,
    offset_delay=33.356e-15,
    offset_loss=2.36e+9,
    fmax=8.5e+9,
    L=[-44.000e-12, 3.7000e-21, -250.00e-33, 5.0000e-42]
)
open_standard = calset.open_standard(
    offset_z0=50.0,
    offset_delay=32.032e-15,
    offset_loss=2.2e+9,
    fmax=8.5e+9,
    C=[-17.500e-15, -2.0000e-24, 140.00e-36, -2.7000e-45]
)
load_standard = calset.load_standard(
    offset_z0=50.0,
    offset_delay=0.0,
    offset_loss=0.0,
    fmax=8.5e+9,
    Zl=50.0
)
adapter = calset.through_standard(
    offset_z0=50.0,
    offset_delay=59.0e-12,
    offset_loss=2.51e+9,
    fmax=8.5e+9
)
embedded_short = short_standard.embed(fixture=adapter)
embedded_open = open_standard.embed(fixture=adapter)
embedded_load = load_standard.embed(fixture=adapter)

%]
%O embed-calibrate.py
%############################ begin calibration ###############################
%# Imports for calibration
from libvna.cal import Calset, CalType, Solver
import numpy as np


# Set the calibration frequency points.
fmin = %{fmin:7.1e%}
fmax = %{fmax:7.1e%}
f_vector = np.logspace(np.log10(fmin), np.log10(fmax), num=%{points%})

# Create the calibration container
calset = Calset()

# Make calibration standards.  These are 3.5mm female.
short_standard = calset.short_standard(
    offset_z0=50.0,
    offset_delay=33.356e-15,
    offset_loss=2.36e+9,
    fmax=8.5e+9,
    L=[-44.000e-12, 3.7000e-21, -250.00e-33, 5.0000e-42]
)
open_standard = calset.open_standard(
    offset_z0=50.0,
    offset_delay=32.032e-15,
    offset_loss=2.2e+9,
    fmax=8.5e+9,
    C=[-17.500e-15, -2.0000e-24, 140.00e-36, -2.7000e-45]
)
load_standard = calset.load_standard(
    offset_z0=50.0,
    offset_delay=0.0,
    offset_loss=0.0,
    fmax=8.5e+9,
    Zl=50.0
)

# The VNA expects male standards, so we insert this male-male adapter,
# keeping the VNA reference plane on the VNA side of the adapter.
adapter = calset.through_standard(
    offset_z0=50.0,
    offset_delay=59.0e-12,
    offset_loss=2.51e+9,
    fmax=8.5e+9
)

# Embed the standards.
embedded_short = short_standard.embed(fixture=adapter)
embedded_open = open_standard.embed(fixture=adapter)
embedded_load = load_standard.embed(fixture=adapter)

# Create the error term solver.
solver = Solver(calset, CalType.E12, rows=1, columns=1,
                frequency_vector=f_vector)

# Add measurement of the short_standard standard.
%[
s = [[embedded_short]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_single_reflect(m, s11=embedded_short)

# Add measurement of the open_standard standard.
%[
s = [[embedded_open]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_single_reflect(m, s11=embedded_open)

# Add measurement of the load_standard standard.
%[
s = [[embedded_load]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_single_reflect(m, s11=embedded_load)

# Solve, add to Calset and save.
solver.solve()
solver.add_to_calset('1-port')
calset.save('embed.vnacal')
%############################# end calibration ################################
%#
%#
%O embed-apply.py
%############################# begin application ##############################
%# Imports for application
from libvna.cal import Calset


# Load the calibration from file.
calset = Calset('embed.vnacal')
calibration = calset.calibrations['1-port']
f_vector = calibration.frequency_vector

# 3.5mm male-to-male adapter to de-embed
adapter = calset.through_standard(
    offset_z0=50.0,
    offset_delay=59.0e-12,
    offset_loss=2.51e+9,
    fmax=8.5e+9
)

%[
# Calculate and save the expected response.
jω = 2j * math.pi * f_vector
expected = ((-2.2819e+18 + 1.3511e+09 * jω - jω**2) /
             (3.4228e+18 + 2.0267e+09 * jω + jω**2))
expected = expected.reshape((len(f_vector), 1, 1))
npd = NPData(PType.S, frequencies=len(f_vector), rows=1, columns=1)
npd.frequency_vector = f_vector
npd.data_array = expected
npd.save('embed-expected.s1p')
%]
# Measured response:
%[
p = calset.vector_parameter(f_vector, expected[:, 0, 0])
s = [[p.embed(fixture=adapter)]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, name="measured", indent=_indent)
npd.data_array = m
npd.save('embed-measured.s1p')
%]

# Apply calibration.
result = calibration.apply(f_vector, measured)

# De-embed the adapter and save in Touchstone format.
result = calset.deembed_npdata(result, fixture=adapter)
result.save('embed-corrected.s1p')
%############################## end application ###############################
