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
fmin = 1.0e+6
fmax = 1.0e+9
points = 10
f_vector = np.logspace(np.log10(fmin), np.log10(fmax), num=points)
eterms = et.RandomErrorTerms(rng, CalType.E12, 1, 1, fmin, fmax)
calset = Calset()

%# Test harness needs the standards, too.
short_std = calset.short_standard(
    offset_z0=51.259595,
    offset_delay=33.340790e-12,
    offset_loss=5.460953e+9,
    fmax=8.5e+9,
    L=[-119.006943e-12, -1.310397249e-21, 1.511773982e-31, -91.480400e-42]
)
open_std = calset.open_standard(
    offset_z0=51.635682,
    offset_delay=37.636850e-12,
    offset_loss=5.611771e+9,
    fmax=8.5e+9,
    C=[-93.936119e-15, 151.860439e-27, -786.852853e-36, 46.121820e-45]
)
load_std = calset.load_standard(
    offset_z0=50.0,
    offset_delay=0.0,
    offset_loss=0.0,
    fmax=8.5e+9,
    Zl=50.0
)
%]
%O 1x1-calibrate.py
%############################ begin calibration ###############################
%# Imports for calibration
from libvna.cal import Calset, CalType, Solver
import numpy as np


# Set the calibration frequency points.
fmin = %{fmin:7.1e%}
fmax = %{fmax:7.1e%}
f_vector = np.logspace(np.log10(fmin), np.log10(fmax), num=%{points%})

# Create the calibration container and error term solver.
calset = Calset()
solver = Solver(calset, CalType.E12, rows=1, columns=1,
                frequency_vector=f_vector)

# Define the calibration standards for out calibration kit.
short_std = calset.short_standard(
    offset_z0=51.259595,
    offset_delay=33.340790e-12,
    offset_loss=5.460953e+9,
    fmax=8.5e+9,
    L=[-119.006943e-12, -1.310397249e-21, 1.511773982e-31, -91.480400e-42]
)
open_std = calset.open_standard(
    offset_z0=51.635682,
    offset_delay=37.636850e-12,
    offset_loss=5.611771e+9,
    fmax=8.5e+9,
    C=[-93.936119e-15, 151.860439e-27, -786.852853e-36, 46.121820e-45]
)
load_std = calset.load_standard(
    offset_z0=50.0,
    offset_delay=0.0,
    offset_loss=0.0,
    fmax=8.5e+9,
    Zl=50.0
)

# Add measurement of the short standard.
%[
s = [[short_std]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_single_reflect(m, s11=short_std)

# Add measurement of the open standard.
%[
s = [[open_std]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_single_reflect(m, s11=open_std)

# Add measurement of the load standard.
%[
s = [[load_std]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_single_reflect(m, s11=load_std)

# Solve, add to Calset and save.
solver.solve()
solver.add_to_calset('1-port')
calset.save('1x1.vnacal')
%############################# end calibration ################################
%#
%#
%O 1x1-apply.py
%############################# begin application ##############################
%# Imports for application
from libvna.cal import Calset


# Load the calibration from file.
calset = Calset('1x1.vnacal')
calibration = calset.calibrations['1-port']
f_vector = calibration.frequency_vector

%[
# Calculate and save the expected response.
jω = 2j * math.pi * f_vector
expected = ((-3.1583e+16 + 1.5895e+08 * jω - jω**2) /
             (4.7374e+16 + 2.3843e+08 * jω + jω**2))
expected = expected.reshape((len(f_vector), 1, 1))
npd = NPData(PType.S, frequencies=len(f_vector), rows=1, columns=1)
npd.frequency_vector = f_vector
npd.data_array = expected
npd.save('1x1-expected.s1p')
%]
# Measured response:
%[
s = [[(f_vector, expected[:, 0, 0])]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, name="measured", indent=_indent)
npd.data_array = m
npd.save('1x1-measured.s1p')
%]

# Apply calibration and save in Touchstone format.
result = calibration.apply(f_vector, measured)
result.save('1x1-corrected.s1p')
%############################## end application ###############################
