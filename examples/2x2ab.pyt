%# Imports for the error term and test data generator.
%[
from libvna.cal import CalType, Calset
from libvna.data import NPData, PType
from libvna.conv import ztos
import math
import numpy as np
import random_error_terms as et


%# Set up random error parameters for the VNA.
# Set up random generator.
rng = np.random.default_rng(seed=4)

# Set frequency range and number of points
fmin = 1.0e+6
fmax = 1.0e+9
points = 10
f_vector = np.logspace(np.log10(fmin), np.log10(fmax), num=points)
eterms = et.RandomErrorTerms(rng, CalType.TE10, 2, 2, fmin, fmax)
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
through_std = calset.through_standard(
    offset_z0=50.280012,
    offset_delay=75.191505e-12,
    offset_loss=6.731522e+9
)
%]
%O 2x2ab-calibrate.py
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
solver = Solver(calset, CalType.TE10, rows=2, columns=2,
                frequency_vector=f_vector)

# Define the standards for our calibration kit.
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
through_std = calset.through_standard(
    offset_z0=50.280012,
    offset_delay=75.191505e-12,
    offset_loss=6.731522e+9
)

# Add measurement of the short-open standard..
%[
s = [[short_std, 0],
     [ 0, open_std]]
a, b = eterms.evaluate(calset, f_vector, s, ab=True)
et.print_matrix(a, file=_file, name='a', indent=_indent)
et.print_matrix(b, file=_file, name='b', indent=_indent)
%]
solver.add_double_reflect(a=a, b=b, s11=short_std, s22=open_std)

# Add measurement of the short-load stnadard.
%[
s = [[short_std, 0],
     [0, load_std]]
a, b = eterms.evaluate(calset, f_vector, s, ab=True)
et.print_matrix(a, file=_file, name='a', indent=_indent)
et.print_matrix(b, file=_file, name='b', indent=_indent)
%]
solver.add_double_reflect(a=a, b=b, s11=short_std, s22=load_std)

# Add measurement of the through standard.  Notice that we have to add
# the calkit standard as a 'line'.
%[
s = through_std
a, b = eterms.evaluate(calset, f_vector, s, ab=True)
et.print_matrix(a, file=_file, name='a', indent=_indent)
et.print_matrix(b, file=_file, name='b', indent=_indent)
%]
solver.add_line(a=a, b=b, s=through_std)

# Solve, add to Calset and save.
solver.solve()
solver.add_to_calset('cal2x2')
calset.save('2x2ab.vnacal')
%############################# end calibration ################################
%#
%#
%O 2x2ab-apply.py
%############################# begin application ##############################
%# Imports for calibration
from libvna.cal import Calset


# Load the calibration from file.
calset = Calset('2x2ab.vnacal')
calibration = calset.calibrations[0]
f_vector = calibration.frequency_vector

%[
# Calculate and save the expected response of an LC low pass filter
# with 250nH and 100pF.
expected = np.empty((points, 2, 2), dtype=complex)
for i, f in enumerate(f_vector):
    jω = 2j * math.pi * f
    z1 = jω * 250e-9
    z2 = 1.0 / (jω * 100e-12)
    z = [[z1+z2, z2],
         [z2,    z2]]
    expected[i, :, :] = ztos(z)
npd = NPData(PType.S, rows=2, columns=2, frequencies=len(f_vector))
npd.frequency_vector = f_vector
npd.data_array = expected
npd.save('2x2ab-expected.s2p')
%]
# Measured response
%[
s = [[calset.vector_parameter(f_vector, expected[:, 0, 0]),
      calset.vector_parameter(f_vector, expected[:, 0, 1])],
     [calset.vector_parameter(f_vector, expected[:, 1, 0]),
      calset.vector_parameter(f_vector, expected[:, 1, 1])]]
a, b = eterms.evaluate(calset, f_vector, s, ab=True)
et.print_matrix(a, file=_file, name='a', indent=_indent)
et.print_matrix(b, file=_file, name='b', indent=_indent)
npd.data_array = b @ np.linalg.inv(a)
npd.save('2x2ab-measured.s2p')
%]

# Apply calibration and save in Touchstone format.
result = calibration.apply(f_vector, a=a, b=b)
result.save('2x2ab-corrected.s2p')
%############################## end application ###############################
