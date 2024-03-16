%# Imports for the error term and test data generator.
%[
from libvna.cal import CalType, Calset, VectorParameter
from libvna.data import NPData, PType
from libvna.conv import ztos
import math
import numpy as np
import random
import random_error_terms as et


%# Set up random error paramters for the VNA.
# Produce consistent results.
random.seed(1)

# Set frequency range and number of points
fmin = 1.0e+6
fmax = 1.0e+9
points = 10
f_vector = np.logspace(np.log10(fmin), np.log10(fmax), num=points)
eterms = et.RandomErrorTerms(CalType.E12, 2, 1, fmin, fmax)
calset = Calset()
%]
%O 2x1-calibrate.py
%############################ begin calibration ###############################
%# Imports for calibration
from libvna.cal import Calset, CalType, Solver
import numpy as np


# Set the calibration frequency points.
fmin = %{fmin:7.1e%}
fmax = %{fmax:7.1e%}
f_vector = np.logspace(np.log10(fmin), np.log10(fmax), num=%{points%})

# Set up libvna.cal's error term solver.
calset = Calset()
solver = Solver(calset, CalType.E12, rows=2, columns=1,
                frequency_vector=f_vector)

# Add measurement of the short standard.
%[
s = [[-1, 0],
     [ 0, 0]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_single_reflect(m, s11=-1)

# Add measurement of the open standard.
%[
s = [[1, 0],
     [0, 0]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_single_reflect(m, s11=1)

# Add measurement of the load standard.
%[
s = [[0, 0],
     [0, 0]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_single_reflect(m, s11=0)

# Add measurement of the through standard.
%[
s = [[0, 1],
     [1, 0]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_through(m)

# Solve, add to Calset and save.
solver.solve()
solver.add_to_calset('cal2x1')
calset.save('2x1.vnacal')
%############################# end calibration ################################
%#
%#
%O 2x1-apply.py
%############################# begin application ##############################
%# Imports for calibration
from libvna.cal import Calset
import numpy as np


# Load the calibration from file.
calset = Calset('2x1.vnacal')
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
npd.save('2x1-expected.s2p')
%]
# Measured response with VNA port1 = DUT port 1, VNA port2 = DUT port 2:
%[
s = [[VectorParameter(calset, f_vector, expected[:, 0, 0]),
      VectorParameter(calset, f_vector, expected[:, 0, 1])],
     [VectorParameter(calset, f_vector, expected[:, 1, 0]),
      VectorParameter(calset, f_vector, expected[:, 1, 1])]]
m1 = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m1, file=_file, name="m1", indent=_indent, asarray=True)
%]

# Measured response with VNA port1 = DUT port 2, VNA port2 = DUT port 1:
%[
s = np.flipud(np.fliplr(s))
m2 = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m2, file=_file, name="m2", indent=_indent, asarray=True)

measured = np.concatenate((m1, np.flip(m2, axis=1)), axis=2)
npd.data_array = measured
npd.save('2x1-measured.s2p')
%]

# Combine the measurements into a vector of 2x2 matrices,
# flipping the rows of m2 in the reversed measurement.
measured = np.concatenate((m1, np.flip(m2, axis=1)), axis=2)

# Apply calibration and save in Touchstone format.
result = calibration.apply(f_vector, measured)
result.save('2x1-corrected.s2p')
%############################## end application ###############################
