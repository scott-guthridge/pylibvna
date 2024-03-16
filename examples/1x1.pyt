%# Imports for the error term and test data generator.
%[
from libvna.cal import CalType, Calset
from libvna.data import NPData, PType
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
eterms = et.RandomErrorTerms(CalType.E12, 1, 1, fmin, fmax)
calset = Calset()
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

# Set up libvna.cal's error term solver.
calset = Calset()
solver = Solver(calset, CalType.E12, rows=1, columns=1,
                frequency_vector=f_vector)

# Add measurement of the short standard.
%[
s = [[-1]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_single_reflect(m, s11=-1)

# Add measurement of the open standard.
%[
s = [[1]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_single_reflect(m, s11=1)

# Add measurement of the load standard.
%[
s = [[0]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_single_reflect(m, s11=0)

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
npd = NPData(PType.S, rows=1, columns=1, frequencies=len(f_vector))
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
