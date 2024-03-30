%# Imports for the error term and test data generator.
%[
from libvna.cal import CalType, Calset
from libvna.conv import ztos
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
eterms1 = et.RandomErrorTerms(CalType.E12, 1, 1, fmin, fmax)
eterms2 = et.RandomErrorTerms(CalType.E12, 1, 1, fmin, fmax)
calset = Calset()
%]
%O 2PR-calibrate.py
%############################ begin calibration ###############################
%# Imports for calibration
from libvna.cal import Calset, CalType, Solver
import numpy as np


# Set the calibration frequency points.
fmin = %{fmin:7.1e%}
fmax = %{fmax:7.1e%}
f_vector = np.logspace(np.log10(fmin), np.log10(fmax), num=%{points%})

# Create two error solvers.
calset = Calset()
solver1 = Solver(calset, CalType.E12, rows=1, columns=1,
                 frequency_vector=f_vector)
solver2 = Solver(calset, CalType.E12, rows=1, columns=1,
                 frequency_vector=f_vector)

# Add measurement of load standard on port 1 and short standard on port 2.
%[
s = [[0]]
m1 = eterms1.evaluate(calset, f_vector, s)
s = [[-1]]
m2 = eterms2.evaluate(calset, f_vector, s)
m = np.concatenate((m1, m2), axis=2)
et.print_matrix(m, file=_file, indent=_indent, asarray=True)
%]
solver1.add_single_reflect(m[:, :, 0].reshape((len(f_vector), 1, 1)), s11=0)
solver2.add_single_reflect(m[:, :, 1].reshape((len(f_vector), 1, 1)), s11=-1)

# Add measurement of short standard on port 1 and open standard on port 2.
%[
s = [[-1]]
m1 = eterms1.evaluate(calset, f_vector, s)
s = [[1]]
m2 = eterms2.evaluate(calset, f_vector, s)
m = np.concatenate((m1, m2), axis=2)
et.print_matrix(m, file=_file, indent=_indent, asarray=True)
%]
solver1.add_single_reflect(m[:, :, 0].reshape((len(f_vector), 1, 1)), s11=-1)
solver2.add_single_reflect(m[:, :, 1].reshape((len(f_vector), 1, 1)), s11=1)

# Add measurement of open standard on port 1 and load standard on port 2.
%[
s = [[1]]
m1 = eterms1.evaluate(calset, f_vector, s)
s = [[0]]
m2 = eterms2.evaluate(calset, f_vector, s)
m = np.concatenate((m1, m2), axis=2)
et.print_matrix(m, file=_file, indent=_indent, asarray=True)
%]
solver1.add_single_reflect(m[:, :, 0].reshape((len(f_vector), 1, 1)), s11=1)
solver2.add_single_reflect(m[:, :, 1].reshape((len(f_vector), 1, 1)), s11=0)

# Solve both calibrations, add to Calset and save.
solver1.solve()
solver2.solve()
solver1.add_to_calset('port 1')
solver2.add_to_calset('port 2')
calset.save('2PR.vnacal')
%############################# end calibration ################################
%#
%#
%O 2PR-apply.py
%############################# begin application ##############################
%# Imports for application
from libvna.cal import Calset
import numpy as np


# Load the calibration from file.
calset = Calset('2PR.vnacal')
calibration1 = calset.calibrations['port 1']
calibration2 = calset.calibrations['port 2']
f_vector = calibration1.frequency_vector[...]

%[
# Calculate and save the expected responses.  On port 1, we are
# measuring a 42pF capacitor; on port 2, we're measureing a 100nF
# inductor.
jω = 2j * math.pi * f_vector
z1 = ((jω * 42.e-12)**-1).reshape((points, 1, 1))
z2 = (jω * 100.e-9).reshape((points, 1, 1))
expected1 = ztos(z1, [[50]]*points)
expected2 = ztos(z2, [[50]]*points)
npd = NPData(PType.S, frequencies=len(f_vector), rows=1, columns=1)
npd.frequency_vector = f_vector
npd.data_array = expected1
npd.save('2PR1-expected.s1p')
npd.data_array = expected2
npd.save('2PR2-expected.s1p')
%]
# Measured response:
%[
s = [[(f_vector, expected1[:, 0, 0])]]
m1 = eterms1.evaluate(calset, f_vector, s)
s = [[(f_vector, expected2[:, 0, 0])]]
m2 = eterms2.evaluate(calset, f_vector, s)
m = np.concatenate((m1, m2), axis=2)
et.print_matrix(m, file=_file, name="measured", indent=_indent, asarray=True)
npd.data_array = m1
npd.save('2PR1-measured.s1p')
npd.data_array = m2
npd.save('2PR2-measured.s1p')
%]

# Apply calibration to each port and save in Touchstone format.
m1 = measured[:, 0, 0].reshape((len(f_vector), 1, 1))
m2 = measured[:, 0, 1].reshape((len(f_vector), 1, 1))
result1 = calibration1.apply(f_vector, m1)
result1.save('2PR1-corrected.s1p')
result2 = calibration2.apply(f_vector, m2)
result2.save('2PR2-corrected.s1p')
%############################## end application ###############################
