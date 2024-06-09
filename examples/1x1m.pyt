%# Imports for the error term and test data generator.
%[
from libvna.cal import Calset, CalType, VectorParameter
from libvna.data import NPData, PType
import math
import numpy as np
import random_error_terms as et


%# Set up random error parameters for the VNA.
# Set up random generator.
rng = np.random.default_rng(seed=2)

# Set frequency range and number of points
fmin = 1.0e+6
fmax = 1.0e+9
points = 10
f_vector = np.logspace(np.log10(fmin), np.log10(fmax), num=points)
eterms = et.RandomErrorTerms(rng, CalType.E12, 1, 1, fmin, fmax)

def parallel(*args):
    '''
    Find parallel combination of impedances.
    '''
    sum = 0.0j
    for arg in args:
        sum += arg**(-1)
    return sum**(-1)

# Set up a different frequency range for the standards.
std_points = 50
f_std = np.logspace(np.log10(100e+3), np.log10(10e+9), num=std_points)
calset = Calset()

# Short standard is (5 ohms + 250nH) || 1pF
d_short = NPData(PType.Z, frequencies=std_points, rows=1, columns=1)
d_short.frequency_vector = f_std
for i, f in enumerate(f_std):
    jω = 2j * math.pi * f
    z = parallel(5.0 + jω * 250e-9, 1.0 / (jω * 1.0e-12))
    d_short.data_array[i, 0, 0] = z
s_short = [[VectorParameter(calset, f_std,
                            d_short.convert(PType.S).data_array[:, 0, 0])]]
d_short.save('short.s1p')

# Open standard is 100pF in series with 2.5nH
d_open = NPData(PType.Z, frequencies=std_points, rows=1, columns=1)
d_open.frequency_vector = f_std
for i, f in enumerate(f_std):
    jω = 2j * math.pi * f
    z = 1 / (jω * 100e-12) + jω * 2.5e-9
    d_open.data_array[i, 0, 0] = z
s_open = [[VectorParameter(calset, f_std,
                           d_open.convert(PType.S).data_array[:, 0, 0])]]
d_open.save('open.s1p')

# Load standard is 75 ohms in parallel with 5pF
d_load = NPData(PType.Z, frequencies=std_points, rows=1, columns=1)
d_load.frequency_vector = f_std
for i, f in enumerate(f_std):
    jω = 2j * math.pi * f
    z = parallel(75.0, 1.0 / (jω  * 5.0e-12))
    d_open.data_array[i, 0, 0] = z
s_load = [[VectorParameter(calset, f_std,
                           d_load.convert(PType.S).data_array[:, 0, 0])]]
d_load.save('load.s1p')
%]
%O 1x1m-calibrate.py
%############################ begin calibration ###############################
%# Imports for calibration
from libvna.cal import Calset, CalType, Solver, VectorParameter
from libvna.data import NPData, PType
import numpy as np


# Set the calibration frequency points.
fmin = %{fmin:7.1e%}
fmax = %{fmax:7.1e%}
f_vector = np.logspace(np.log10(fmin), np.log10(fmax), num=%{points%})

# Set up libvna.cal's error term solver.
calset = Calset()
solver = Solver(calset, CalType.E12, rows=1, columns=1,
                frequency_vector=f_vector)

# Load the parameters of the short standard, convert to S-parameters
# (if not already in S), and form into a VectorParameter.
short = NPData(filename='short.s1p', ptype=PType.S)
s11 = VectorParameter(calset,
                      short.frequency_vector,
                      short.data_array[:, 0, 0])

# Add measurement of the short standard.
%[
m = eterms.evaluate(calset, f_vector, s_short)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_single_reflect(m, s11)

# Load the parameters of the open standard, convert to S-parameters,
# and form into a VectorParameter.
open = NPData(filename='open.s1p', ptype=PType.S)
s11 = VectorParameter(calset,
                      open.frequency_vector,
                      open.data_array[:, 0, 0])

# Add measurement of the open standard.
%[
m = eterms.evaluate(calset, f_vector, s_open)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_single_reflect(m, s11)

# Load the parameters of the load standard, convert to S-parameters,
# and form into a VectorParameter.
load = NPData(filename='load.s1p', ptype=PType.S)
s11 = VectorParameter(calset,
                      load.frequency_vector,
                      load.data_array[:, 0, 0])

# Add measurement of the load standard.
%[
m = eterms.evaluate(calset, f_vector, s_load)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_single_reflect(m, s11)

# Solve, add to Calset and save.
solver.solve()
solver.add_to_calset('1-port')
calset.save('1x1m.vnacal')
%############################# end calibration ################################
%#
%#
%O 1x1m-apply.py
%############################# begin application ##############################
%# Imports for application
from libvna.cal import Calset


# Load the calibration from file.
calset = Calset('1x1m.vnacal')
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
npd.save('1x1m-expected.s1p')
%]

# Measured response:
%[
s = [[(f_vector, expected[:, 0, 0])]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, name="measured", indent=_indent)
npd.data_array = expected
npd.save('1x1m-measured.s1p')
%]

# Apply calibration and save in Touchstone format.
result = calibration.apply(f_vector, measured)
result.save('1x1m-corrected.s1p')
%############################## end application ###############################
