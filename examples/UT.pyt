%# Imports for the error term and test data generator.
%[
from libvna.cal import CalType, Calset, VectorParameter
from libvna.data import NPData, PType
from libvna.conv import ztos
import math
import numpy as np
import random_error_terms as et


%# Set up random error parameters for the VNA.
# Set up random generator.
rng = np.random.default_rng(seed=7)

# Set frequency range and number of points
fmin = 1.0e+6
fmax = 1.0e+9
points = 10
f_vector = np.logspace(np.log10(fmin), np.log10(fmax), num=points)

# Generate random error terms
eterms = et.RandomErrorTerms(rng, CalType.TE10, 2, 2, fmin, fmax)
calset = Calset()

# Find the actual through standard modelling it as a length of
# transmission line with both metal and dielectric losses.
length = 0.02			# length in meters
vf = 2.0 / 3.0			# velocity factor
lm = 2.50e-4			# metal loss Np/m/Hz**(1/2)
ld = 5.86e-9			# dielectric loss (Np/m/Hz)
line_loss = 0.02		# dB per mm per sqrt(f_GHz)
mm_per_m = 1000.0		# mm per meter
np_per_db = 0.11512925		# neper to dB
c = 2.9979246e+8		# speed of light (m/s)
gl = [(lm * math.sqrt(f) + ld * f + 2.0j * math.pi * f / (vf * c)) * length
      for f in f_vector]
γ_through_actual = np.exp(-np.asarray(gl))
T_actual = VectorParameter(calset, f_vector, γ_through_actual)
%]
%O UT-calibrate.py
%############################ begin calibration ###############################
%# Imports for calibration
from libvna.cal import Calset, CalType, Solver, UnknownParameter
from libvna.data import NPData, PType
import math
import numpy as np


# Set the calibration frequency points.
fmin = %{fmin:7.1e%}
fmax = %{fmax:7.1e%}
f_vector = np.logspace(np.log10(fmin), np.log10(fmax), num=%{points%})

# Set up libvna.cal's error term solver.
calset = Calset()
solver = Solver(calset, CalType.TE10, rows=2, columns=2,
                frequency_vector=f_vector)

# Create an unknown parameter for our through standard.  As long as the
# phase shift in the through remains safely below 90 degrees, we can
# give 1 as the estimate.  For a larger phase shift, we would need to
# give a more accurate estimate for each frequency.
T = UnknownParameter(calset, 1)

# Add measurement of the short-open standard.
%[
s = [[-1, 0],
     [0, 1]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_double_reflect(m, s11=-1, s22=1)

# Add measurement of the open-match standard.
%[
s = [[1, 0],
     [0, 0]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_double_reflect(m, s11=1, s22=0)

# Add measurement of the match-short standard.
%[
s = [[0, 0],
     [0, -1]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_double_reflect(m, s11=0, s22=-1)


# Add measurement of the unknown through standard.  Note that we
# have to add it as a "line".
%[
s = [[0, T_actual],
     [T_actual, 0]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_line(m, [[0, T], [T, 0]])

# Solve, add to Calset and save.
solver.solve()
solver.add_to_calset('mycal')
calset.save('UT.vnacal')

# Save the solved T
t_array = np.asarray(T.get_value(f_vector)).reshape((len(f_vector), 1, 1))
t_data = NPData(PType.S, len(f_vector), 1, 1)
t_data.frequency_vector = f_vector
t_data.data_array = t_array
t_data.format = "Sma"
t_data.save('UT-T.s1p')
%############################# end calibration ################################
%#
%#
%O UT-apply.py
%############################# begin application ##############################
%# Imports for calibration
from libvna.cal import Calset


# Load the calibration from file.
calset = Calset('UT.vnacal')
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
npd = NPData(PType.S, frequencies=len(f_vector), rows=2, columns=2)
npd.frequency_vector = f_vector
npd.data_array = expected
npd.save('UT-expected.s2p')
%]
# Measured response
%[
s = [[VectorParameter(calset, f_vector, expected[:, 0, 0]),
      VectorParameter(calset, f_vector, expected[:, 0, 1])],
     [VectorParameter(calset, f_vector, expected[:, 1, 0]),
      VectorParameter(calset, f_vector, expected[:, 1, 1])]]
m = eterms.evaluate(calset, f_vector, s)
npd.data_array = m
npd.save('UT-measured.s2p')
et.print_matrix(m, file=_file, indent=_indent)
%]

# Apply calibration and save in Touchstone format.
result = calibration.apply(f_vector, m)
result.save('UT-corrected.s2p')
%############################## end application ###############################
