%# Imports for the error term and test data generator.
%[
from libvna.cal import CalType, Calset, VectorParameter
from libvna.data import NPData, PType
from libvna.conv import ztos
import math
import numpy as np
import random_error_terms as et


# Set up random generator.
rng = np.random.default_rng(seed=8)

#
# Set frequency range and number of points.
#
fmin = 1.0e+6
fmax = 1.0e+8
points = 10
f_vector = np.linspace(fmin, fmax, num=points)

#
# Make measurement error and connection nonrepeatability error vectors.
# These should really be measured values.
#
nf_vector = np.linspace(1.0e-4, 2.0e-4, num=len(f_vector))
tr_vector = np.linspace(1.0e-3, 2.0e-3, num=len(f_vector))
σ_vector = np.linspace(0.01, 0.02, points)

#
# Generate random error terms.
#
eterms = et.RandomErrorTerms(rng, CalType.T16, 2, 2, fmin, fmax,
 	                     nf_vector=nf_vector, tr_vector=tr_vector)
calset = Calset()

#
# Generate the actual parameters of the standards, including random
# connection non-repeatability errors.
#
actual_M1 = [ 0 + et.random_complex(rng, σ) for σ in σ_vector]
actual_S1 = [-1 + et.random_complex(rng, σ) for σ in σ_vector]
actual_O1 = [+1 + et.random_complex(rng, σ) for σ in σ_vector]
actual_M2 = [ 0 + et.random_complex(rng, σ) for σ in σ_vector]
actual_S2 = [-1 + et.random_complex(rng, σ) for σ in σ_vector]
actual_O2 = [+1 + et.random_complex(rng, σ) for σ in σ_vector]
actual_M3 = [ 0 + et.random_complex(rng, σ) for σ in σ_vector]
actual_S3 = [-1 + et.random_complex(rng, σ) for σ in σ_vector]
actual_O3 = [+1 + et.random_complex(rng, σ) for σ in σ_vector]
actual_M4 = [ 0 + et.random_complex(rng, σ) for σ in σ_vector]
actual_M5 = [ 0 + et.random_complex(rng, σ) for σ in σ_vector]
M1a = VectorParameter(calset, f_vector, actual_M1)
S1a = VectorParameter(calset, f_vector, actual_S1)
O1a = VectorParameter(calset, f_vector, actual_O1)
M2a = VectorParameter(calset, f_vector, actual_M2)
S2a = VectorParameter(calset, f_vector, actual_S2)
O2a = VectorParameter(calset, f_vector, actual_O2)
M3a = VectorParameter(calset, f_vector, actual_M3)
S3a = VectorParameter(calset, f_vector, actual_S3)
O3a = VectorParameter(calset, f_vector, actual_O3)
M4a = VectorParameter(calset, f_vector, actual_M4)
M5a = VectorParameter(calset, f_vector, actual_M5)

#
# Find the actual through standard, modelling it as a length of
# transmission line with both metal and dielectric losses.
#
length = 0.01			# length in meters
vf = 2.0 / 3.0			# velocity factor
lm = 2.50e-4			# metal loss Np/m/Hz**(1/2)
ld = 5.86e-9			# dielectric loss (Np/m/Hz)
line_loss = 0.02		# dB per mm per sqrt(f_GHz)
mm_per_m = 1000.0		# mm per meter
np_per_db = 0.11512925		# neper per dB
c = 2.9979246e+8		# speed of light (m/s)
gl = [(lm * math.sqrt(f) + ld * f + 2.0j * math.pi * f / (vf * c)) * length
      for f in f_vector]
γ_through_actual = np.exp(-np.asarray(gl))
T_actual = VectorParameter(calset, f_vector, γ_through_actual)

#
# Create an NPData object for saving the measurements to files.
#
npd = NPData(ptype=PType.S, frequencies=points, rows=2, columns=2)
npd.frequency_vector = f_vector
%]
%O T16-EM-calibrate.py
%############################ begin calibration ###############################
%# Imports for calibration
from libvna.cal import (Calset, CalType, Solver, CorrelatedParameter,
	                UnknownParameter)
from libvna.data import NPData, PType
import math
import numpy as np


#
# Set the calibration frequency points.
#
fmin = %{fmin:7.1e%}
fmax = %{fmax:7.1e%}
f_vector = np.linspace(fmin, fmax, num=%{points%})

#
# Set up libvna.cal's error term solver.
#
calset = Calset()
solver = Solver(calset, CalType.T16, rows=2, columns=2,
                frequency_vector=f_vector)

#
# Set the measurement error.  First column is the noise floor of the
# receiver with no input signal; second column is noise in the generator,
# proportional to the received amplitude.
#
%[
print('m_error = np.asarray([', file=_file)
for i in range(len(f_vector)):
    c = ',' if i < len(f_vector) - 1 else ''
    print(f'    ({nf_vector[i]:9.7f}, {tr_vector[i]:9.7f}){c}', file=_file)
print('])', file=_file)
%]
solver.set_m_error(f_vector, m_error[:, 0], m_error[:, 1])
solver.et_tolerance = 1.0e-5	# set 10x smaller than error
solver.p_tolerance = 1.0e-5

#
# Get the (measured) connection non-repeatability error for each frequency
# and standard.  For simplicity, we assume the error is the same for all
# our calibration standards; in practice, it's probably not.  A better
# model would use different sigma values for each standard.
#
%[
print('σ_vector = np.asarray([', file=_file)
for i in range(len(f_vector)):
    print(f'    {σ_vector[i]:9.7f},', file=_file)
print('])', file=_file)
%]

#
# Create correlated parameters for each connection of the reflect
# standards.  Here, we assume that the standards are perfect on average.
# The second parameter can be a VectorParameter instead of a constant
# if we have a more trusted model of the standards.
#
M1 = CorrelatedParameter(calset,  0, f_vector, σ_vector)
S1 = CorrelatedParameter(calset, -1, f_vector, σ_vector)
O1 = CorrelatedParameter(calset, +1, f_vector, σ_vector)
M2 = CorrelatedParameter(calset,  0, f_vector, σ_vector)
S2 = CorrelatedParameter(calset, -1, f_vector, σ_vector)
O2 = CorrelatedParameter(calset, +1, f_vector, σ_vector)
M3 = CorrelatedParameter(calset,  0, f_vector, σ_vector)
S3 = CorrelatedParameter(calset, -1, f_vector, σ_vector)
O3 = CorrelatedParameter(calset, +1, f_vector, σ_vector)

#
# Create two impedance match errors and an unknown parameter for the
# through standard.  If the phase shift in the through remains safely
# below 90 degrees, we can give 1 as the estimate.  For a larger phase
# shift, we would have to provide a VectorParameter instead of 1 as the
# initial guess and give a more accurate estimate for each frequency.
#
M4 = CorrelatedParameter(calset, 0, f_vector, σ_vector)
M5 = CorrelatedParameter(calset, 0, f_vector, σ_vector)
T = UnknownParameter(calset, 1)

#
# Start with match on port 1 and short on port 2 (M-S).
#
%[
s = [[M1a, 0],
     [0, S1a]]
m = eterms.evaluate(calset, f_vector, s)
npd.data_array = m
npd.save('T16-EM-MS1.s2p')
%]
m = NPData(ptype=PType.S, filename='T16-EM-MS1.s2p')
solver.add_double_reflect(m.data_array, s11=M1, s22=S1)

#
# Change port 1 to open (O-S).
#
%[
s = [[O1a, 0],
     [0, S1a]]
m = eterms.evaluate(calset, f_vector, s)
npd.data_array = m
npd.save('T16-EM-OS2.s2p')
%]
m = NPData(ptype=PType.S, filename='T16-EM-OS2.s2p')
solver.add_double_reflect(m.data_array, s11=O1, s22=S1)

#
# Change port 2 to match (O-M).
#
%[
s = [[O1a, 0],
     [0, M2a]]
m = eterms.evaluate(calset, f_vector, s)
npd.data_array = m
npd.save('T16-EM-OM3.s2p')
%]
m = NPData(ptype=PType.S, filename='T16-EM-OM3.s2p')
solver.add_double_reflect(m.data_array, s11=O1, s22=M2)

#
# Change port 1 to short (S-M).
#
%[
s = [[S2a, 0],
     [0, M2a]]
m = eterms.evaluate(calset, f_vector, s)
npd.data_array = m
npd.save('T16-EM-SM4.s2p')
%]
m = NPData(ptype=PType.S, filename='T16-EM-SM4.s2p')
solver.add_double_reflect(m.data_array, s11=S2, s22=M2)

#
# Change port 2 to open (S-O).
#
%[
s = [[S2a, 0],
     [0, O2a]]
m = eterms.evaluate(calset, f_vector, s)
npd.data_array = m
npd.save('T16-EM-SO5.s2p')
%]
m = NPData(ptype=PType.S, filename='T16-EM-SO5.s2p')
solver.add_double_reflect(m.data_array, s11=S2, s22=O2)

#
# Change port 1 to match (M-O).
#
%[
s = [[M3a, 0],
     [0, O2a]]
m = eterms.evaluate(calset, f_vector, s)
npd.data_array = m
npd.save('T16-EM-MO6.s2p')
%]
m = NPData(ptype=PType.S, filename='T16-EM-MO6.s2p')
solver.add_double_reflect(m.data_array, s11=M3, s22=O2)

#
# Change port 2 to short (M-S).
#
%[
s = [[M3a, 0],
     [0, S3a]]
m = eterms.evaluate(calset, f_vector, s)
npd.data_array = m
npd.save('T16-EM-MS7.s2p')
%]
m = NPData(ptype=PType.S, filename='T16-EM-MS7.s2p')
solver.add_double_reflect(m.data_array, s11=M3, s22=S3)

#
# Change port 1 to open (O-S).
#
%[
s = [[O3a, 0],
     [0, S3a]]
m = eterms.evaluate(calset, f_vector, s)
npd.data_array = m
npd.save('T16-EM-OS8.s2p')
%]
m = NPData(ptype=PType.S, filename='T16-EM-OS8.s2p')
solver.add_double_reflect(m.data_array, s11=O3, s22=S3)

#
# Add measurement of the through standard.
#
%[
s = [[M4a, T_actual],
     [T_actual, M5a]]
m = eterms.evaluate(calset, f_vector, s)
npd.data_array = m
npd.save('T16-EM-T.s2p')
%]
m = NPData(ptype=PType.S, filename='T16-EM-T.s2p')
solver.add_line(m.data_array, [[M4, T], [T, M5]])

#
# Solve, add to Calset and save.
#
solver.solve()
solver.add_to_calset('mycal')
calset.save('T16-EM.vnacal')
%############################# end calibration ################################
%#
%#
%O T16-EM-apply.py
%############################# begin application ##############################
%# Imports for calibration
from libvna.cal import Calset


# Load the calibration from file.
calset = Calset('T16-EM.vnacal')
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
npd.save('T16-EM-expected.s2p')
%]
# Measured response
%[
s = [[VectorParameter(calset, f_vector, expected[:, 0, 0]),
      VectorParameter(calset, f_vector, expected[:, 0, 1])],
     [VectorParameter(calset, f_vector, expected[:, 1, 0]),
      VectorParameter(calset, f_vector, expected[:, 1, 1])]]
m = eterms.evaluate(calset, f_vector, s)
npd.data_array = m
npd.save('T16-EM-measured.s2p')
et.print_matrix(m, file=_file, indent=_indent)
%]

# Apply calibration and save in Touchstone format.
result = calibration.apply(f_vector, m)
result.save('T16-EM-corrected.s2p')
%############################## end application ###############################
