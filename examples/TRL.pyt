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
rng = np.random.default_rng(seed=7)

# Set frequency range and number of points
fmin = 1.0e+9
fmax = 8.0e+9
points = 10
f_vector = np.linspace(fmin, fmax, num=points)

# Generate random error terms
eterms = et.RandomErrorTerms(rng, CalType.TE10, 2, 2, fmin, fmax)
calset = Calset()

# The actual reflect standard is 5 ohms in series with 700pH.  Calculate
# the reflection coefficient.
r_vector = []
for i, f in enumerate(f_vector):
    jω = 2.0j * math.pi * f
    zr = 5.0 + 700e-12 * jω 
    r_vector.append(ztos([[zr]])[0, 0])
R_actual = calset.vector_parameter(f_vector, r_vector)

# For the line standard, the effective permittivity is less than we
# estimate, and there is some loss.
line_length = 0.01		# length in meters
εr_eff = 2.50			# actual effective permittivity
vf = 1 / math.sqrt(εr_eff)	# velocity factor
c = 2.9979246e+8		# speed of light (m/s)
line_loss = 0.02		# dB per mm per sqrt(f_GHz)
mm_per_m = 1000.0		# mm per meter
np_per_db = 0.11512925		# neper to dB
gl = [(np_per_db + mm_per_m * line_loss * math.sqrt(f / 1.0e+9)
      + 2.0j * math.pi / (c * vf) * f) * line_length for f in f_vector]
γ_line_actual = np.exp(-np.asarray(gl))
L_actual = calset.vector_parameter(f_vector, γ_line_actual)

calset.parameter_matrix([[R_actual]]).to_npdata(f_vector).save("TRL-Ra.s1p")
calset.parameter_matrix([[0, L_actual], [L_actual, 0]]) \
    .to_npdata(f_vector).save("TRL-La.s2p")
%]
%O TRL-calibrate.py
%############################ begin calibration ###############################
%# Imports for calibration
from libvna.cal import Calset, CalType, Solver
from libvna.data import NPData, PType
import math
import numpy as np


# Set the calibration frequency points.
fmin = %{fmin:7.1e%}
fmax = %{fmax:7.1e%}
f_vector = np.linspace(fmin, fmax, num=%{points%})

# Create the calibration container and error term solver.
calset = Calset()
solver = Solver(calset, CalType.TE10, rows=2, columns=2,
                frequency_vector=f_vector)

# Our line is a 1cm long trace of microstrip on FR4 at 50 ohms.
# We estimate εr_eff to be 2.75.  From these, we can estimate the
# transmission parameter.
#
line_length = 0.01		# length in meters
εr_eff = 2.75			# effective permittivity
vf = 1 / math.sqrt(εr_eff)	# velocity factor
c = 2.9979246e+8		# speed of light (m/s)
line_estimated = np.exp(-2.0j * math.pi / (c * vf) * line_length * f_vector)

# Create unknown parameters for the reflect and line, giving
# estimated values for each.
R = calset.unknown_parameter(-1)
L = calset.unknown_parameter((f_vector, line_estimated))

# Add measurement of the through standard.
%[
s = [[0, 1],
     [1, 0]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_through(m)

# Add measurement of the unknown reflect standard.
%[
s = [[R_actual, 0],
     [0, R_actual]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_double_reflect(m, s11=R, s22=R)

# Add measurement of the unknown line standard.
%[
s = [[0, L_actual],
     [L_actual, 0]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_line(m, [[0, L], [L, 0]])

# Solve, add to Calset and save.
solver.solve()
solver.add_to_calset('mycal')
calset.save('TRL.vnacal')

# Save the solved R
calset.parameter_matrix([[R]]).to_npdata(f_vector).save("TRL-R.s1p")

# Save the solved L
calset.parameter_matrix([[0, L], [L, 0]]).to_npdata(f_vector).save("TRL-L.s2p")
%############################# end calibration ################################
%#
%#
%O TRL-apply.py
%############################# begin application ##############################
%# Imports for calibration
from libvna.cal import Calset


# Load the calibration from file.
calset = Calset('TRL.vnacal')
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
npd.save('TRL-expected.s2p')
%]
# Measured response
%[
s = [[calset.vector_parameter(f_vector, expected[:, 0, 0]),
      calset.vector_parameter(f_vector, expected[:, 0, 1])],
     [calset.vector_parameter(f_vector, expected[:, 1, 0]),
      calset.vector_parameter(f_vector, expected[:, 1, 1])]]
m = eterms.evaluate(calset, f_vector, s)
npd.data_array = m
npd.save('TRL-measured.s2p')
et.print_matrix(m, file=_file, indent=_indent)
%]

# Apply calibration and save in Touchstone format.
result = calibration.apply(f_vector, m)
result.save('TRL-corrected.s2p')
%############################## end application ###############################
