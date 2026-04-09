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
fmin = 1.0e+6
fmax = 1.0e+9
points = 10
f_vector = np.logspace(np.log10(fmin), np.log10(fmax), num=points)

# Generate random error terms
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
T_actual = calset.vector_parameter(f_vector, γ_through_actual)

calset.parameter_matrix([[0, T_actual], [T_actual, 0]]) \
    .to_npdata(f_vector).save("UT-Ta.s2p")
%]
%O UT-calibrate.py
%############################ begin calibration ###############################
%# Imports for calibration
from libvna.cal import Calset, CalType, Solver
from libvna.data import NPData, PType
import math
import numpy as np


# Set the calibration frequency points.
fmin = %{fmin:7.1e%}
fmax = %{fmax:7.1e%}
f_vector = np.logspace(np.log10(fmin), np.log10(fmax), num=%{points%})

# Create the calibration container and error term solver.
calset = Calset()
solver = Solver(calset, CalType.TE10, rows=2, columns=2,
                frequency_vector=f_vector)

# Define the short, open and load standards from our calibration kit.
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

# Create an unknown parameter for our through standard.  As long as the
# phase shift in the through remains safely below 90 degrees, we can
# give 1 as the estimate.  For a larger phase shift, we would need to
# give a more accurate estimate for each frequency.
T = calset.unknown_parameter(1)

# Add measurement of the short-open standard.
%[
s = [[short_std, 0],
     [0, open_std]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_double_reflect(m, s11=short_std, s22=open_std)

# Add measurement of the open-load standard.
%[
s = [[open_std, 0],
     [0, load_std]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_double_reflect(m, s11=open_std, s22=load_std)

# Add measurement of the load-short standard.
%[
s = [[load_std, 0],
     [0, short_std]]
m = eterms.evaluate(calset, f_vector, s)
et.print_matrix(m, file=_file, indent=_indent)
%]
solver.add_double_reflect(m, s11=load_std, s22=short_std)


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
calset.parameter_matrix([[0, T], [T, 0]]).to_npdata(f_vector).save("UT-T.s2p")
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
s = [[calset.vector_parameter(f_vector, expected[:, 0, 0]),
      calset.vector_parameter(f_vector, expected[:, 0, 1])],
     [calset.vector_parameter(f_vector, expected[:, 1, 0]),
      calset.vector_parameter(f_vector, expected[:, 1, 1])]]
m = eterms.evaluate(calset, f_vector, s)
npd.data_array = m
npd.save('UT-measured.s2p')
et.print_matrix(m, file=_file, indent=_indent)
%]

# Apply calibration and save in Touchstone format.
result = calibration.apply(f_vector, m)
result.save('UT-corrected.s2p')
%############################## end application ###############################
