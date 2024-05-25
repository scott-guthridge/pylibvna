#!/usr/bin/python3
from libvna.cal import (Calset, CalType, Solver, CorrelatedParameter,
                        ScalarParameter, UnknownParameter)
from libvna.data import NPData, PType
import math
import matplotlib.pyplot as plt
import numpy as np
import random_error_terms as et
from scipy.ndimage import gaussian_filter1d
import scipy.stats as stats
import sys
import yaml


Trials = 10000
noise_floor = 0.002
tracking_error = 0.02
connection_error = 0.02

#
# Set up random number generator.
#
rng = np.random.default_rng(seed=3)

remove_spaces = np.vectorize(lambda s: s.replace(' ', ''))
change_none_to_zero = np.vectorize(lambda v: 0.0j if v is None else v)

def str_to_complex_matrix(s):
    temp = np.asarray(s, dtype=str)
    return np.asarray(remove_spaces(temp), dtype=complex)


def run_trial(rng, nf_error, tr_error, conn_error, model_errors):
    frequency = 1.0e+9      # arbitrary test frequency
    f_vector = [frequency]
    nf_vector = [nf_error]
    tr_vector = [tr_error]
    eterms = et.RandomErrorTerms(rng, CalType.T16, 2, 2, frequency, frequency,
                                 nf_vector=nf_vector, tr_vector=tr_vector)
    calset = Calset()
    solver = Solver(calset, CalType.T16, rows=2, columns=2,
                    frequency_vector=f_vector)
    if model_errors:
        solver.set_m_error(f_vector, nf_vector, tr_vector)
        tolerance = nf_error / 10.0     # noise floor plus 1 digit
        solver.et_tolerance = tolerance
        solver.p_tolerance  = tolerance

    # Form the actual calibration parameters.
    L1a = ScalarParameter(calset, et.random_complex(rng, conn_error))
    S1a = ScalarParameter(calset, et.random_complex(rng, conn_error) - 1)
    O1a = ScalarParameter(calset, et.random_complex(rng, conn_error) + 1)
    L2a = ScalarParameter(calset, et.random_complex(rng, conn_error))
    S2a = ScalarParameter(calset, et.random_complex(rng, conn_error) - 1)
    O2a = ScalarParameter(calset, et.random_complex(rng, conn_error) + 1)
    L3a = ScalarParameter(calset, et.random_complex(rng, conn_error))
    S3a = ScalarParameter(calset, et.random_complex(rng, conn_error) - 1)
    O3a = ScalarParameter(calset, et.random_complex(rng, conn_error) + 1)
    L4a = ScalarParameter(calset, et.random_complex(rng, conn_error))
    L5a = ScalarParameter(calset, et.random_complex(rng, conn_error))
    t_actual = (rng.uniform(0.75, 1) *
                1j ** (rng.uniform(-0.2, 0.2) * math.pi))
    Ta  = ScalarParameter(calset, t_actual)

    # Form the operational calibration parameters.
    if model_errors:
        L1 = CorrelatedParameter(calset, +0, f_vector, [conn_error])
        S1 = CorrelatedParameter(calset, -1, f_vector, [conn_error])
        O1 = CorrelatedParameter(calset, +1, f_vector, [conn_error])
        L2 = CorrelatedParameter(calset, +0, f_vector, [conn_error])
        S2 = CorrelatedParameter(calset, -1, f_vector, [conn_error])
        O2 = CorrelatedParameter(calset, +1, f_vector, [conn_error])
        L3 = CorrelatedParameter(calset, +0, f_vector, [conn_error])
        S3 = CorrelatedParameter(calset, -1, f_vector, [conn_error])
        O3 = CorrelatedParameter(calset, +1, f_vector, [conn_error])
        L4 = CorrelatedParameter(calset, +0, f_vector, [conn_error])
        L5 = CorrelatedParameter(calset, +0, f_vector, [conn_error])
    else:
        L1 = ScalarParameter(calset, +0)
        S1 = ScalarParameter(calset, -1)
        O1 = ScalarParameter(calset, +1)
        L2 = L1
        S2 = S1
        O2 = O1
        L3 = L1
        S3 = S1
        O3 = O1
        L4 = L1
        L5 = L1
    T = UnknownParameter(calset, 1)

    # Add measurement of the load-short standard.
    s = [[L1a, 0],
         [0, S1a]]
    m = eterms.evaluate(calset, f_vector, s)
    solver.add_double_reflect(m, s11=L1, s22=S1)

    # Change left side to open.
    s = [[O1a, 0],
         [0, S1a]]
    m = eterms.evaluate(calset, f_vector, s)
    solver.add_double_reflect(m, s11=O1, s22=S1)

    # Change right side to load.
    s = [[O1a, 0],
         [0, L2a]]
    m = eterms.evaluate(calset, f_vector, s)
    solver.add_double_reflect(m, s11=O1, s22=L2)

    # Change left side to short.
    s = [[S2a, 0],
         [0, L2a]]
    m = eterms.evaluate(calset, f_vector, s)
    solver.add_double_reflect(m, s11=S2, s22=L2)

    # Change right side to open.
    s = [[S2a, 0],
         [0, O2a]]
    m = eterms.evaluate(calset, f_vector, s)
    solver.add_double_reflect(m, s11=S2, s22=O2)

    # Change left side to load.
    s = [[L3a, 0],
         [0, O2a]]
    m = eterms.evaluate(calset, f_vector, s)
    solver.add_double_reflect(m, s11=L3, s22=O2)

    # Change right side to short.
    s = [[L3a, 0],
         [0, S3a]]
    m = eterms.evaluate(calset, f_vector, s)
    solver.add_double_reflect(m, s11=L3, s22=S3)

    # Change left side to open.
    s = [[O3a, 0],
         [0, S3a]]
    m = eterms.evaluate(calset, f_vector, s)
    solver.add_double_reflect(m, s11=O3, s22=S3)

    # Add measurement of the unknown through standard.
    s = [[L4a, Ta],
         [Ta, L5a]]
    m = eterms.evaluate(calset, f_vector, s)
    solver.add_line(m, [[L4, T], [T, L5]])

    # Solve, add to calset and save.
    solver.solve()
    solver.add_to_calset('mycal')
    calset.save('T16-m.vnacal')

    # Load the solved error parameters.
    with open('T16-m.vnacal', 'r') as vnacal_file:
        vnacal_file.readline() # discard
        vnacal = yaml.safe_load(vnacal_file)

    # Convert from T parameters to E parameters and normalize et11 to 1.
    data = vnacal['calibrations'][0]['data'][0]
    ts = str_to_complex_matrix(data['ts'])
    ti = str_to_complex_matrix(data['ti'])
    tx = str_to_complex_matrix(data['tx'])
    tm = str_to_complex_matrix(data['tm'])
    c_et = np.linalg.inv(tm)
    c_el = ti @ c_et
    c_er = ts - ti @ c_et @ tx
    c_em = -c_et @ tx
    k = c_et[0, 0]
    c_et /= k
    c_er *= k

    # Get the actual error terms and normalize et11 to 1.
    a_el, a_er, a_et, a_em = eterms.get_eterms(0, f_vector[0])
    k = a_et[0, 0]
    a_et /= k
    a_er *= k

    # Find the RMS difference between computed and actual error terms.
    d_el = c_el - a_el
    d_er = c_er - a_er
    d_et = c_et - a_et
    d_em = c_em - a_em
    s = np.sum((d_el * np.conj(d_el)).real)
    s += np.sum((d_er * np.conj(d_er)).real)
    s += np.sum((d_et * np.conj(d_et)).real)
    s += np.sum((d_em * np.conj(d_em)).real)
    s /= 15

    return math.sqrt(s)


#
# Run trials without and with measurement error modeling.
#
errors_without = []
for trial in range(Trials):
    try:
        errors_without.append(run_trial(rng, noise_floor, tracking_error,
                                        connection_error, False))
    except Exception as e:
        pass

errors_with = []
for trial in range(Trials):
    try:
        errors_with.append(run_trial(rng, noise_floor, tracking_error,
                                     connection_error, True))
    except Exception as e:
        pass

#
# Show the success rates and stop if either is too low.
#
success_rate1 = len(errors_without) / Trials
success_rate2 = len(errors_with) / Trials
if success_rate1 < 0.9:
    print('Invalid result: success without error modeling is only '
          f'{success_rate1 * 100:5.1f}%.')
    sys.exit(1)
if success_rate2 < 0.9:
    print('Invalid result: success with error modeling is only '
          f'{success_rate2 * 100:5.1f}%.')
    sys.exit(1)

#
# Bin, smooth and calculate CDF
#
bin_count = 1000
smoothness = 7
sorted1 = np.sort(errors_without)
sorted2 = np.sort(errors_with)
cut1 = int(0.01 * len(sorted1)) # discard the highest 1%
cut2 = int(0.01 * len(sorted2))
data1 = sorted1[:-cut1]
data2 = sorted2[:-cut2]
hist1, bin_edges1 = np.histogram(data1, bins=bin_count, density=True)
hist2, bin_edges2 = np.histogram(data2, bins=bin_count, density=True)
bin_centers1 = (bin_edges1[:-1] + bin_edges1[1:]) / 2
bin_centers2 = (bin_edges2[:-1] + bin_edges2[1:]) / 2
smoothed1 = gaussian_filter1d(hist1, sigma=smoothness)
smoothed2 = gaussian_filter1d(hist2, sigma=smoothness)
smoothed1 /= np.sum(smoothed1)
smoothed2 /= np.sum(smoothed2)
cdf1 = np.cumsum(smoothed1)
cdf2 = np.cumsum(smoothed2)

#
# Show statistics.
#
mean1 = np.mean(errors_without)
p25_1 = np.percentile(errors_without, 25)
median1 = np.median(errors_without)
p75_1 = np.percentile(errors_without, 75)
mean2 = np.mean(errors_with)
p25_2 = np.percentile(errors_with, 25)
median2 = np.median(errors_with)
p75_2 = np.percentile(errors_with, 75)
print('Error Modeling   Mean    25 Pctile  Median   75 Pctile Success Rate')
print('-------------- --------- --------- --------- --------- ------------')
print(f'      No       {mean1:9.5f} {p25_1:9.5f} {median1:9.5f} {p75_1:9.5f} '
      f'{success_rate1 * 100:11.2f}%')
print(f'      Yes      {mean2:9.5f} {p25_2:9.5f} {median2:9.5f} {p75_2:9.5f} '
      f'{success_rate2 * 100:11.2f}%')

#
# Plot the CDFs
#
plt.figure(figsize=(10, 6))
plt.plot(bin_centers2, cdf2, label='With Error Modeling')
plt.plot(bin_centers1, cdf1, label='Without Error Modeling')
plt.xlim(0, max(max(bin_centers1), max(bin_centers2)))
plt.ylim(0, 1)
plt.grid(True)
plt.xlabel('RMS Error of Calibration Error Terms')
plt.ylabel('Cumulative Probability')
plt.title('Cumulative Probability Distribution')
plt.legend()
plt.savefig('error-modeling-cdfs.png')
plt.show()
plt.close()
