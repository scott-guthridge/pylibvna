#!/usr/bin/python3
"""
Test the libvna.data module.
"""
from libvna.cal import CalType, Calset, Solver, VectorParameter
from libvna.conv import ytos, ztos
from libvna.data import NPData, PType
import math
import numpy as np
import os
import unittest
import sys

sys.path.append('examples')
sys.path.append('../examples')
import random_error_terms as ret

rng = np.random.default_rng(seed=7)

fmin = 1.0e+9
fmax = 10.0e+9
points = 10
f_vector = np.linspace(fmin, fmax, num=points)


class TestModule(unittest.TestCase):
    def test_SOLT(self):
        eterms = ret.RandomErrorTerms(rng, CalType.E12, 2, 1, fmin, fmax)
        calset = Calset()
        solver = Solver(calset, CalType.E12, rows=2, columns=1,
                        frequency_vector=f_vector)
        delay_short   = np.random.normal() / fmax
        delay_open    = np.random.normal() / fmax
        delay_match   = np.random.normal() / fmax
        delay_through = np.random.normal() / fmax
        delay_dut1    = np.random.normal() / fmax
        delay_dut2    = np.random.normal() / fmax

        r_short    = np.exp(-4.0j * math.pi * delay_short   * f_vector)
        r_open     = np.exp(-4.0j * math.pi * delay_open    * f_vector)
        r_match    = np.exp(-4.0j * math.pi * delay_match   * f_vector)
        r_through  = np.exp(-2.0j * math.pi * delay_through * f_vector)
        r_dut1     = np.exp(-2.0j * math.pi * delay_dut1    * f_vector)
        r_dut2     = np.exp(-2.0j * math.pi * delay_dut2    * f_vector)

        # Add short standard.
        s11 = VectorParameter(calset, f_vector, r_short * -1.0)
        s = [[s11, 0.0],
             [0.0, 0.0]]
        m = eterms.evaluate(calset, f_vector, s)
        solver.add_single_reflect(m, -1.0, delay=delay_short)

        # Add open standard.
        s11 = VectorParameter(calset, f_vector, r_open * 1.0)
        s = [[s11, 0.0],
             [0.0, 0.0]]
        m = eterms.evaluate(calset, f_vector, s)
        solver.add_single_reflect(m, 1.0, delay=delay_open)

        # Add match standard.
        s11 = 0.0
        s = [[s11, 0.0],
             [0.0, 0.0]]
        m = eterms.evaluate(calset, f_vector, s)
        solver.add_single_reflect(m, 0.0, delay=delay_match)

        # Add through standard.
        s12 = VectorParameter(calset, f_vector, r_through)
        s = [[0.0, s12],
             [s12, 0.0]]
        m = eterms.evaluate(calset, f_vector, s)
        solver.add_through(m, delay=delay_through)

        solver.solve()
        solver.add_to_calset('cal')

        # Add DUT as LC LPF
        l = 2.5e-9
        c = 1.0e-12
        expected = np.empty((points, 2, 2), dtype=complex)
        delayed = np.empty((points, 2, 2), dtype=complex)
        for i, f in enumerate(f_vector):
            jω = 2.0j * math.pi * f
            z1 = jω * l
            z2 = 1.0 / (jω * c)
            z = [[z1 + z2, z2],
                 [z2,      z2]]
            s = ztos(z)
            expected[i, ...] = s
            r1 = r_dut1[i]
            r2 = r_dut2[i]
            delayed[i, ...] = s * [[r1 * r1, r1 * r2],
                                   [r1 * r2, r2 * r2]]
        s = [[VectorParameter(calset, f_vector, delayed[:, 0, 0]),
              VectorParameter(calset, f_vector, delayed[:, 0, 1])],
             [VectorParameter(calset, f_vector, delayed[:, 1, 0]),
              VectorParameter(calset, f_vector, delayed[:, 1, 1])]]
        m1 = eterms.evaluate(calset, f_vector, s)
        m2 = eterms.evaluate(calset, f_vector, np.flipud(np.fliplr(s)))
        m = np.concatenate((m1, np.flip(m2, axis=1)), axis=2)
        calibration = calset.calibrations[0]
        result = calibration.apply(f_vector, m,
                                   delay_vector=[delay_dut1, delay_dut2])
        self.assertTrue(np.allclose(result.data_array, expected,
                                    rtol=1.0e-5, atol=1.0e-5))

    def test_TE10(self):
        eterms = ret.RandomErrorTerms(rng, CalType.TE10, 2, 2, fmin, fmax)
        calset = Calset()
        solver = Solver(calset, CalType.TE10, rows=2, columns=2,
                        frequency_vector=f_vector)
        delay_short   = np.random.normal() / fmax
        delay_open    = np.random.normal() / fmax
        delay_match   = np.random.normal() / fmax
        delay_through = np.random.normal() / fmax
        delay_dut1    = np.random.normal() / fmax
        delay_dut2    = np.random.normal() / fmax

        r_short    = np.exp(-4.0j * math.pi * delay_short   * f_vector)
        r_open     = np.exp(-4.0j * math.pi * delay_open    * f_vector)
        r_match    = np.exp(-4.0j * math.pi * delay_match   * f_vector)
        r_through  = np.exp(-2.0j * math.pi * delay_through * f_vector)
        r_dut1     = np.exp(-2.0j * math.pi * delay_dut1    * f_vector)
        r_dut2     = np.exp(-2.0j * math.pi * delay_dut2    * f_vector)

        # Add short-open standard.
        s11 = VectorParameter(calset, f_vector, r_short * -1.0)
        s22 = VectorParameter(calset, f_vector, r_open  *  1.0)
        s = [[s11, 0.0],
             [0.0, s22]]
        m = eterms.evaluate(calset, f_vector, s)
        solver.add_double_reflect(m, -1.0, 1.0,
                                  delay1=delay_short, delay2=delay_open)

        # Add match-short standard.
        s11 = 0
        s22 = VectorParameter(calset, f_vector, r_short * -1.0)
        s = [[s11, 0.0],
             [0.0, s22]]
        m = eterms.evaluate(calset, f_vector, s)
        solver.add_double_reflect(m, 0.0, -1.0,
                                  delay1=delay_match, delay2=delay_short)

        # Add through standard (using line).
        s12 = VectorParameter(calset, f_vector, r_through)
        s = [[0.0, s12],
             [s12, 0.0]]
        m = eterms.evaluate(calset, f_vector, s)
        solver.add_line(m, [[0.0, 1.0], [1.0, 0.0]],
                        delay2=delay_through)

        solver.solve()
        solver.add_to_calset('cal')

        # Add DUT as LC HPF
        l = 2.5e-9
        c = 1.0e-12
        expected = np.empty((points, 2, 2), dtype=complex)
        delayed = np.empty((points, 2, 2), dtype=complex)
        for i, f in enumerate(f_vector):
            jω = 2.0j * math.pi * f
            z1 = 1.0 / (jω * c)
            z2 = jω * l
            z = [[z1 + z2, z2],
                 [z2,      z2]]
            s = ztos(z)
            expected[i, ...] = s
            r1 = r_dut1[i]
            r2 = r_dut2[i]
            delayed[i, ...] = s * [[r1 * r1, r1 * r2],
                                   [r1 * r2, r2 * r2]]
        s = [[VectorParameter(calset, f_vector, delayed[:, 0, 0]),
              VectorParameter(calset, f_vector, delayed[:, 0, 1])],
             [VectorParameter(calset, f_vector, delayed[:, 1, 0]),
              VectorParameter(calset, f_vector, delayed[:, 1, 1])]]
        m = eterms.evaluate(calset, f_vector, s)
        calibration = calset.calibrations[0]
        result = calibration.apply(f_vector, m,
                                   delay_vector=[delay_dut1, delay_dut2])
        self.assertTrue(np.allclose(result.data_array, expected,
                                    rtol=1.0e-5, atol=1.0e-5))

    def test_T16_3x3(self):
        eterms = ret.RandomErrorTerms(rng, CalType.T16, 3, 3, fmin, fmax)
        calset = Calset()
        solver = Solver(calset, CalType.T16, rows=3, columns=3,
                        frequency_vector=f_vector)
        delay11    = np.random.normal() / fmax
        delay12    = np.random.normal() / fmax
        delay13    = np.random.normal() / fmax
        delay21    = np.random.normal() / fmax
        delay22    = np.random.normal() / fmax
        delay23    = np.random.normal() / fmax
        delay31    = np.random.normal() / fmax
        delay32    = np.random.normal() / fmax
        delay33    = np.random.normal() / fmax
        delay41    = np.random.normal() / fmax
        delay42    = np.random.normal() / fmax
        delay43    = np.random.normal() / fmax
        delay51    = np.random.normal() / fmax
        delay52    = np.random.normal() / fmax
        delay53    = np.random.normal() / fmax
        delay_dut1 = np.random.normal() / fmax
        delay_dut2 = np.random.normal() / fmax
        delay_dut3 = np.random.normal() / fmax

        r11 = np.exp(-2.0j * math.pi * delay11    * f_vector)
        r12 = np.exp(-2.0j * math.pi * delay12    * f_vector)
        r13 = np.exp(-2.0j * math.pi * delay13    * f_vector)
        r21 = np.exp(-2.0j * math.pi * delay21    * f_vector)
        r22 = np.exp(-2.0j * math.pi * delay22    * f_vector)
        r23 = np.exp(-2.0j * math.pi * delay23    * f_vector)
        r31 = np.exp(-2.0j * math.pi * delay31    * f_vector)
        r32 = np.exp(-2.0j * math.pi * delay32    * f_vector)
        r33 = np.exp(-2.0j * math.pi * delay33    * f_vector)
        r41 = np.exp(-2.0j * math.pi * delay41    * f_vector)
        r42 = np.exp(-2.0j * math.pi * delay42    * f_vector)
        r43 = np.exp(-2.0j * math.pi * delay43    * f_vector)
        r51 = np.exp(-2.0j * math.pi * delay51    * f_vector)
        r52 = np.exp(-2.0j * math.pi * delay52    * f_vector)
        r53 = np.exp(-2.0j * math.pi * delay53    * f_vector)
        r_dut1 = np.exp(-2.0j * math.pi * delay_dut1 * f_vector)
        r_dut2 = np.exp(-2.0j * math.pi * delay_dut2 * f_vector)
        r_dut3 = np.exp(-2.0j * math.pi * delay_dut3 * f_vector)

        #
        # Add random 3-port standard 1.
        #
        s_matrix = np.asarray([
            [+0.4034-0.2752j, +0.7065+1.0564j, -0.0052-0.2987j],
            [-0.3819-0.2956j, -0.0749-0.5035j, -0.8099-2.1261j],
            [+0.1404-1.8447j, +0.4446-0.9603j, -0.3800-0.3787j]
        ])
        s_delayed = np.empty((3, 3), dtype=object)
        r = [r11, r12, r13]
        for i in range(3):
            for j in range(3):
                s_delayed[i, j] = VectorParameter(calset, f_vector,
                                                  r[i] * r[j] * s_matrix[i, j])
        m = eterms.evaluate(calset, f_vector, s_delayed)
        # map = [2, 3, 1]
        map = [1, 2, 3]
        inverse_indices = [None, None, None]
        for i, p in enumerate(map):
            inverse_indices[p-1] = i
        s_permuted = s_matrix[inverse_indices, :][:, inverse_indices]
        assert np.allclose(s_matrix, s_permuted)
        solver.add_mapped_matrix(m, s_permuted,
                                 delay_vector=[delay11, delay12, delay13],
                                 port_map=map)

        #
        # Add random 3-port standard 2.
        #
        s_matrix = np.asarray([
            [+0.8030+0.5993j, +0.8734+0.2789j, -0.3054-0.5189j],
            [+0.9983+0.8756j, +0.4573-2.0263j, -2.4729+0.1327j],
            [-2.0109+0.3117j, -0.9881-0.5521j, -1.9446+1.3456j]
        ])
        s_delayed = np.empty((3, 3), dtype=object)
        r = [r21, r22, r23]
        for i in range(3):
            for j in range(3):
                s_delayed[i, j] = VectorParameter(calset, f_vector,
                                                  r[i] * r[j] * s_matrix[i, j])
        m = eterms.evaluate(calset, f_vector, s_delayed)
        # map = [3, 1, 2]
        map = [1, 2, 3]
        inverse_indices = [None, None, None]
        for i, p in enumerate(map):
            inverse_indices[p-1] = i
        s_permuted = s_matrix[inverse_indices, :][:, inverse_indices]
        assert np.allclose(s_matrix, s_permuted)
        solver.add_mapped_matrix(m, s_permuted,
                                 delay_vector=[delay21, delay22, delay23],
                                 port_map=map)

        #
        # Add random 3-port standard 3.
        #
        s_matrix = np.asarray([
            [-0.8374-0.4580j, -0.4390-1.7829j, +0.2250-0.5186j],
            [+0.4076+0.0577j, +0.7048+1.2195j, -1.2111+0.2743j],
            [-0.4238-0.4229j, -0.7292-1.0598j, -0.3230+0.0595j]
        ])
        s_delayed = np.empty((3, 3), dtype=object)
        r = [r31, r32, r33]
        for i in range(3):
            for j in range(3):
                s_delayed[i, j] = VectorParameter(calset, f_vector,
                                                  r[i] * r[j] * s_matrix[i, j])
        m = eterms.evaluate(calset, f_vector, s_delayed)
        # map = [3, 2, 1]
        map = [1, 2, 3]
        inverse_indices = [None, None, None]
        for i, p in enumerate(map):
            inverse_indices[p-1] = i
        s_permuted = s_matrix[inverse_indices, :][:, inverse_indices]
        assert np.allclose(s_matrix, s_permuted)
        solver.add_mapped_matrix(m, s_permuted,
                                 delay_vector=[delay31, delay32, delay33],
                                 port_map=map)

        #
        # Add random 3-port standard 4.
        #
        s_matrix = np.asarray([
            [-0.6309+0.1204j, +0.4045-0.0695j, -0.9811-0.8745j],
            [-1.1644+0.0586j, +0.3026+0.6481j, -1.2244-1.0485j],
            [+0.3459-1.6778j, -0.2581-0.1265j, -0.4891+0.2133j]
        ])
        s_delayed = np.empty((3, 3), dtype=object)
        r = [r41, r42, r43]
        for i in range(3):
            for j in range(3):
                s_delayed[i, j] = VectorParameter(calset, f_vector,
                                                  r[i] * r[j] * s_matrix[i, j])
        m = eterms.evaluate(calset, f_vector, s_delayed)
        # map = [2, 1, 3]
        map = [1, 2, 3]
        inverse_indices = [None, None, None]
        for i, p in enumerate(map):
            inverse_indices[p-1] = i
        s_permuted = s_matrix[inverse_indices, :][:, inverse_indices]
        assert np.allclose(s_matrix, s_permuted)
        solver.add_mapped_matrix(m, s_permuted,
                                 delay_vector=[delay41, delay42, delay43],
                                 port_map=map)
        solver.solve()
        solver.add_to_calset('cal')

        #
        # Add random 3-port standard 5.
        #
        s_matrix = np.asarray([
            [+0.8902-1.0414j, -0.6846+0.5228j, -0.5762-0.6446j],
            [-0.1297+0.1755j, +1.8558-0.2308j, -0.5737-2.2328j],
            [+0.7810-0.3430j, -0.1914-0.0231j, +0.6078+0.2086j]
        ])
        s_delayed = np.empty((3, 3), dtype=object)
        r = [r51, r52, r53]
        for i in range(3):
            for j in range(3):
                s_delayed[i, j] = VectorParameter(calset, f_vector,
                                                  r[i] * r[j] * s_matrix[i, j])
        m = eterms.evaluate(calset, f_vector, s_delayed)
        # map = [2, 1, 3]
        map = [1, 2, 3]
        inverse_indices = [None, None, None]
        for i, p in enumerate(map):
            inverse_indices[p-1] = i
        s_permuted = s_matrix[inverse_indices, :][:, inverse_indices]
        assert np.allclose(s_matrix, s_permuted)
        solver.add_mapped_matrix(m, s_permuted,
                                 delay_vector=[delay51, delay52, delay53],
                                 port_map=map)
        solver.solve()
        solver.add_to_calset('cal')
        calset.save('mycal.vnacal')

        #
        # Add DUT.  Our DUT is a delta configuration with:
        #     capacitor, c, between ports 1 and 2
        #     inductor, l, between ports 2 and 3
        #     resistor, r, between ports 3 and 1
        #
        l = 2.5e-9
        c = 1.0e-12
        r = 50
        g = 1.0 / r
        expected = np.empty((points, 3, 3), dtype=complex)
        delayed = np.empty((points, 3, 3), dtype=complex)
        for i, f in enumerate(f_vector):
            jω = 2.0j * math.pi * f
            yl = 1.0 / (jω * l)
            yc = jω * c
            yr = g
            y = [
                [yc + yr, -yc, -yr],
                [-yc, yc + yl, -yl],
                [-yr, -yl, yl + yr]
            ]
            s = ytos(y)
            expected[i, ...] = s
            r1 = r_dut1[i]
            r2 = r_dut2[i]
            r3 = r_dut3[i]
            delayed[i, ...] = s * [[r1 * r1, r1 * r2, r1 * r3],
                                   [r2 * r1, r2 * r2, r2 * r3],
                                   [r3 * r1, r3 * r2, r3 * r3]]

        s_delayed = np.empty((3, 3), dtype=object)
        for i in range(3):
            for j in range(3):
                s_delayed[i, j] = VectorParameter(calset, f_vector,
                                                  delayed[:, i, j])
        m = eterms.evaluate(calset, f_vector, s_delayed)
        calibration = calset.calibrations[0]
        result = calibration.apply(f_vector, m, delay_vector=[delay_dut1,
                                                              delay_dut2,
                                                              delay_dut3])
        self.assertTrue(np.allclose(result.data_array, expected,
                                    rtol=1.0e-5, atol=1.0e-5))


if __name__ == '__main__':
    unittest.main()
