#!/usr/bin/python3
"""
Test the libvna.data module.
"""
import os
import unittest
import numpy as np
import libvna.data as data

M_SQRT1_2 = 1.41421356237309504880

def crandn(*args, **kwargs):
    """
    Return a gaussian complex random number
    """
    return M_SQRT1_2 * (np.random.normal(*args, **kwargs)
                        + 1j * np.random.normal(*args, **kwargs))

class TestNPD:
    """
    Test data to compare against.
    """
    __test__ = False

    def __init__(self,
                 ptype=data.PType.ANY,
                 frequencies=0,
                 rows=0,
                 columns=0,
                 has_fz0=False):
        """
        Construct zero data with the given dimensions.
        """
        self.ptype = ptype
        self.frequencies = frequencies
        self.rows = rows
        self.columns = columns
        self.has_fz0 = has_fz0
        ports = max(rows, columns)
        self.ports = ports
        self.frequency_vector = np.zeros(shape=(frequencies,),
                                         dtype=np.double, order="C")
        self.data_array = np.zeros(shape=(frequencies, rows, columns),
                                   dtype=np.complex128, order="C")
        if has_fz0:
            self.fz0_array = np.zeros(shape=(frequencies, ports),
                                      dtype=np.complex128, order="C")
            self.fz0_array.fill(50.0)
        else:
            self.z0_vector = np.zeros(shape=(ports,),
                                      dtype=np.complex128, order="C")
            self.z0_vector.fill(50.0)

    def init(self):
        """
        Init to default values
        """
        self.frequency_vector.fill(0.0)
        self.data_array.fill(0.0+0.0j)
        if self.has_fz0:
            self.fz0_array.fill(50.0)
        else:
            self.z0_vector.fill(50.0)

    def randomize(self):
        """
        Fill with random values
        """
        for findex in range(self.frequencies):
            self.frequency_vector[findex] = np.random.normal()
            for row in range(self.rows):
                for column in range(self.columns):
                    self.data_array[findex, row, column] = crandn()
        if self.has_fz0:
            for findex in range(self.frequencies):
                for port in range(self.ports):
                    self.fz0_array[findex, port] = crandn()
        else:
            for port in range(self.ports):
                self.z0_vector[port] = crandn()

    def resize(self, ptype, frequencies, rows, columns, has_fz0):
        """
        Resize this object.
        """
        self.frequency_vector.resize((frequencies,))

        dv = np.zeros((frequencies, rows, columns), dtype=complex, order="C")
        for findex in range(min(self.frequencies, frequencies)):
            m = np.array(self.data_array[findex])
            m.resize(rows, columns)
            dv[findex, :, :] = m
        self.data_array = dv

        ports = max(rows, columns)
        if has_fz0:
            fz0v = np.zeros((frequencies, ports))
            if self.has_fz0:
                # preserve values
                for findex in range(min(self.ports, ports)):
                    v = np.array(self.fz0_array[findex])
                    v.resize((ports,))
                    fz0v[findex] = v
                    for port in range(self.ports, ports):
                        fz0v[findex, port] = 50.0
                self.fz0_array = fz0v

            else:
                # duplicate z0_vector into new matrix
                for findex in range(min(self.ports, ports)):
                    v = np.array(self.fz0_array)
                    v.resize((ports,))
                    fz0v[findex] = v
                    for port in range(self.ports, ports):
                        fz0v[findex, port] = 50.0
                self.fz0_array = fz0v

            self.z0_vector = None

        else:
            if self.has_fz0:
                # no preservation in this case
                z0v = np.zeros((ports,))
                z0v.fill(50)
                self.z0_vector = z0v

            else:
                # resize and init any new elements to 50.0
                self.z0_vector.resize((ports,))
                for port in range(self.ports, ports):
                    self.fz0_array[findex, port] = 50.0

            self.fz0_array = None

        self.ptype = ptype
        self.frequencies = frequencies
        self.rows = rows
        self.columns = columns
        self.has_fz0 = has_fz0
        self.ports = ports

    def fill(self, d):
        """
        Resize and fill in d to look like this object.
        """
        d.init(self.ptype,
               self.frequencies,
               self.rows,
               self.columns)
        d.frequency_vector = self.frequency_vector
        d.data_array = self.data_array
        if self.has_fz0:
            d.fz0_array = self.fz0_array
        else:
            d.z0_vector = self.z0_vector

    def isequal(self, d):
        """
        Test if d is equal to this object.
        """
        if d.ptype != self.ptype:
            print(f"ptype: {d.ptype} != {self.ptype}")
            return False
        if d.rows != self.rows:
            print(f"rows: {d.rows} != {self.rows}")
            return False
        if d.columns != self.columns:
            print(f"columns: {d.columns} != {self.columns}")
            return False
        if d.frequencies != self.frequencies:
            print(f"frequencies: {d.frequencies} != {self.frequencies}")
            return False
        if d.has_fz0 != self.has_fz0:
            print(f"has_fz0: {d.has_fz0} != {self.has_fz0}")
            return False
        if not np.array_equal(d.frequency_vector,
                              self.frequency_vector):
            print("frequency_vector does not match")
            print('d.frequency_vector', d.frequency_vector)
            print('self.frequency_vector', self.frequency_vector)
            return False
        if not np.array_equal(d.data_array, self.data_array):
            print("data_array does not match")
            print('d.data_array:', np.asarray(d.data_array))
            print('self.data_array:', self.data_array)
            return False
        if self.has_fz0:
            if not np.array_equal(d.fz0_array, self.fz0_array):
                print("fz0_array does not match")
                return False
        else:
            if not np.array_equal(d.z0_vector, self.z0_vector):
                print("z0_vector does not match")
                return False
        return True

def make_lc_abcd(l, c, f_vector):
    """
    make_lc_abcd(l, c, f_vector)
        Return a vector of matrices containing ABCD parameters for
        an L-C divider at the given frequencies.

    Parameters
    ----------
    l: double
        The value of the series inductor

    c: double
        The value of the shunt capacitor

    f_vector: vector of double
        The frequency points to evaluate.

    Returns
    -------
    a: (n, 2, 2) complex array
        A vector of 2x2 matrices giving the ABCD parameters
        for each frequency.
    """
    n = len(f_vector)
    a = np.empty((n, 2, 2), dtype=np.complex128)
    for findex in range(n):
        s = 2 * np.pi * 1j * f_vector[findex]
        ls = l * s
        cs = c * s
        a1 = np.asarray([[1, ls], [0, 1]]) # series inductor
        a2 = np.asarray([[1, 0], [cs, 1]]) # shunt capacitor
        a[findex, ...] = np.matmul(a1, a2)
    return a

class TestModule(unittest.TestCase):
    def test_zero_dim(self):
        """
        Test that the default constructor produces all dimensions zero
        without frequency dependent reference impedances.
        """
        d = data.NPData()
        self.assertEqual(d.frequencies, 0)
        self.assertEqual(d.rows, 0)
        self.assertEqual(d.columns, 0)
        self.assertEqual(len(d.frequency_vector), 0)
        self.assertEqual(len(d.data_array), 0)
        self.assertEqual(len(d.z0_vector), 0)
        self.assertFalse(d.has_fz0)

    def test_zero_values(self):
        d = data.NPData(data.PType.ANY, 5, 2, 3)
        self.assertEqual(d.frequencies, 5)
        self.assertEqual(d.rows, 2)
        self.assertEqual(d.columns, 3)
        self.assertEqual(len(d.frequency_vector), 5)
        self.assertEqual(len(d.data_array), 5)
        self.assertEqual(np.asarray(d.data_array).shape, (5, 2, 3))
        self.assertEqual(len(d.z0_vector), 3)
        self.assertFalse(d.has_fz0)
        t = TestNPD(data.PType.ANY, 5, 2, 3)
        self.assertTrue(t.isequal(d))

    def test_fill1(self):
        t = TestNPD(data.PType.ANY, 3, 1, 2, False)
        t.randomize()
        d = data.NPData()
        t.fill(d)
        self.assertTrue(t.isequal(d))

    def test_fill2(self):
        t = TestNPD(data.PType.ANY, 3, 1, 2, True)
        t.randomize()
        d = data.NPData()
        t.fill(d)
        self.assertTrue(t.isequal(d))

    def test_init(self):
        td1 = TestNPD(data.PType.ANY, 5, 3, 4, True)
        td1.randomize()
        d = data.NPData()
        td1.fill(d)
        d.init(data.PType.S, 4, 3, 3)
        td2 = TestNPD(data.PType.S, 4, 3, 3)
        self.assertTrue(td2.isequal(d))

    def test_get_frequency_vector(self):
        t = TestNPD(data.PType.ANY, 3, 5, 2, False)
        t.randomize()
        d = data.NPData()
        t.fill(d)
        self.assertEqual(len(d.frequency_vector), 3)
        self.assertTrue(np.array_equal(d.frequency_vector, t.frequency_vector))
        self.assertEqual(d.frequency_vector[0], t.frequency_vector[0])
        self.assertEqual(d.frequency_vector[1], t.frequency_vector[1])
        self.assertEqual(d.frequency_vector[2], t.frequency_vector[2])
        self.assertEqual(d.frequency_vector[-1], t.frequency_vector[-1])
        self.assertEqual(d.frequency_vector[-2], t.frequency_vector[-2])
        self.assertEqual(d.frequency_vector[-3], t.frequency_vector[-3])
        self.assertTrue(np.array_equal(d.frequency_vector[...],
                                       t.frequency_vector[...]))
        self.assertTrue(np.array_equal(d.frequency_vector[:],
                                       t.frequency_vector[:]))
        self.assertTrue(np.array_equal(d.frequency_vector[::-1],
                                       t.frequency_vector[::-1]))
        self.assertTrue(np.array_equal(d.frequency_vector[0:2],
                                       t.frequency_vector[0:2]))
        self.assertTrue(np.array_equal(d.frequency_vector[1:3],
                                       t.frequency_vector[1:3]))
        # Test forward iterator
        for i, f in enumerate(d.frequency_vector):
            self.assertEqual(f, t.frequency_vector[i])
            i += 1

        # Test reverse iterator
        i = 2
        for f in reversed(d.frequency_vector):
            self.assertEqual(f, t.frequency_vector[i])
            i -= 1

        with self.assertRaises(IndexError):
            _ = d.frequency_vector[3]
        with self.assertRaises(IndexError):
            _ = d.frequency_vector[-4]

    def test_set_frequency_vector(self):
        t = TestNPD(data.PType.ANY, 3, 5, 2, False)
        t.randomize()
        d = data.NPData()
        t.fill(d)

        d.frequency_vector[1] = 666.0
        t.frequency_vector[1] = 666.0
        self.assertTrue(t.isequal(d))

        d.frequency_vector[-3] = 125.0
        t.frequency_vector[-3] = 125.0
        self.assertTrue(t.isequal(d))

        d.frequency_vector[0:2] = [10.0, 20.0]
        t.frequency_vector[0:2] = [10.0, 20.0]
        self.assertTrue(t.isequal(d))

        d.frequency_vector[...] = [1, 2, 3]
        t.frequency_vector[...] = [1.0, 2.0, 3.0]
        self.assertTrue(t.isequal(d))

        with self.assertRaises(IndexError):
            d.frequency_vector[3] = 1.0
        with self.assertRaises(IndexError):
            d.frequency_vector[-4] = 1.0
        with self.assertRaises(ValueError):
            d.frequency_vector[...] = [1, 2, 3, 4]

        d.frequency_vector = [1, 2, 3, 4]
        self.assertEqual(len(d.frequency_vector), 4)
        self.assertEqual(d.frequencies, 4)
        self.assertTrue(np.array_equal(d.frequency_vector,
                                       [1.0, 2.0, 3.0, 4.0]))

    def test_get_data_vector(self):
        t = TestNPD(data.PType.ANY, 11, 4, 3, False)
        t.randomize()
        d = data.NPData()
        t.fill(d)

        self.assertEqual(len(d.data_array), 11)
        self.assertEqual(d.data_array[10, 3, 2],
                         t.data_array[10, 3, 2])
        self.assertEqual(d.data_array[-11, -4, -3],
                         t.data_array[-11, -4, -3])
        self.assertTrue(np.array_equal(d.data_array[5, 2, ::-2],
                                       t.data_array[5, 2, ::-2]))
        self.assertTrue(np.array_equal(d.data_array[5, 2, ...],
                                       t.data_array[5, 2, ...]))
        self.assertTrue(np.array_equal(d.data_array[5, ::-2, 2],
                                       t.data_array[5, ::-2, 2]))
        self.assertTrue(np.array_equal(d.data_array[5, ..., 1],
                                       t.data_array[5, ..., 1]))
        self.assertTrue(np.array_equal(d.data_array[5, ::-3, ::-2],
                                       t.data_array[5, ::-3, ::-2]))
        self.assertTrue(np.array_equal(d.data_array[5, ...],
                                       t.data_array[5, ...]))
        self.assertTrue(np.array_equal(d.data_array[::-3, 0, 0],
                                       t.data_array[::-3, 0, 0]))
        self.assertTrue(np.array_equal(d.data_array[..., 0, 0],
                                       t.data_array[..., 0, 0]))
        self.assertTrue(np.array_equal(d.data_array[0:5, 0, 0:2],
                                       t.data_array[0:5, 0, 0:2]))
        self.assertTrue(np.array_equal(d.data_array[0:5, 0, ...],
                                       t.data_array[0:5, 0, ...]))
        self.assertTrue(np.array_equal(d.data_array[..., 0, 0:2],
                                       t.data_array[..., 0, 0:2]))
        self.assertTrue(np.array_equal(d.data_array[::-1, 1:3, 0],
                                       t.data_array[::-1, 1:3, 0]))
        self.assertTrue(np.array_equal(d.data_array[..., 1:3, 0],
                                       t.data_array[..., 1:3, 0]))
        self.assertTrue(np.array_equal(d.data_array[..., 0],
                                       t.data_array[..., 0]))
        self.assertTrue(np.array_equal(d.data_array[::-1, ::-1, ::-1],
                                       t.data_array[::-1, ::-1, ::-1]))
        self.assertTrue(np.array_equal(d.data_array[...],
                                       t.data_array[...]))

        # Test forward iterator
        for i, x in enumerate(d.data_array):
            self.assertTrue(np.array_equal(x, t.data_array[i]))
            i += 1

        # Test reverse iterator
        i = 10
        for x in reversed(d.data_array):
            self.assertTrue(np.array_equal(x, t.data_array[i]))
            i -= 1

    def test_set_data_vector(self):
        t = TestNPD(data.PType.ANY, 11, 4, 3, False)
        t.randomize()
        d = data.NPData()
        t.fill(d)

        v = crandn()
        d.data_array[10, 3, 2] = v
        t.data_array[10, 3, 2] = v
        self.assertTrue(t.isequal(d))

        v = crandn()
        d.data_array[-11, -4, -3] = v
        t.data_array[-11, -4, -3] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(2,))
        d.data_array[5, 2, ::-2] = v
        t.data_array[5, 2, ::-2] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(3,))
        d.data_array[5, 2, ...] = v
        t.data_array[5, 2, ...] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(2,))
        d.data_array[5, ::-2, 2] = v
        t.data_array[5, ::-2, 2] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(4,))
        d.data_array[5, ..., 1] = v
        t.data_array[5, ..., 1] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(1, 2))
        d.data_array[5, ::-3, ::-2] = v
        t.data_array[5, ::-3, ::-2] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(4, 3))
        d.data_array[5, ...] = v
        t.data_array[5, ...] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(4,))
        d.data_array[::-3, 0, 0] = v
        t.data_array[::-3, 0, 0] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(11,))
        d.data_array[..., 0, 0] = v
        t.data_array[..., 0, 0] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(5, 2))
        d.data_array[0:5, 0, 0:2] = v
        t.data_array[0:5, 0, 0:2] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(5, 3))
        d.data_array[0:5, 0, ...] = v
        t.data_array[0:5, 0, ...] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(11, 2))
        d.data_array[..., 0, 0:2] = v
        t.data_array[..., 0, 0:2] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(11, 2))
        d.data_array[::-1, 1:3, 0] = v
        t.data_array[::-1, 1:3, 0] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(11, 2))
        d.data_array[..., 1:3, 0] = v
        t.data_array[..., 1:3, 0] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(11, 4))
        d.data_array[..., 0] = v
        t.data_array[..., 0] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(11, 4, 3))
        d.data_array[::-1, ::-1, ::-1] = v
        t.data_array[::-1, ::-1, ::-1] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(11, 4, 3))
        d.data_array[...] = v
        t.data_array[...] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(5, 3, 3))
        d.data_array = v
        self.assertEqual(d.frequencies, 5)
        self.assertEqual(d.rows, 3)
        self.assertEqual(d.columns, 3)
        t.resize(data.PType.ANY, 5, 3, 3, False)
        t.data_array = v
        self.assertTrue(t.isequal(d))

    def test_get_z0_vector(self):
        t = TestNPD(data.PType.ANY, 2, 3, 5, False)
        t.randomize()
        d = data.NPData()
        t.fill(d)

        self.assertEqual(len(d.z0_vector), 5)
        self.assertTrue(np.array_equal(d.z0_vector, t.z0_vector))
        for i in range(5):
            self.assertEqual(d.z0_vector[i], t.z0_vector[i])
        for i in range(-5, 0):
            self.assertEqual(d.z0_vector[i], t.z0_vector[i])
        self.assertTrue(np.array_equal(d.z0_vector[...],
                                       t.z0_vector[...]))
        self.assertTrue(np.array_equal(d.z0_vector[:],
                                       t.z0_vector[:]))
        self.assertTrue(np.array_equal(d.z0_vector[::-1],
                                       t.z0_vector[::-1]))
        self.assertTrue(np.array_equal(d.z0_vector[1:4],
                                       t.z0_vector[1:4]))
        # Test forward iterator
        for i, f in enumerate(d.z0_vector):
            self.assertEqual(f, t.z0_vector[i])
            i += 1

        # Test reverse iterator
        i = 4
        for f in reversed(d.z0_vector):
            self.assertEqual(f, t.z0_vector[i])
            i -= 1

        with self.assertRaises(IndexError):
            _ = d.z0_vector[5]
        with self.assertRaises(IndexError):
            _ = d.z0_vector[-6]

    def test_set_z0_vector(self):
        t = TestNPD(data.PType.ANY, 2, 5, 3, False)
        t.randomize()
        d = data.NPData()
        t.fill(d)

        v = crandn()
        d.z0_vector[0] = v
        t.z0_vector[0] = v
        self.assertTrue(t.isequal(d))

        v = crandn()
        d.z0_vector[4] = v
        t.z0_vector[4] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(2,))
        d.z0_vector[1:3] = v
        t.z0_vector[1:3] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(5,))
        d.z0_vector[...] = v
        t.z0_vector[...] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(6,))
        with self.assertRaises(ValueError):
            d.z0_vector[...] = v

        # Test that resize with conflicting data present raises
        v = crandn(size=(3,))
        with self.assertRaises(ValueError):
            d.z0_vector = v

        # Test that resize from (0 x 0) doesn't raise on same v
        d.resize(data.PType.ANY, 2, 0, 0)
        self.assertEqual(d.rows, 0)
        self.assertEqual(d.columns, 0)
        d.z0_vector = v
        self.assertEqual(d.rows, 3)
        self.assertEqual(d.columns, 3)

    def test_get_fz0_vector(self):
        t = TestNPD(data.PType.ANY, 5, 3, 4, True)
        t.randomize()
        d = data.NPData()
        t.fill(d)
        self.assertTrue(t.isequal(d))

        self.assertEqual(d.fz0_array[0, 0],
                         t.fz0_array[0, 0])
        self.assertEqual(d.fz0_array[4, 3],
                         t.fz0_array[4, 3])

        self.assertTrue(np.array_equal(d.fz0_array[4, ::-1],
                                       t.fz0_array[4, ::-1]))

        self.assertTrue(np.array_equal(d.fz0_array[0, ...],
                                       t.fz0_array[0, ...]))

        self.assertTrue(np.array_equal(d.fz0_array[::-2, 2],
                                       t.fz0_array[::-2, 2]))

        self.assertTrue(np.array_equal(d.fz0_array[..., 1],
                                       t.fz0_array[..., 1]))

        self.assertTrue(np.array_equal(d.fz0_array[::-1, 1:],
                                       t.fz0_array[::-1, 1:]))

        self.assertTrue(np.array_equal(d.fz0_array[...],
                                       t.fz0_array[...]))

        # Test forward iterator
        for i, x in enumerate(d.fz0_array):
            self.assertTrue(np.array_equal(x, t.fz0_array[i, :]))
            i += 1

        # Test reverse iterator
        i = 4
        for x in reversed(d.fz0_array):
            self.assertTrue(np.array_equal(x, t.fz0_array[i, :]))
            i -= 1

        with self.assertRaises(IndexError):
            _ = d.fz0_array[5]
        with self.assertRaises(IndexError):
            _ = d.fz0_array[-6]

    def test_set_fz0_vector(self):
        t = TestNPD(data.PType.ANY, 5, 4, 3, True)
        t.randomize()
        d = data.NPData()
        t.fill(d)
        self.assertTrue(t.isequal(d))

        v = crandn()
        d.fz0_array[0, 0] = v
        t.fz0_array[0, 0] = v
        self.assertTrue(t.isequal(d))

        v = crandn()
        d.fz0_array[4, 3] = v
        t.fz0_array[4, 3] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(2,))
        d.fz0_array[1, ::-2] = v
        t.fz0_array[1, ::-2] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(3,))
        d.fz0_array[::-2, 1] = v
        t.fz0_array[::-2, 1] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(5, 4))
        d.fz0_array[...] = v
        t.fz0_array[...] = v
        self.assertTrue(t.isequal(d))

        v = crandn(size=(6,))
        with self.assertRaises(ValueError):
            d.fz0_array[...] = v

        # Test that resize with conflicting data present raises
        v = crandn(size=(4, 3))
        with self.assertRaises(ValueError):
            d.fz0_array = v

        # Test that resize from (0 x 0) doesn't raise on same v
        d.resize(data.PType.ANY, 0, 0, 0)
        self.assertEqual(d.rows, 0)
        self.assertEqual(d.columns, 0)
        d.fz0_array = v
        self.assertEqual(d.frequencies, 4)
        self.assertEqual(d.rows, 3)
        self.assertEqual(d.columns, 3)

    def test_convert(self):
        #
        # 2nd order L-C low-pass filter with cut-off frequency of 1 MHz
        #
        l = 11.25395395196383e-06
        c = 2.250790790392765e-09
        frequency_vector = np.logspace(3, 9, 7)

        #
        # Expected S parameters to 16 places generated using Mathematica
        #
        s_expected = np.asarray([
            [[-1.249997343750000e-07 +3.535536999523635e-04j,
              +9.999993749998281e-01 -1.060660039197100e-03j],
             [+9.999993749998281e-01 -1.060660039197100e-03j,
              +8.749996406248000e-07 +3.535526392923243e-04j]],
            [[-1.249734375196000e-05 +3.535843252443067e-03j,
              +9.999374982814278e-01 -1.060646911041820e-02j],
             [+9.999374982814278e-01 -1.060646911041820e-02j,
              +8.749640607618300e-05 +3.534782605532025e-03j]],
            [[-1.223440113854840e-03 +3.566342158739672e-02j,
              +9.937329904371927e-01 -1.059309552100893e-01j],
             [+9.937329904371927e-01 -1.059309552100893e-01j,
              +8.713889790517090e-03 +3.460411203529583e-02j]],
            [[+9.090909090909090e-02 +5.142594772265801e-01j,
              +3.636363636363636e-01 -7.713892158398700e-01j],
             [+3.636363636363636e-01 -7.713892158398700e-01j,
              +4.545454545454545e-01 -2.571297386132900e-01j]],
            [[+9.896558583648299e-01 +1.420684004373211e-01j,
              -1.949472846628208e-02 -4.219853478336270e-03j],
             [-1.949472846628208e-02 -4.219853478336270e-03j,
              -9.598169882633778e-01 -2.799169473963059e-01j]],
            [[+9.998999650057511e-01 +1.414284212947320e-02j,
              -1.999499945022751e-04 -4.242428396002300e-06j],
             [-1.999499945022751e-04 -4.242428396002300e-06j,
              -9.995999800169999e-01 -2.828144183055030e-02j]],
            [[+9.999989999965000e-01 +1.414214269473900e-03j,
              -1.999994999994500e-06 -4.242638565783000e-09j],
             [-1.999994999994500e-06 -4.242638565783000e-09j,
              -9.999959999980000e-01 -2.828424296309200e-03j]]])

        #
        # Make ABCD parameters for the filter and store in
        # libvna.data.NPData object
        #
        A = data.NPData(data.PType.A, len(frequency_vector), 2, 2)
        A.frequency_vector = frequency_vector
        A.data_array = make_lc_abcd(l, c, frequency_vector)

        #
        # Convert to S parameters and check
        #
        S = A.convert(data.PType.S)
        self.assertTrue(np.allclose(S.data_array, s_expected))

    def test_renormalize(self):
        #
        # We have S parameters for an input matching network consisting of
        # a 31.8pF capacitor in parallel with the primary of a transformer
        # with L1=79.6nH, L2=2.91nH, M=9.42nH, and secondary to match
        # the input of an MRF314 RF transistor.  The network is tuned to
        # have a perfect match at 100MHz.
        #
        S = data.NPData(data.PType.S, 6, 2, 2)
        S.frequency_vector = [30e+6, 50e+6, 70e+6, 100e+6, 150e+6, 200e+6]
        S.data_array = [
            [[-0.8017510606625184 + 0.5930606677473840j,
               0.02366743348740367 + 0.07010088223177492j],
             [+0.02366743348740367 + 0.07010088223177492j,
              -0.9971567737488537 + 0.0142859185642236j]],
            [[-0.3814067753351794 + 0.9151217531805256j,
               0.07373248742615983 + 0.10791000362418888j],
             [+0.07373248742615983 + 0.10791000362418888j,
              -0.9911622349606514 + 0.0227238949915087j]],
            [[+0.3007047211854493 + 0.9346986201701367j,
              +0.1546680445869756 + 0.1095123558894842j],
             [+0.1546680445869756 + 0.1095123558894842j,
              -0.9815116323597084 + 0.0268282581200657j]],
            [[+0.9723892731275819 + 0.0002722951368089j,
              +0.2333532329237111 - 0.0023013139341589j],
             [+0.2333532329237111 - 0.0023013139341589j,
              -0.9721947766012899 + 0.0194496526292040j]],
            [[+0.1775023682932040 - 0.9674633141308717j,
              +0.1375759274824572 - 0.1165354419846409j],
             [+0.1375759274824572 - 0.1165354419846409j,
              -0.9834823431001016 + 0.0159635756997033j]],
            [[-0.3812098616018271 - 0.9152042788264460j,
              +0.07102205911898073 - 0.10970887180980763j],
             [+0.07102205911898073 - 0.10970887180980763j,
              -0.9910597513127330 + 0.0268402662554413j]]]
        #
        # New reference impedance table.  First column has the reference
        # impedances at the input that remain at 50 ohms.  Second column
        # has the input impedances to the MRF314 transistor from the
        # data sheet.
        #
        new_z0 = [
            [50.0, 2.4-3.4j],
            [50.0, 1.6-2.6j],
            [50.0, 0.8-0.8j],
            [50.0, 0.7-0.5j],
            [50.0, 0.9-0.9j],
            [50.0, 1.3-1.2j]]

        #
        # Expected S parameters after renormalization.
        #
        s_expected = [
            [[-0.7700183510870420 + 0.6034465948626218j,
              -0.1119063196176172 + 0.1743586010712961j],
             [-0.1119063196176172 + 0.1743586010712961j,
              +0.2279004088509418 - 0.9513871228409989j]],
            [[-0.2351367496948131 + 0.8709384171839692j,
              -0.1077514364756572 + 0.4178116924566300j],
             [-0.1077514364756572 + 0.4178116924566300j,
              +0.2153970026752552 - 0.8760291940469079j]],
            [[+0.1492544520597341 + 0.2448640692738425j,
              +0.7187665931789708 + 0.6333555720516708j],
             [+0.7187665931789708 + 0.6333555720516708j,
              -0.2616990227637168 - 0.1172558116491983j]],
            [[+0, # perfect match here
              +1.000000000000000],
             [+1.000000000000000,
              +0]],
            [[-0.1231028667790894 - 0.4789207419898820j,
              +0.8209413848834480 - 0.2855448295246776j],
             [+0.8209413848834480 - 0.2855448295246776j,
              -0.2006745362563238 - 0.4519392917369309j]],
            [[-0.3678539359473212 - 0.6548937443879078j,
              +0.5253660971586890 - 0.3997350739250895j],
             [+0.5253660971586890 - 0.3997350739250895j,
              -0.5330739299610895 - 0.5291828793774319j]]]

        #
        # Renormalize and check.
        #
        S.convert(data.PType.S, new_z0=new_z0, inplace=True)
        self.assertTrue(np.allclose(S.data_array, s_expected))
        self.assertTrue(np.allclose(S.fz0_array, new_z0))


    def test_load_save(self):
        n = 10
        t = TestNPD(data.PType.S, n, 2, 2, False)
        t.frequency_vector = np.logspace(6, 10, n)
        t.z0_vector[:] = [50, 75]
        t.data_array = crandn(size=(n, 2, 2))
        d = data.NPData()
        t.fill(d)
        d.format = "Sdb"
        d.fprecision = 8
        d.dprecision = 8
        d.save("test.ts")
        d = None

        d = data.NPData()
        d.load("test.ts")
        self.assertEqual(d.ptype, data.PType.S)
        self.assertTrue(np.allclose(d.frequency_vector, t.frequency_vector))
        self.assertTrue(np.allclose(d.z0_vector, t.z0_vector))
        self.assertTrue(np.allclose(d.data_array, t.data_array))

        os.remove("test.ts")


if __name__ == '__main__':
    unittest.main()
