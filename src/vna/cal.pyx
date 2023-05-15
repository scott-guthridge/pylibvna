#cython: language_level=3
#
# Python Bindings for Vector Network Analyzer Library
# Copyright © 2023 D Scott Guthridge <scott_guthridge@rompromity.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

from cpython.exc cimport PyErr_SetFromErrno
from cpython.pycapsule cimport PyCapsule_GetPointer
import errno
from libc.stdio cimport FILE, fdopen, fclose
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy
import numpy as np
cimport numpy as npc
from threading import local
import warnings
import vna.data
from vna.data import Data

npc.import_array()   # otherwise, PyArray_SimpleNewFromData segfaults

cpdef enum CalType:
    # """
    # Calibration error term type
    # """
    T8   = VNACAL_T8
    U8   = VNACAL_U8
    TE10 = VNACAL_TE10
    UE10 = VNACAL_UE10
    T16  = VNACAL_T16
    U16  = VNACAL_U16
    UE14 = VNACAL_UE14
    E12  = VNACAL_E12


cdef void _error_fn(const char *message, void *error_arg,
                    vnaerr_category_t category) noexcept:
    # """
    # C callback function for vnaerr
    # """
    self = <CalSet>error_arg
    if self._thread_local._vna_cal_exception is not None:
        return
    umessage = (<bytes>message).decode("utf8")
    if   category == VNAERR_SYSTEM:
        if errno.errorcode == errno.ENOMEM:
            self._thread_local._vna_cal_exception = MemoryError(umessage)
        else:
            self._thread_local._vna_cal_exception = OSError(errno.errorcode,
                                                            umessage)
    elif category == VNAERR_USAGE:
        self._thread_local._vna_cal_exception = ValueError(umessage)
    elif category == VNAERR_VERSION:
        self._thread_local._vna_cal_exception = ValueError(umessage)
    elif category == VNAERR_SYNTAX:
        self._thread_local._vna_cal_exception = SyntaxError(umessage, None)
    elif category == VNAERR_WARNING:
        if self._thread_local._vna_cal_warning is None:
            self._thread_local._vna_cal_warning = umessage
    elif category == VNAERR_MATH:
        self._thread_local._vna_cal_exception = ArithmeticError(umessage)
    elif category == VNAERR_INTERNAL:
        self._thread_local._vna_cal_exception = AssertionError(umessage)
    else:
        self._thread_local._vna_cal_exception = Exception(umessage)


cdef object _property_to_py(vnaproperty_t *root):
    # """
    # Convert a vnaproperty_t tree to tree of python objects.
    # """
    cdef int i
    cdef int count
    cdef int ptype = vnaproperty_type(root, ".")
    cdef const char *value
    cdef const char **keys = NULL
    cdef char *quoted = NULL
    cdef vnaproperty_t *subtree

    if root == NULL:
        return None

    if ptype == ord("s"):
        value = vnaproperty_get(root, ".")
        return value.decode("UTF-8")

    if ptype == ord("m"):
        result = {}
        keys = vnaproperty_keys(root, "{}")
        assert(keys != NULL)
        try:
            i = 0
            while keys[i] != NULL:
                quoted = vnaproperty_quote_key(keys[i])
                if quoted == NULL:
                    raise MemoryError()
                subtree = vnaproperty_get_subtree(root, "%s", quoted)
                free(quoted)
                quoted = NULL
                strkey = keys[i].decode("UTF-8")
                result[strkey] = _property_to_py(subtree)
                i += 1
        finally:
            free(quoted)
            free(keys)
        return result

    if ptype == ord("l"):
        result = []
        count = vnaproperty_count(root, "[]")
        assert(count >= 0)
        for i in range(count):
            subtree = vnaproperty_get_subtree(root, "[%d]", i)
            result.append(_property_to_py(subtree))
        return result

    raise AssertionError(f"_property_to_py: unknown type {ptype}")


cdef void _py_to_property(object root, vnaproperty_t **rootptr):
    # """
    # Convert a tree of python objects to vnaproperty_t tree
    # """
    cdef int i
    cdef int rc
    cdef char *quoted = NULL
    cdef vnaproperty_t **subtreeptr
    cdef const unsigned char[:] cstring

    if root is None:
        rootptr[0] = NULL
        return

    if isinstance(root, str):
        if isinstance(root, unicode):
            root = (<unicode>root).encode("UTF-8")
        cstring = root
        rc = vnaproperty_set(rootptr, ".=%s", <const char *>&cstring[0])
        if rc == -1:
            raise OSError(errno.errorcode)
        return

    if isinstance(root, bytes):
        cstring = root
        rc = vnaproperty_set(rootptr, ".=%s", <const char *>&cstring[0])
        if rc == -1:
            raise OSError(errno.errorcode)
        return

    if isinstance(root, dict):
        if vnaproperty_set_subtree(rootptr, "{}") == NULL:
            raise MemoryError()
        for key, value in root.items():
            if isinstance(key, unicode):
                key  = (<unicode>key).encode("UTF-8")
            elif not isinstance(key, bytes):
                raise ValueError("dictionary key must be str or bytes")
            cstring = key
            quoted = vnaproperty_quote_key(<const char *>&cstring[0])
            if quoted == NULL:
                raise MemoryError()
            subtreeptr = vnaproperty_set_subtree(rootptr, "%s", quoted)
            if subtreeptr == NULL:
                free(quoted)
                quoted = NULL
                raise MemoryError()
            free(quoted)
            quoted = NULL
            _py_to_property(value, subtreeptr)
        return

    if isinstance(root, list):
        if vnaproperty_set_subtree(rootptr, "[]") == NULL:
            raise MemoryError()
        for i, value in enumerate(root):
            subtreeptr = vnaproperty_set_subtree(rootptr, "[%d]", i)
            if subtreeptr == NULL:
                raise MemoryError()
            _py_to_property(value, subtreeptr)
        return

    raise AssertionError(f"_py_to_property: invalid type {type(root)}")


cdef class Parameter:
    """
    A known or unknown parameter of a calibration standard
    """
    cdef CalSet calset
    cdef int pindex

    def __cinit__(self):
        self.calset = None
        self.pindex = -1

    def __init__(self):
        raise TypeError("This class cannot be instantiated directly.")

    def __dealloc__(self):
        cdef CalSet calset = self.calset
        cdef pindex = self.pindex
        if calset is not None and pindex >= 3:
            vnacal_delete_parameter(calset.vcp, pindex)

    def get_value(self, frequency):
        cdef vnacal_t *vcp = self.calset.vcp
        cdef complex result
        result = vnacal_get_parameter_value(vcp, self.pindex, frequency)
        self.calset._handle_error(0)
        return result


cdef object _prepare_C_array(object array, object name, int frequencies,
        double complex ***clfppp, int *rows, int *columns):
    # """
    # Given a (rows x columns x frequencies) array-like object, return
    # a flattened rows x columns matrix of pointers to frequencies long
    # vectors of values expected by the C code.  Because the resulting
    # C array points into the python ndarray, we have to make sure the
    # ndarray remains referenced while the C array is in-use, thus array
    # is passed by reference.
    # """
    array = np.asarray(array, dtype=np.complex128, order="C")
    if array.ndim != 3:
        raise ValueError(f"{name} must be a (rows x columns x frequencies) "
                         f"array")
    cdef int r = array.shape[0]
    cdef int c = array.shape[1]
    if array.shape[2] != frequencies:
        raise ValueError(f"third dimension of {name} must be {frequencies}")
    cdef double complex[:, :, ::1] v = array
    cdef double complex **clfpp = <double complex **>malloc(
            r * c * sizeof(double complex *))
    if clfpp == NULL:
        raise MemoryError()
    cdef int i
    cdef int j
    cdef int k
    try:
        k = 0
        for i in range(r):
            for j in range(c):
                clfpp[k] = &v[i, j, 0]
                k += 1

        clfppp[0] = clfpp
        clfpp = NULL    # prevent free at finally
        rows[0] = r
        columns[0] = c
        return array

    finally:
        free(<void *>clfpp)


cdef class Solver:
    cdef CalSet calset
    cdef int frequencies
    cdef vnacal_new_t *vnp
    cdef double _pvalue_limit
    cdef double _et_tolerance
    cdef double _p_tolerance
    cdef int _iteration_limit

    def __cinit__(self):
        calset = NULL
        frequencies = 0
        vnp = NULL
        _pvalue_limit = 0.001
        _et_tolerance = 1.0e-6
        _p_tolerance = 1.0e-6
        _iteration_limit = 30

    def __init__(self):
        raise TypeError("This class cannot be instantiated directly.")

    def __dealloc__(self):
        # Note that the reference on calset is released after this
        # function returns.
        vnacal_new_free(self.vnp)
        self.vnp = NULL

    def add_single_reflect(self, a, b, s11, int port):
        """
        Add the measurement of a single reflect standard with parameter
        s11 on the given VNA port.  VNA port numbers start at 1.
        """
        cdef int a_rows = 0
        cdef int a_columns = 0
        cdef int b_rows
        cdef int b_columns
        cdef double complex **a_clfpp = NULL
        cdef double complex **b_clfpp = NULL
        cdef Parameter c_s11
        cdef int rc

        try:
            #
            # Prepare a and b arrays.
            #
            if a is not None:
                a = _prepare_C_array(a, "a", self.frequencies,
                                     &a_clfpp, &a_rows, &a_columns)
            b = _prepare_C_array(b, "b", self.frequencies,
                                 &b_clfpp, &b_rows, &b_columns)

            #
            # Prepare S parameters
            #
            c_s11 = self.calset._generalize_parameter(s11)

            #
            # Call add function
            #
            if a_clfpp != NULL:
                rc = vnacal_new_add_single_reflect(self.vnp,
                                                   a_clfpp, a_rows, a_columns,
                                                   b_clfpp, b_rows, b_columns,
                                                   c_s11.pindex, port)
            else:
                rc = vnacal_new_add_single_reflect_m(self.vnp,
                                                     b_clfpp, b_rows, b_columns,
                                                     c_s11.pindex, port)
            self.calset._handle_error(rc)

            return

        finally:
            free(<void *>b_clfpp)
            free(<void *>a_clfpp)

    def add_double_reflect(self, a, b, s11, s22, int port1, int port2):
        """
        Add the measurement of a double reflect standard with parameters
        s11 and s22 on the given VNA ports.  VNA port numbers start at 1.
        """
        cdef int a_rows = 0
        cdef int a_columns = 0
        cdef int b_rows
        cdef int b_columns
        cdef double complex **a_clfpp = NULL
        cdef double complex **b_clfpp = NULL
        cdef Parameter c_s11
        cdef Parameter c_s22
        cdef int rc

        try:
            #
            # Prepare a and b arrays.
            #
            if a is not None:
                a = _prepare_C_array(a, "a", self.frequencies,
                                     &a_clfpp, &a_rows, &a_columns)
            b = _prepare_C_array(b, "b", self.frequencies,
                                 &b_clfpp, &b_rows, &b_columns)

            #
            # Prepare S parameters
            #
            c_s11 = self.calset._generalize_parameter(s11)
            c_s22 = self.calset._generalize_parameter(s22)

            #
            # Call add function
            #
            if a_clfpp != NULL:
                rc = vnacal_new_add_double_reflect(self.vnp,
                                                   a_clfpp, a_rows, a_columns,
                                                   b_clfpp, b_rows, b_columns,
                                                   c_s11.pindex, c_s22.pindex,
                                                   port1, port2)
            else:
                rc = vnacal_new_add_double_reflect_m(self.vnp,
                                                     b_clfpp, b_rows, b_columns,
                                                     c_s11.pindex, c_s22.pindex,
                                                     port1, port2)
            self.calset._handle_error(rc)

            return

        finally:
            free(<void *>b_clfpp)
            free(<void *>a_clfpp)

    def add_through(self, a, b, int port1, int port2):
        """
        Add the measurement of a perfect through standard between
        port1 and port2.  VNA port numbers start at 1.
        """
        cdef int a_rows = 0
        cdef int a_columns = 0
        cdef int b_rows
        cdef int b_columns
        cdef double complex **a_clfpp = NULL
        cdef double complex **b_clfpp = NULL
        cdef int rc

        try:
            #
            # Prepare a and b arrays.
            #
            if a is not None:
                a = _prepare_C_array(a, "a", self.frequencies,
                                     &a_clfpp, &a_rows, &a_columns)
            b = _prepare_C_array(b, "b", self.frequencies,
                                 &b_clfpp, &b_rows, &b_columns)

            #
            # Call add function
            #
            if a_clfpp != NULL:
                rc = vnacal_new_add_through(self.vnp,
                                            a_clfpp, a_rows, a_columns,
                                            b_clfpp, b_rows, b_columns,
                                            port1, port2)
            else:
                rc = vnacal_new_add_through_m(self.vnp,
                                              b_clfpp, b_rows, b_columns,
                                              port1, port2)
            self.calset._handle_error(rc)

            return

        finally:
            free(<void *>b_clfpp)
            free(<void *>a_clfpp)

    def add_line(self, a, b, s, int port1, int port2):
        """
        Add the measurement of an arbitrary two-port standard with S
        parameter matrix s on the given VNA ports.  VNA port numbers
        start at 1.
        """
        cdef int a_rows = 0
        cdef int a_columns = 0
        cdef int b_rows
        cdef int b_columns
        cdef double complex **a_clfpp = NULL
        cdef double complex **b_clfpp = NULL
        cdef int i
        cdef int j
        cdef Parameter c_parameter
        cdef int si[2][2]
        cdef int rc

        try:
            #
            # Prepare a and b arrays.
            #
            if a is not None:
                a = _prepare_C_array(a, "a", self.frequencies,
                                     &a_clfpp, &a_rows, &a_columns)
            b = _prepare_C_array(b, "b", self.frequencies,
                                 &b_clfpp, &b_rows, &b_columns)

            #
            # Prepare S parameters
            #
            s = np.asarray(s, dtype=object)
            if s.ndim < 2 or s.shape[0] != 2 or s.shape[1] != 2:
                raise ValueError("s must be a 2x2 matrix of parameters")
            for i in range(2):
                for j in range(2):
                    c_parameter = self.calset._generalize_parameter(s[i, j])
                    s[i, j] = c_parameter
                    si[i][j] = c_parameter.pindex

            #
            # Call add function
            #
            if a_clfpp != NULL:
                rc = vnacal_new_add_line(self.vnp,
                                         a_clfpp, a_rows, a_columns,
                                         b_clfpp, b_rows, b_columns,
                                         &si[0][0], port1, port2)
            else:
                rc = vnacal_new_add_line_m(self.vnp,
                                           b_clfpp, b_rows, b_columns,
                                           &si[0][0], port1, port2)
            self.calset._handle_error(rc)

            return

        finally:
            free(<void *>b_clfpp)
            free(<void *>a_clfpp)

    def add_mapped_matrix(self, a, b, s, port_map):
        """
        Add the measurement of an arbitrary n-port standard with S
        parameter matrix s, and map of ports of the standard to ports
        of the VNA in port_map.  VNA port numbers start at 1.
        """
        cdef int a_rows = 0
        cdef int a_columns = 0
        cdef int b_rows
        cdef int b_columns
        cdef int s_rows
        cdef int s_columns
        cdef int s_ports
        cdef double complex **a_clfpp = NULL
        cdef double complex **b_clfpp = NULL
        cdef int i
        cdef int j
        cdef int k
        cdef Parameter c_parameter
        cdef int *sip = NULL
        cdef int [::1] iv
        cdef int rc

        try:
            #
            # Prepare a and b arrays.
            #
            if a is not None:
                a = _prepare_C_array(a, "a", self.frequencies,
                                     &a_clfpp, &a_rows, &a_columns)
            b = _prepare_C_array(b, "b", self.frequencies,
                                 &b_clfpp, &b_rows, &b_columns)

            #
            # Prepare S parameters
            #
            s = np.asarray(s, dtype=object)
            if s.ndim < 2:
                raise ValueError("s must be a matrix of parameters")
            s_rows = s.shape[0]
            s_columns = s.shape[1]
            s_ports = s_rows if s_rows >= s_columns else s_columns
            sip = <int *>malloc(s.rows * s.columns * sizeof(int))
            if sip == NULL:
                raise MemoryError()
            k = 0
            for i in range(s_rows):
                for j in range(s_columns):
                    c_parameter = self.calset._generalize_parameter(s[i, j])
                    s[i, j] = c_parameter
                    sip[k] = c_parameter.pindex
                    k += 1

            #
            # Prepare port map.
            #
            port_map = np.asarray(port_map, dtype=np.intc, order="C")
            if port_map.ndim != 1 or port_map.shape[0] < s_ports:
                raise ValueError(f"port_map must be a length {s_ports} vector")
            iv = port_map

            #
            # Call add function
            #
            if a_clfpp != NULL:
                rc = vnacal_new_add_mapped_matrix(self.vnp,
                                                  a_clfpp, a_rows, a_columns,
                                                  b_clfpp, b_rows, b_columns,
                                                  sip, s_rows, s_columns,
                                                  &iv[0])
            else:
                rc = vnacal_new_add_mapped_matrix_m(self.vnp,
                                                    b_clfpp, b_rows, b_columns,
                                                    sip, s_rows, s_columns,
                                                    &iv[0])
            self.calset._handle_error(rc)
            return

        finally:
            free(<void *>b_clfpp)
            free(<void *>a_clfpp)
            free(<void *>sip)

    def set_m_error(self, frequency_vector, noise_floor,
                    tracking_error = None):
        """
        Enable measurement error modeling.  Specifying measurement
        errors with this function can significanlty improve accuracy,
        especially for well overdetermined systems as the 16 term models
        typically are.

        Parameters:
            frequency_vector:
                vector of ascending frequencies spanning at least the
                calibration frequency range

            noise_floor:
                vector of noise floor root-power measurements at the
                VNA detectors at each frequency when no signal is applied

            tracking_error:
                optional vector describing of additional root-power
                noise source proportional to the measured signal

        Both noise sources are assumed Gaussian and independent.
        """
        cdef double [::1] fv
        cdef double [::1] nfv
        cdef double [::1] trv
        cdef double *trp = NULL
        cdef int rc

        frequency_vector = np.asarray(frequency_vector, dtype=np.double,
                                      order="C")
        if frequency_vector.ndim != 1:
            raise ValueError("frequency_vector must be a one-dimensional "
                             "array")
        cdef int n = len(frequency_vector)
        fv = frequency_vector
        if noise_floor is None:
            raise ValueError("noise_floor cannot be None")
        noise_floor = np.asarray(noise_floor, dtype=np.double, order="C")
        if len(noise_floor) != n:
            raise ValueError("length of noise_floor vector must match "
                             "that of frequency_vector")
        nfv = noise_floor
        if tracking_error is not None:
            tracking_error = np.asarray(tracking_error, dtype=np.double,
                                        order="C")
            if len(tracking_error) != n:
                raise ValueError("length of tracking_error vector must match "
                                 "that of frequency_vector")
            trv = tracking_error
            trp = &trv[0]

        rc = vnacal_new_set_m_error(self.vnp, &fv[0], n, &nfv[0], trp)
        self.calset._handle_error(rc)
        return

    @property
    def pvalue_limit(self):
        """
        p-value, below which to reject the null hypothesis that
        the measurement errors are distributed according to the
        values given in set_m_error.

        This parameter has no effect if set_m_error is not called
        to enable measurement error modeling.
        """
        return self._pvalue_limit

    @pvalue_limit.setter
    def pvalue_limit(self, double value):
        # no docstring for setter
        cdef int rc
        rc = vnacal_new_set_pvalue_limit(self.vnp, value)
        self.calset._handle_error(rc)
        self._pvalue_limit = value

    @property
    def et_tolerance(self):
        """
        Degree of change in the root-mean-squared of the error terms
        sufficiently small to stop iteration.  Default is 1.0e-6.
        """
        return self._et_tolerance

    @et_tolerance.setter
    def et_tolerance(self, double value):
        # no docstring for setter
        cdef int rc
        rc = vnacal_new_set_et_tolerance(self.vnp, value)
        self.calset._handle_error(rc)
        self._et_tolerance = value

    @property
    def p_tolerance(self):
        """
        Degree of change in the root-mean-squared of the unknown
        parameters sufficiently small to stop iteration.  This parameter
        has no effect if there are no unknown parameters in the S matrix.
        Default is 1.0e-6.
        """
        return self._et_tolerance

    @p_tolerance.setter
    def p_tolerance(self, double value):
        # no docstring for setter
        cdef int rc
        rc = vnacal_new_set_p_tolerance(self.vnp, value)
        self.calset._handle_error(rc)
        self._p_tolerance = value

    @property
    def iteration_limit(self):
        """
        For iterative solutions, this parameter controls the maximum
        number of iterations permitted to reach convergence.
        """
        return self._iteration_limit

    @iteration_limit.setter
    def iteration_limit(self, int value):
        # no docstring for setter
        cdef int rc
        rc = vnacal_new_set_iteration_limit(self.vnp, value)
        self.calset._handle_error(rc)
        self._iteration_limit = value

    def solve(self):
        """
        Solve for the error terms.
        """
        cdef int rc = vnacal_new_solve(self.vnp)
        self.calset._handle_error(rc)


cdef class Calibration:
    cdef CalSet calset
    cdef int ci

    @property
    def name(self) -> str:
        """
        Name of the calibration
        """
        cdef vnacal_t *vcp = self.calset.vcp
        cdef int ci = self.ci
        cdef const char *cp
        cp = vnacal_get_name(vcp, ci)
        assert(cp != NULL)
        return cp.decode("UTF-8")

    @property
    def ctype(self) -> CalType:
        """
        Type of the calibration
        """
        cdef vnacal_t *vcp = self.calset.vcp
        cdef int ci = self.ci
        cdef vnacal_type_t ctype = vnacal_get_type(vcp, ci)
        return <CalType>ctype

    @property
    def rows(self) -> int:
        """
        Number of rows in the calibration
        """
        cdef vnacal_t *vcp = self.calset.vcp
        cdef int ci = self.ci
        cdef int rows = vnacal_get_rows(vcp, ci)
        assert(rows != -1)
        return rows

    @property
    def columns(self) -> int:
        """
        Number of columns in the calibration
        """
        cdef vnacal_t *vcp = self.calset.vcp
        cdef int ci = self.ci
        cdef int columns = vnacal_get_columns(vcp, ci)
        assert(columns != -1)
        return columns

    @property
    def frequencies(self) -> int:
        """
        Number of frequencies in the calibration
        """
        cdef vnacal_t *vcp = self.calset.vcp
        cdef int ci = self.ci
        cdef int frequencies = vnacal_get_frequencies(vcp, ci)
        assert(frequencies != -1)
        return frequencies

    @property
    def frequency_vector(self):
        """
        Vector of calibration frequencies
        """
        cdef vnacal_t *vcp = self.calset.vcp
        cdef int ci = self.ci
        cdef int frequencies = vnacal_get_frequencies(vcp, ci)
        assert(frequencies != -1)
        result = np.empty((frequencies,), dtype=np.double, order="C")
        cdef double[::1] v = result
        cdef const double *lfp = vnacal_get_frequency_vector(vcp, ci)
        memcpy(&v[0], lfp, frequencies * sizeof(double))
        return result

    @property
    def z0(self):
        """
        Reference frequency all all ports in the calibration.
        """
        cdef vnacal_t *vcp = self.calset.vcp
        cdef int ci = self.ci
        cdef double complex z0 = vnacal_get_z0(vcp, ci)
        return z0

    def apply(self, f, a, b):
        """
        Apply the calibration correction to measured data

        Parameters:
            f: vector of frequencies at which the measurements were
               made, or None if measured at the calibration frequencies
            a: matrix of vectors (one entry per frequency) of incident
               voltages on each DUT port
            b: matrix of vectors (one entry per frequency) of reflected
               voltages from each DUT port

        Return:
            vna.data.Data object containing the corrected parameters
               
        """
        cdef vnacal_t *vcp = self.calset.vcp
        cdef int ci = self.ci
        cdef int frequencies
        cdef const double [::1] f_view
        cdef const double *frequency_vector
        cdef int a_rows = 0
        cdef int a_columns = 0
        cdef int b_rows
        cdef int b_columns
        cdef double complex **a_clfpp = NULL
        cdef double complex **b_clfpp = NULL
        cdef vnadata_t *vdp
        cdef int rc

        #
        # Process frequency vector.  If f is None, default to the
        # calibration frequencies.
        #
        if f is None:
            frequencies = vnacal_get_frequencies(vcp, ci)
            frequency_vector = vnacal_get_frequency_vector(vcp, ci)
            assert(frequency_vector != NULL)
        else:
            f_array = np.asarray(f, dtype=np.double, order="C")
            if f_array.ndim != 1:
                raise ValueError("f must be a one-dimensional array or None")
            frequencies = f_array.shape[0]
            f_view = f_array
            frequency_vector = &f_view[0]

        try:
            #
            # Process the input arrays.
            #
            if a is not None:
                a = _prepare_C_array(a, "a", frequencies,
                                     &a_clfpp, &a_rows, &a_columns)
            b = _prepare_C_array(b, "b", frequencies,
                                 &b_clfpp, &b_rows, &b_columns)

            #
            # Prepare result object
            #
            result = vna.data.Data()
            vdp = <vnadata_t *>PyCapsule_GetPointer(result._get_vdp(), NULL)

            #
            # Call apply
            #
            if a_clfpp != NULL:
                rc = vnacal_apply(vcp, ci,
                                  frequency_vector, frequencies,
                                  a_clfpp, a_rows, a_columns,
                                  b_clfpp, b_rows, b_columns,
                                  vdp)
            else:
                rc = vnacal_apply_m(vcp, ci,
                                  frequency_vector, frequencies,
                                  b_clfpp, b_rows, b_columns,
                                  vdp)
            self.calset._handle_error(rc)
            return result

        finally:
            free(<void *>b_clfpp)
            free(<void *>a_clfpp)

    @property
    def properties(self):
        """
        Reading from this property returns a tree of dictionaries, lists,
        scalars and None representing the properties of this calibration.
        Writing this this property replaces the property tree.
        """
        cdef vnacal_t *vcp = self.calset.vcp
        cdef int ci = self.ci
        cdef vnaproperty_t *root = vnacal_property_get_subtree(vcp, ci, ".")
        return _property_to_py(root)

    @properties.setter
    def properties(self, value):
        # no docstring for setter
        cdef vnacal_t *vcp = self.calset.vcp
        cdef int ci = self.ci
        cdef vnaproperty_t **rootptr
        cdef vnaproperty_t *new_root = NULL
        success = False
        try:
            _py_to_property(value, &new_root)
            rootptr = vnacal_property_set_subtree(vcp, ci, ".")
            rootptr[0] = new_root
            new_root = NULL
            success = True
        finally:
            if not success:
                vnaproperty_delete(&new_root, ".")


cdef class _CalHelper:
    # """
    # Internal class that provides __getitem__, __delitem__, etc on
    # calibrations.
    # """
    cdef CalSet calset

    def _find_ci_by_index_or_name(self, index):
        # """
        # Find the vnacal calibration index (ci) from the dense
        # index or name.
        # """
        cdef CalSet calset = self.calset
        cdef vnacal_t *vcp = calset.vcp
        cdef Calibration calibration
        if isinstance(index, int):
            return calset._index_to_ci[index]

        if isinstance(index, str):
            return calset.index(index)

        raise IndexError("index must have type int or str")

    def __getitem__(self, index):
        """
        Return a Calibration by name or position.
        """
        cdef CalSet calset = self.calset
        cdef int ci = self._find_ci_by_index_or_name(index)
        cdef Calibration calibration
        calibration = Calibration()
        calibration.calset = calset
        calibration.ci = ci
        return calibration

    def __len__(self):
        cdef CalSet calset = self.calset
        return len(calset._index_to_ci)

    def __delitem__(self, index):
        """
        Delete a calibration.
        """
        cdef CalSet calset = self.calset
        cdef vnacal_t *vcp = calset.vcp
        cdef int ci = self._find_ci_by_index_or_name(index)
        cdef int rc = vnacal_delete_calibration(vcp, ci)
        assert(rc == 0)
        del calset._index_to_ci[index]

    def __iter__(self):
        """
        Iterate over the Calibrations
        """
        cdef CalSet calset = self.calset
        cdef int i
        for i in range(len(calset._index_to_ci)):
            yield self[i]

    def  __reversed__(self):
        """
        Iterate reversed over the Calibrations
        """
        cdef CalSet calset = self.calset
        cdef int i
        for i in reversed(range(len(calset._index_to_ci))):
            yield self[i]

    def __string__(self):
        """
        Return a string representation of the Calibrations array
        """
        d = {}
        cdef int i
        for i, name in enumerate(self):
            d[i] = name
        return "Calibrations: " + str(d)


cdef class CalSet:
    """
    One or more vector network analyzer calibrations with a
    common save file
    """
    cdef vnacal_t *vcp
    cdef object _thread_local
    cdef object _index_to_ci

    def __cinit__(self, filename = None):
        cdef int ci
        cdef int ci_end
        cdef vnacal_type_t ctype
        cdef const unsigned char[:] cfilename
        self.vcp = NULL
        self._thread_local = local()
        self._thread_local._vna_cal_exception = None
        self._thread_local._vna_cal_warning = None
        if filename is not None:
            if isinstance(filename, unicode):
                filename = (<unicode>filename).encode("UTF-8")
            cfilename = filename
            self.vcp = vnacal_load(<const char *>&cfilename[0],
                                   <vnaerr_error_fn_t *>&_error_fn,
                                   <void *>self)
        else:
            self.vcp = vnacal_create(<vnaerr_error_fn_t *>&_error_fn,
                                     <void *>self)
        self._handle_error(0 if self.vcp != NULL else -1)

        #
        # Build a dense index of the saved calibrations.
        #
        self._index_to_ci = []
        ci_end = vnacal_get_calibration_end(self.vcp)
        for ci in range(ci_end):
            ctype = vnacal_get_type(self.vcp, ci)
            if ctype != VNACAL_NOTYPE:
                self._index_to_ci.append(ci)

    def __dealloc__(self):
        vnacal_free(self.vcp)

    @staticmethod
    def create():
        """
        Return an empty CalSet.
        """
        return CalSet();

    @staticmethod
    def load(filename):
        """
        Load a calibration file and return a CalSet object.
        """
        return CalSet(filename)

    def save(self, filename):
        """
        Save the CalSet to a file.
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode("UTF-8")
        cdef const unsigned char[:] cfilename = filename
        cdef int rc = vnacal_save(self.vcp, <const char *>&cfilename[0])
        self._handle_error(rc)

    def make_solver(self, CalType ctype, frequency_vector,
                    int rows, int columns, double complex z0 = 50.0):
        """
        Create and return a new error term solver.
        """
        frequency_vector = np.asarray(frequency_vector, dtype=np.double,
                                      order="C")
        if frequency_vector.ndim != 1:
            raise ValueError("frequency_vector must be 1 dimensional "
                             "array-like")
        cdef vnacal_t *vcp = self.vcp
        cdef vnacal_new_t *vnp = NULL
        vnp = vnacal_new_alloc(vcp, <vnacal_type_t>ctype, rows, columns,
                               len(frequency_vector))
        if vnp == NULL:
            self._handle_error(-1)
        cdef int rc
        cdef const double [::1] fv_view = frequency_vector
        rc = vnacal_new_set_frequency_vector(vnp, &fv_view[0])
        if rc == -1:
            vnacal_new_free(vnp)
            self._handle_error(-1)
        if z0 != 50.0:
            rc = vnacal_new_set_z0(vnp, z0)
            if rc == -1:
                vnacal_new_free(vnp)
                self._handle_error(-1)
        cdef Solver solver = Solver.__new__(Solver)
        solver.calset = self
        solver.frequencies = len(frequency_vector)
        solver.vnp = vnp
        return solver

    def add(self, Solver solver, name) -> int:
        """
        Add a new solved calibration to the CalSet.
        """
        cdef vnacal_t *vcp = self.vcp
        cdef vnacal_new_t *vnp = solver.vnp
        if isinstance(name, unicode):
            name = (<unicode>name).encode("UTF-8")
        cdef const unsigned char[:] c_name = name
        cdef int rc = vnacal_add_calibration(vcp, <const char *>&c_name[0],
                                             vnp)
        self._handle_error(rc)
        self._index_to_ci.append(rc)

    def index(self, name) -> int:
        """
        Return the calibration index with the given name.
        """
        cdef vnacal_t *vcp = self.vcp
        cdef int i
        cdef int ci
        cdef const char *cp
        for i, ci in enumerate(self._index_to_ci):
            cp = vnacal_get_name(vcp, ci)
            assert(cp != NULL)
            if cp.decode("UTF-8") == name:
                return i
        raise IndexError(f"calibration {name} not found")

    cdef Parameter _generalize_parameter(self, parameter):
        # """
        # If a number is found where a Parameter is expected, automatically
        # convert it into a scalar parameter.
        # """
        if isinstance(parameter, Parameter):
            return parameter
        cdef double complex value
        if isinstance(parameter, complex) or isinstance(parameter, float) \
                or isinstance(parameter, int):
            return self.make_scalar(parameter)
        raise ValueError("parameter must be class Parameter or a number")

    def make_scalar(self, double complex gamma):
        """
        Make a frequency-independent S parameter with value gamma.
        """
        cdef vnacal_t *vcp = self.vcp
        cdef int rc
        if gamma == 0.0:        # special-case these values
            rc = VNACAL_MATCH
        elif gamma == +1.0:
            rc = VNACAL_OPEN
        elif gamma == -1.0:
            rc = VNACAL_SHORT
        else:
            rc = vnacal_make_scalar_parameter(vcp, gamma)
            self._handle_error(rc)
        cdef Parameter parameter = Parameter.__new__(Parameter)
        parameter.calset = self
        parameter.pindex = rc
        return parameter

    def make_vector(self, frequency_vector, gamma_vector):
        """
        Make a frequency-dependent S parameter with values in
        gamma_vector.  The frequencies given here must cover the entire
        span of the calibration frequency range, but do not have to
        coincide with the calibration frequencies -- the function uses
        rational function interpolation as needed to interpolate between
        points.
        """
        cdef vnacal_t *vcp = self.vcp
        frequency_vector = np.asarray(frequency_vector, dtype=np.double,
                                      order="C")
        if frequency_vector.ndim != 1:
            raise ValueError("frequency_vector must be a one-dimensional "
                             "array")
        gamma_vector = np.asarray(gamma_vector, dtype=np.complex128, order="C")
        if gamma_vector.ndim != 1:
            raise ValueError("gamma_vector must be a one-dimensional array")
        if len(gamma_vector) != len(frequency_vector):
            raise ValueError("gamma_vector must be same length as "
                             "frequency_vector")
        cdef const double [::1] fv_view = frequency_vector
        cdef const double complex [::1] gv_view = gamma_vector
        cdef int rc = vnacal_make_vector_parameter(vcp, &fv_view[0],
                                                   len(frequency_vector),
                                                   &gv_view[0])
        self._handle_error(rc)
        cdef Parameter parameter = Parameter.__new__(Parameter)
        parameter.calset = self
        parameter.pindex = rc
        return parameter

    def make_unknown(self, other):
        """
        Make an unknown S parameter that the Solver must determine
        using the given starting value.
        """
        cdef vnacal_t *vcp = self.vcp
        cdef Parameter c_other = self._generalize_parameter(other)
        cdef int rc = vnacal_make_unknown_parameter(vcp, c_other.pindex)
        self._handle_error(rc)
        cdef Parameter parameter = Parameter.__new__(Parameter)
        parameter.calset = self
        parameter.pindex = rc
        return parameter

    def make_correlated(self, other, frequency_vector, sigma_vector):
        """
        Make an unknown S parameter that is known to be correlated
        with another (possibly unknown) parameter with per-frequency
        standard deviation of the difference given by frequency_vector
        and sigma_vector.  The frequencies given here must cover the
        entire span of the calibration frequency range, but do not have
        to coincide with the calibration frequencies -- the function
        uses natural cubic spline interpolation as needed to interpolate
        between points.
        """
        cdef vnacal_t *vcp = self.vcp
        cdef Parameter c_other = self._generalize_parameter(other)
        frequency_vector = np.asarray(frequency_vector, dtype=np.double,
                                      order="C")
        if frequency_vector.ndim != 1:
            raise ValueError("frequency_vector must be a one-dimensional "
                             "array")
        sigma_vector = np.asarray(sigma_vector, dtype=np.complex128, order="C")
        if sigma_vector.ndim != 1:
            raise ValueError("sigma_vector must be a one-dimensional array")
        if len(sigma_vector) != len(frequency_vector):
            raise ValueError("sigma_vector must be same length as "
                             "frequency_vector")
        cdef const double [::1] fv_view = frequency_vector
        cdef const double [::1] sv_view = sigma_vector
        cdef int rc = vnacal_make_correlated_parameter(vcp, c_other.pindex,
                                                       &fv_view[0],
                                                       len(frequency_vector),
                                                       &sv_view[0])
        self._handle_error(rc)
        cdef Parameter parameter = Parameter.__new__(Parameter)
        parameter.calset = self
        parameter.pindex = rc
        return parameter

    @property
    def calibrations(self):
        """
        Vector of calibrations
        """
        helper = _CalHelper()
        helper.calset = self
        return helper

    @property
    def properties(self):
        """
        Reading from this property returns a tree of dictionaries,
        lists, scalars and None representing the global properties.
        Writing this this property replaces the global property tree.
        """
        cdef vnacal_t *vcp = self.vcp
        cdef vnaproperty_t *root = vnacal_property_get_subtree(vcp, -1, ".")
        return _property_to_py(root)

    @properties.setter
    def properties(self, value):
        # no docstring for setter
        cdef vnacal_t *vcp = self.vcp
        cdef vnaproperty_t **rootptr
        cdef vnaproperty_t *new_root = NULL
        success = False
        try:
            _py_to_property(value, &new_root)
            rootptr = vnacal_property_set_subtree(vcp, -1, ".")
            rootptr[0] = new_root
            new_root = NULL
            success = True
        finally:
            if not success:
                vnaproperty_delete(&new_root, ".")


    ######################################################################
    # Internal Functions
    ######################################################################

    cdef _handle_error(self, int rc):
        # """
        # Check if _error_fn has saved an exception for this thread or
        # if rc is -1.  In either case, raise an exception or warning.
        #
        # Parameters:
        #    self:  vna.data.Data class reference
        #    rc:    return value from C function
        #
        # Raises:
        #    See exceptions in _error_fn.
        # """
        exception = self._thread_local._vna_cal_exception
        warning = self._thread_local._vna_cal_warning
        self._thread_local._vna_cal_exception = None
        self._thread_local._vna_cal_warning = None
        if rc == -1 and exception is None:
            PyErr_SetFromErrno(OSError)
        if exception is not None:
            raise exception
        if warning is not None:
            warnings.warn(warning)
