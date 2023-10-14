#cython: language_level=3
#cython: binding=True
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

"""
This module uses measurements of calibration standards to find error
terms for vector network analyzers (VNAs), then applies the resulting
calibration to measurements of devices under test to correct for the
errors.

The basic sequence for generating a new calibration is:

#. Create a *Calset* object, then a *Solver* object.
#. Make measurements of calibration standards and add them using the
   Solver.add_* methods.
#. Solve using `Solver.solve()`
#. Use `Solver.add_to_calset()` to add the new calibration to the Calset
#. Use `Calset.save()` to save the calibration to a file

The sequence for applying a calibration to a device measurement is:

#. Load the calibration by passing the calibration filename to the
   *Calset* constructor.
#. Get the *Calibration* object from the *Calset.calibrations* attribute.
#. Measure the device under test.
#. Use `Calibration.apply()` to apply the calibration correction to
   the measurements.
"""

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
import libvna.data
from libvna.data import NPData

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
    self = <Calset>error_arg
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
    An element of the S-parameter matrix for a calibration standard.
    The Parameter is an abstraction that can represent a single complex
    scalar such as -1 for short, a vector of (frequency, gamma) tuples
    representing a reflect with complex impedance or the through component
    of a transmission line, or an unknown parameter the library must
    solve, e.g. the R or L parameters in TRL.

    Note: this class cannot be instantiated directly: use
    :class:`ScalarParameter`, :class:`VectorParameter`,
    :class:`UnknownParameter` and :class:`CorrelatedParameter`.
    """
    cdef Calset calset
    cdef int pindex

    def __cinit__(self):
        self.calset = None
        self.pindex = -1

    def __init__(self):
        raise TypeError("Can't instantiate abstract class Parameter")

    def __dealloc__(self):
        cdef Calset calset = self.calset
        cdef pindex = self.pindex
        if calset is not None and pindex >= 3:
            vnacal_delete_parameter(calset.vcp, pindex)

    @staticmethod
    cdef Parameter _generalize(Calset calset, parameter):
        # """
        # If a number is found where a Parameter is expected, automatically
        # convert it into a scalar parameter.  If a tuple(f_vector, g_vector)
        # is found, automatically convert it to a vector parameter.
        # """
        assert calset is not None
        cdef double complex value
        if isinstance(parameter, Parameter):
            return parameter
        if isinstance(parameter, complex) or isinstance(parameter, float) \
                or isinstance(parameter, int):
            return ScalarParameter(calset, parameter)
        if isinstance(parameter, tuple):
            if len(parameter) != 2:
                raise ValueError("expected Parameter, number or "
                                 "tuple(frequency_vector, gamma_vector)")
            return VectorParameter(calset, *parameter)
        raise ValueError("parameter must be class Parameter, a number, or "
                         "tuple(frequency_vector, gamma_vector)")

    def get_value(self, frequency) -> complex:
        """
        Parameters:
            frequency (float):
                Frequency at which to evaluate the value.

        Return the value of the parameter at a given frequency.
        For a scalar parameter, the function ignores *frequency*
        and simply returns the fixed gamma value.  For a vector
        parameter, it returns the gamma value at the given frequency,
        interpolating as necessary.  If the parameter is unknown
        and :func:`Solver.solve` has completed successfully,
        :func:`get_value` returns the solved value, again
        interpolating as necessary.
        """
        cdef vnacal_t *vcp = self.calset.vcp
        cdef complex result
        result = vnacal_get_parameter_value(vcp, self.pindex, frequency)
        self.calset._handle_error(0)
        return result


cdef class ScalarParameter(Parameter):
    """
    An element of the S parameter matrix of a calibration standard that
    has a constant value at all frequencies, e.g. -1 for short.

    Parameters:
        calset (Calset):
            The associated calibration set.
        gamma (complex):
            an element of the S parameter matrix of the standard
            that doesn't depend on frequency, e.g. 0 for match.
    """
    def __cinit__(self, Calset calset, gamma):
        if calset is None:
            raise ValueError("calset cannot be None")
        cdef vnacal_t *vcp = calset.vcp
        cdef int rc
        if gamma == 0.0:        # special-case these values
            rc = VNACAL_MATCH
        elif gamma == +1.0:
            rc = VNACAL_OPEN
        elif gamma == -1.0:
            rc = VNACAL_SHORT
        else:
            rc = vnacal_make_scalar_parameter(vcp, gamma)
            calset._handle_error(rc)
        self.calset = calset
        self.pindex = rc

    def __init__(self, *args, **kwargs):
        pass


cdef class VectorParameter(Parameter):
    """
    An element of the S parameter matrix of a calibration standard
    that varies with frequency.

    Parameters:
        calset (Calset):
            The associated calibration set.
        frequency_vector (vector of float):
            monotonically increasing list of frequencies
        gamma_vector (vector of complex):
            list or array of complex values corresponding
            to each frequency in *frequency_vector*

    The frequencies must cover the entire span of the calibration
    frequency range, but do not have to coincide with the calibration
    frequencies -- the library uses rational function interpolation
    as needed to interpolate between frequency points.
    """
    def __cinit__(self, Calset calset, frequency_vector, gamma_vector):
        if calset is None:
            raise ValueError("calset cannot be None")
        cdef vnacal_t *vcp = calset.vcp
        frequency_vector = np.asarray(frequency_vector, dtype=np.double,
                                      order="C")
        if frequency_vector.ndim != 1:
            raise ValueError("frequency_vector must be a one-dimensional "
                             "array")
        gamma_vector = np.asarray(gamma_vector, dtype=np.complex128,
                                  order="C")
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
        calset._handle_error(rc)
        self.calset = calset
        self.pindex = rc

    def __init__(self, *args, **kwargs):
        pass


cdef class UnknownParameter(Parameter):
    """
    An element of the S parameter matrix of a calibration standard that
    is only approximately known and that the library must determine.

    Parameters:
        calset (Calset):
            The associated calibration set.
        initial_guess (complex, list of (frequency, gamma) tuples or Parameter):
            Approximate value of the unknown parameter
    """
    def __cinit__(self, Calset calset, initial_guess):
        if calset is None:
            raise ValueError("calset cannot be None")
        cdef vnacal_t *vcp = calset.vcp
        cdef Parameter c_other = Parameter._generalize(calset, initial_guess)
        cdef int rc = vnacal_make_unknown_parameter(vcp, c_other.pindex)
        calset._handle_error(rc)
        self.calset = calset
        self.pindex = rc

    def __init__(self, *args, **kwargs):
        pass


cdef class CorrelatedParameter(Parameter):
    """
    An element of the S parameter matrix of a calibration standard
    that is known to be correlated with another (possibly unknown)
    Parameter.  This type of parameter is useful for modeling connection
    non-repeatability.

    Parameters:
        calset (Calset):
            The associated calibration set.
        other (Parameter):
            Another Parameter to which this
            Parameter is known to be correlated
        frequency_vector (vector of float):
            monotonically increasing list of frequencies
        sigma_vector (vector of float):
            standard deviations of the differences between
            this Parameter and its correlate at
            each frequency

    The frequencies must cover at least the entire span of the calibration
    frequency range, but do not have to coincide with the calibration
    frequencies -- the library uses natural cubic spline interpolation
    as needed to interpolate between points.
    """
    def __cinit__(self, Calset calset, other, frequency_vector, sigma_vector):
        if calset is None:
            raise ValueError("calset cannot be None")
        cdef vnacal_t *vcp = calset.vcp
        cdef Parameter c_other = Parameter._generalize(calset, other)
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
        calset._handle_error(rc)
        self.calset = calset
        self.pindex = rc

    def __init__(self, *args, **kwargs):
        pass


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
    """
    Error Term Solver: solve for VNA error terms from measurements
    of calibration standards.

    Args:
        calset (Calset):
            The associated calibration set.

        ctype (CalType):
            The calibration type determines which error terms the
            library corrects.  Valid values are:

            T8, U8:
                8-term T or U parameters: correct for directivity,
                reflection / transmission tracking, and port
                match errors on each VNA port.  At least three
                standards, (e.g. short-open, short-match, through)
                are needed to solve the 2x2 T8 or U8 calibration.

            TE10, UE10:
                8-term T or U parameters plus 2 leakage: correct
                for the same errors as T8 and U8, but also correct
                for leakage within the VNA from the driving port
                to the other ports.  At least three standards
                (e.g. short-open, short-match, through) are needed
                to solve the 2x2 TE10 or UE10 calibration.

            T16, U16:
                16-term T or U parameters: correct for the same
                errors as TE10 and UE10, but add the remaining
                leakage terms including leakage between the DUT
                ports in the test fixture.  At least five standards
                (e.g. short-open, short-match, open-match, open-short,
                through) are needed to solve the 2x2 T16 or U16
                calibration.

            UE14:
                Correct for the same errors as TE10 and UE10, except treat
                each column (driving port) as an independent calibration.
                This type produces separate error parameters for the
                forward and reverse directions, and for this reason,
                it's able to correct for errors in the forward-reverse
                switch without reference (*a* matrix) measurements,
                even for a switch that lies between the detectors and
                the DUT.  At least four standards (e.g. short-open,
                match-open, match-short, through) are needed to solve
                the 2x2 UE14 calibration.

            E12:
                Generalization of classic SOLT: the library uses UE14
                terms internally to solve this calibration and this type
                corrects for exactly the same errors.  The difference
                is only in the representation of the saved error terms:
                after finding the error terms, the library converts
                from inverse scattering transfer (U) to scattering (E)
                at the end of the calibration procedure.  At least four
                standards (e.g. short-open, match-open, match-short,
                through) are needed to solve the 2x2 E12 calibration.

        rows (int):
            Number of rows in the calibration, where *rows* is the
            number of VNA ports that detect signal.  Normally, *rows*
            and *columns* are both simply the number of VNA ports.
            Some simple VNAs, however, measure only a subset of the
            S parameters such as :math:`S_{11}` and :math:`S_{21}`.
            If the VNA measures :math:`S_{11}` and :math:`S_{21}`
            only, set *rows* to 2 and *columns* to 1.  If the VNA
            measures :math:`S_{11}` and :math:`S_{12}` only, set
            *rows* to 1 and *columns* to 2.  In general, if *rows* >
            *columns*, U parameters must be used; if *rows* < *columns*,
            T parameters must be used.  For square calibrations, either
            T or U parameters may be used.

        columns (int):
            Number columns in the calibration, where *columns* is
            the number of VNA ports that transmit signal.  See *rows*.

        frequency_vector (list/array of float):
            Vector of frequency points to be used in the calibration.
            Must be monotonically increasing.

        z0 (complex, optional):
            Reference impedance of the VNA ports.  All ports must
            have the same reference impedance.  If not specified,
            *z0* defaults to 50 ohms.

    Note that the calibration method used, e.g. LRL, LRM, LRRL, LRRM,
    LXYZ, SOLT, TRD, TRL, TRM, TXYZ, UXYZ, etc., does not have to be
    specified.  The library automatically determines the method based
    on the standards given, or uses a general solver if it does not
    have a special solver for the given set of standards.  The main
    requirement is that the standards used must provide a sufficient
    number of conditions for the number of unknowns to be solved.
    """
    cdef Calset calset
    cdef int frequencies
    cdef vnacal_new_t *vnp
    cdef double pvalue_limit
    cdef double et_tolerance
    cdef double p_tolerance
    cdef int iteration_limit

    def __cinit__(self, Calset calset, CalType ctype, int rows, int columns,
                  frequency_vector, double complex z0 = 50.0):
        if calset is None:
            raise ValueError("calset cannot be None")

        frequency_vector = np.asarray(frequency_vector, dtype=np.double,
                                      order="C")
        if frequency_vector.ndim != 1:
            raise ValueError("frequency_vector must be 1 dimensional "
                             "array-like")
        cdef vnacal_t *vcp = calset.vcp
        cdef vnacal_new_t *vnp = NULL
        vnp = vnacal_new_alloc(vcp, <vnacal_type_t>ctype, rows, columns,
                               len(frequency_vector))
        if vnp == NULL:
            calset._handle_error(-1)
        cdef int rc
        cdef const double [::1] fv_view = frequency_vector
        rc = vnacal_new_set_frequency_vector(vnp, &fv_view[0])
        if rc == -1:
            vnacal_new_free(vnp)
            calset._handle_error(-1)
        if z0 != 50.0:
            rc = vnacal_new_set_z0(vnp, z0)
            if rc == -1:
                vnacal_new_free(vnp)
                calset._handle_error(-1)
        self.calset = calset
        self.frequencies = len(frequency_vector)
        self.vnp = vnp
        self.pvalue_limit = 0.001
        self.et_tolerance = 1.0e-6
        self.p_tolerance = 1.0e-6
        self.iteration_limit = 30

    def __dealloc__(self):
        # Note that the reference on calset is released after this
        # function returns.
        vnacal_new_free(self.vnp)
        self.vnp = NULL

    def add_single_reflect(self, a, b, s11, int port=1):
        """
        Add the measurement of a single reflect standard with parameter
        *s11* on the given VNA port.

        Parameters:
            a (matrix of vectors of complex):
                incident root power into each DUT port, or None if
                not available
            b (matrix of vectors of complex):
                reflected root power from each DUT port
            s11 (complex, list of (frequency, gamma) tuples, or Parameter):
                :math:`S_{11}` parameter of the the calibration standard
            port (int, optional):
                VNA port number connected to the standard.  If not given,
                defaults to 1.
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
            c_s11 = Parameter._generalize(self.calset, s11)

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

    def add_double_reflect(self, a, b, s11, s22, int port1=1, int port2=2):
        """
        Add the measurement of a double reflect standard with parameters
        *s11* and *s22* on the given VNA ports, assuming :math:`S_{12} =
        S_{21} = 0`.

        Parameters:
            a (matrix of vectors of complex):
                incident root power into each DUT port, or None if
                not available
            b (matrix of vectors of complex):
                reflected root power from each DUT port
            s11 (complex, list of (frequency, gamma) tuples, or Parameter):
                the :math:`S_{11}` parameter of the the calibration standard
            s22 (complex, list of (frequency, gamma) tuples, or Parameter):
                the :math:`S_{22}` parameter of the calibration standard
            port1 (int, optional):
                VNA port number connected to port 1 of the calibration
                standard.  If not given, defaults to 1.
            port2 (int, optional):
                VNA port number connected to port 2 of the calibration
                standard.  If not given, defaults to 2.
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
            c_s11 = Parameter._generalize(self.calset, s11)
            c_s22 = Parameter._generalize(self.calset, s22)

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

    def add_through(self, a, b, int port1=1, int port2=2):
        """
        Add the measurement of a perfect through standard between *port1*
        and *port2*, i.e.
        :math:`S_{12} = S_{21} = 1` and :math:`S_{11} = S_{22} = 0`.

        Parameters:
            a (matrix of vectors of complex):
                incident root power into each DUT port, or None if
                not available
            b (matrix of vectors of complex):
                reflected root power from each DUT port
            port1 (int, optional):
                First VNA port connected to the through standard.
                If not given, defaults to 1.
            port2 (int, optional):
                Second VNA port connected to the through standard.
                If not given, defaults to 2.
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

    def add_line(self, a, b, s, int port1=1, int port2=2):
        """
        Add the measurement of an arbitrary two-port standard with S
        parameter matrix, *s*, on the given VNA ports.

        Parameters:
            a (matrix of vectors of complex):
                incident root power into each DUT port, or None if
                not available
            b (matrix of vectors of complex):
                reflected root power from each DUT port
            s (2x2 matrix):
                S-parameter matrix of the standard, where each element
                of the matrix can be a complex scalar, a vector of
                (frequency, gamma) tuples, or a Parameter
            port1 (int, optional):
                VNA port number connected to port 1 of the calibration
                standard.  If not given, defaults to 1.
            port2 (int, optional):
                VNA port number connected to port 2 of the calibration
                standard.  If not given, defaults to 2.
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
                    c_parameter = Parameter._generalize(self.calset, s[i, j])
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

    def add_mapped_matrix(self, a, b, s, port_map = None):
        """
        Add the measurement of an arbitrary n-port standard with S
        parameter matrix, *s*, and a map of ports of the standard to
        ports of the VNA in *port_map*.

        Parameters:
            a (matrix of vectors of complex):
                incident root power into each DUT port, or None if
                not available
            b (matrix of vectors of complex):
                reflected root power from each DUT port
            s (matrix):
                S-parameter matrix of the standard, where each element
                of the matrix can be a complex scalar, a vector of
                (frequency, gamma) tuples, or a Parameter
            port_map (vector of int, optional):
                List of the VNA port numbers attached to each port of
                the standard in order.  Optional if the standard has
                the same number of ports as the VNA and the ports of
                the VNA are attached to the corresponding port numbers
                of the standard.  VNA port numbers start with 1.
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
        cdef int *mip = NULL
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
                    c_parameter = Parameter._generalize(self.calset, s[i, j])
                    s[i, j] = c_parameter
                    sip[k] = c_parameter.pindex
                    k += 1

            #
            # Prepare port map.
            #
            if port_map is not None:
                port_map = np.asarray(port_map, dtype=np.intc, order="C")
                if port_map.ndim != 1 or port_map.shape[0] < s_ports:
                    raise ValueError(f"port_map must be a length {s_ports} "
                                     f"vector")
                iv = port_map
                mip = &iv[0]

            #
            # Call add function
            #
            if a_clfpp != NULL:
                rc = vnacal_new_add_mapped_matrix(self.vnp,
                                                  a_clfpp, a_rows, a_columns,
                                                  b_clfpp, b_rows, b_columns,
                                                  sip, s_rows, s_columns,
                                                  mip)
            else:
                rc = vnacal_new_add_mapped_matrix_m(self.vnp,
                                                    b_clfpp, b_rows, b_columns,
                                                    sip, s_rows, s_columns,
                                                    mip)
            self.calset._handle_error(rc)
            return

        finally:
            free(<void *>b_clfpp)
            free(<void *>a_clfpp)
            free(<void *>sip)

    def set_m_error(self, frequency_vector, noise_floor,
                    tracking_error = None):
        """
        Enable measurement error modeling.

        Parameters:
            frequency_vector (vector of float):
                vector of ascending frequencies spanning at least the
                full calibration frequency range

            noise_floor (vector of float):
                vector of noise floor root-power measurements at the
                VNA detectors at each frequency when no signal is applied

            tracking_error (vector of float, optional):
                optional vector describing an additional root-power
                noise source proportional to the amplitude of the
                measured signal

        Both noise sources are assumed to be Gaussian and independent.
        Specifying measurement errors with this function can significanlty
        improve accuracy, especially for significantly overdetermined
        systems, as the 16 term models typically are.
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
        p-value, below which to reject the null hypothesis that the
        measurement errors are less than or equal to the values given
        in set_m_error (float).

        This parameter has no effect if measurement error modeling has not
        been enabled through :func:`set_m_error`.  The default is 0.001.
        """
        return self.pvalue_limit

    @pvalue_limit.setter
    def pvalue_limit(self, double value):
        # no docstring for setter
        cdef int rc
        rc = vnacal_new_set_pvalue_limit(self.vnp, value)
        self.calset._handle_error(rc)
        self.pvalue_limit = value

    @property
    def et_tolerance(self):
        """
        For iterative solution methods, this parameter controls the
        degree of change in the root-mean-squared of the error terms
        sufficiently small to stop iteration (float).  Default is 1.0e-6.
        """
        return self.et_tolerance

    @et_tolerance.setter
    def et_tolerance(self, double value):
        # no docstring for setter
        cdef int rc
        rc = vnacal_new_set_et_tolerance(self.vnp, value)
        self.calset._handle_error(rc)
        self.et_tolerance = value

    @property
    def p_tolerance(self):
        """
        For iterative solution methods, this parameter controls the degree
        of change in the root-mean-squared of the unknown parameters
        sufficiently small to stop iteration (float).  This parameter
        has no effect if there are no unknown parameters in the S matrix.
        Default is 1.0e-6.
        """
        return self.et_tolerance

    @p_tolerance.setter
    def p_tolerance(self, double value):
        # no docstring for setter
        cdef int rc
        rc = vnacal_new_set_p_tolerance(self.vnp, value)
        self.calset._handle_error(rc)
        self.p_tolerance = value

    @property
    def iteration_limit(self):
        """
        For iterative solution methods, this parameter controls the
        maximum number of iterations permitted to reach convergence (int).
        The default is 30.
        """
        return self.iteration_limit

    @iteration_limit.setter
    def iteration_limit(self, int value):
        # no docstring for setter
        cdef int rc
        rc = vnacal_new_set_iteration_limit(self.vnp, value)
        self.calset._handle_error(rc)
        self.iteration_limit = value

    def solve(self):
        """
        Solve for the error terms.  Note: if this function raises an
        exception due to an insufficient number of standards, you may
        add additional standards and try again.
        """
        cdef int rc = vnacal_new_solve(self.vnp)
        self.calset._handle_error(rc)

    def add_to_calset(self, name) -> int:
        """
        Add the solved calibration to the Calset.  If *name* matches
        an existing calibration, the existing calibration is replaced;
        if *name* is unique, then the new calibration is appended to
        the *Calset.calibrations* array.

        Parameters:
            name:
                Name for the calibration.

        Returns:
            Index of the new entry in Calset.calibrations
        """
        cdef vnacal_t *vcp = self.calset.vcp
        cdef vnacal_new_t *vnp = self.vnp
        if isinstance(name, unicode):
            name = (<unicode>name).encode("UTF-8")
        cdef const unsigned char[:] c_name = name
        cdef int rc = vnacal_add_calibration(vcp, <const char *>&c_name[0],
                                             vnp)
        self.calset._handle_error(rc)
        self.calset._index_to_ci.append(rc)
        return len(self.calset._index_to_ci) - 1


cdef class Calibration:
    """
    A solved calibration.  Note: you cannot instantiate this class
    directly: use ``Calset.calibrations`` to access instances of this
    class.
    """
    cdef Calset calset
    cdef int ci

    @property
    def name(self) -> str:
        """
        Name of the calibration (str, readonly)
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
        Type of the calibration (CalType, readonly)
        """
        cdef vnacal_t *vcp = self.calset.vcp
        cdef int ci = self.ci
        cdef vnacal_type_t ctype = vnacal_get_type(vcp, ci)
        return <CalType>ctype

    @property
    def rows(self) -> int:
        """
        Number of rows in the calibration (int, readonly)
        """
        cdef vnacal_t *vcp = self.calset.vcp
        cdef int ci = self.ci
        cdef int rows = vnacal_get_rows(vcp, ci)
        assert(rows != -1)
        return rows

    @property
    def columns(self) -> int:
        """
        Number of columns in the calibration (int, readonly)
        """
        cdef vnacal_t *vcp = self.calset.vcp
        cdef int ci = self.ci
        cdef int columns = vnacal_get_columns(vcp, ci)
        assert(columns != -1)
        return columns

    @property
    def frequencies(self) -> int:
        """
        Number of frequencies in the calibration (int, readonly)
        """
        cdef vnacal_t *vcp = self.calset.vcp
        cdef int ci = self.ci
        cdef int frequencies = vnacal_get_frequencies(vcp, ci)
        assert(frequencies != -1)
        return frequencies

    @property
    def frequency_vector(self):
        """
        Vector of calibration frequencies (array of float, readonly)
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
    def z0(self) -> complex:
        """
        Reference frequency of all ports in the calibration (complex,
        readonly)
        """
        cdef vnacal_t *vcp = self.calset.vcp
        cdef int ci = self.ci
        cdef double complex z0 = vnacal_get_z0(vcp, ci)
        return z0

    def apply(self, f, a, b) -> NPData:
        """
        Apply the calibration correction to measured data.  The
        calibration must have dimensions 2x1, 1x2, or NxN.

        Parameters:
            f (array of float or None):
                vector of frequencies at which the measurements were made,
                or None if measured at the calibration frequencies
            a (matrix of vectors of complex):
                incident root power into each DUT port, or None if
                not available
            b (matrix of vectors of complex):
                reflected root power from each DUT port

        Return:
            libvna.data.NPData object containing the corrected parameters
               
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
            result = libvna.data.NPData()
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
        A tree of nested dictionaries, lists, scalars and None values
        representing arbitrary user-defined properties of the calibration.
        Some examples are the calibration date, cable lengths, and test
        set used.
        """
        cdef int ci = self.ci
        return self.calset._get_properties(ci)

    @properties.setter
    def properties(self, value):
        cdef int ci = self.ci
        self.calset._set_properties(ci, value)

    def __str__(self):
        return f"name=\"{self.name}\", ctype={self.ctype.name}, " \
               f"rows={self.rows}, columns={self.columns}, " \
               f"fmin={self.frequency_vector[0]:.3e}, " \
               f"fmax={self.frequency_vector[-1]:.3e}"


cdef class _CalHelper:
    # """
    # Internal class that provides __getitem__, __delitem__, etc on
    # calibrations.
    # """
    cdef Calset calset

    def _find_ci_by_index_or_name(self, index):
        # """
        # Find the vnacal calibration index (ci) from the dense
        # index or name.
        # """
        cdef Calset calset = self.calset
        cdef vnacal_t *vcp = calset.vcp
        cdef Calibration calibration
        if isinstance(index, int):
            return calset._index_to_ci[index]

        if isinstance(index, str):
            return calset.index(index)

        raise IndexError("index must have type int or str")

    def __getitem__(self, index) -> Calibration:
        """
        Return a Calibration object by name or index.
        """
        cdef Calset calset = self.calset
        cdef int ci = self._find_ci_by_index_or_name(index)
        cdef Calibration calibration
        calibration = Calibration()
        calibration.calset = calset
        calibration.ci = ci
        return calibration

    def __len__(self):
        cdef Calset calset = self.calset
        return len(calset._index_to_ci)

    def __delitem__(self, index):
        """
        Delete a calibration object by name or index.
        """
        cdef Calset calset = self.calset
        cdef vnacal_t *vcp = calset.vcp
        cdef int ci = self._find_ci_by_index_or_name(index)
        cdef int rc = vnacal_delete_calibration(vcp, ci)
        assert(rc == 0)
        del calset._index_to_ci[index]

    def __iter__(self):
        """
        Iterate over the Calibrations
        """
        cdef Calset calset = self.calset
        cdef int i
        for i in range(len(calset._index_to_ci)):
            yield self[i]

    def  __reversed__(self):
        """
        Iterate reversed over the Calibrations
        """
        cdef Calset calset = self.calset
        cdef int i
        for i in reversed(range(len(calset._index_to_ci))):
            yield self[i]

    def __str__(self):
        """
        Return a string representation of the Calibrations array
        """
        result = ""
        cdef int i
        for i, cal in enumerate(self.calset.calibrations):
            if i > 0:
                result += "\n"
            result += str(cal)
        return result


cdef class Calset:
    """
    The Calset loads, saves and manages solved VNA calibrations.  It's
    the central data structure of VNA calibration, needed by other
    classes in the module.

    Though the Calset usually holds just one calibration, it can hold
    any number of related calibrations, for example, covering different
    frequency bands or test set configurations.

    Args:
        filename (str, optional):
            Load the Calset from the given file.  Note that the
            recommended .vnacal extension is not added automatically
            and should be included in filename.
    """
    cdef vnacal_t *vcp
    cdef object _properties
    cdef object _thread_local
    cdef object _index_to_ci

    def __cinit__(self, filename = None):
        cdef int ci
        cdef int ci_end
        cdef vnacal_type_t ctype
        cdef const unsigned char[:] cfilename
        self.vcp = NULL
        self._properties = list()
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

    def _extend_properties(self, int ci):
        # """
        # Grow the _properties list as needed to hold index ci + 1.
        # """
        cdef old_len = len(self._properties)
        if old_len < ci + 2:
            self._properties += [...] * (ci + 2 - old_len)

    def _get_properties(self, int ci):
        # """
        # Populate _properties[ci + 1] with properties from the C vnacal_t
        # structure.  The index is + 1 so that we store the global properties
        # at index 0, the properties for ci=0 at index 1, etc.  The special
        # value, ..., is used to indicate that a given slot doesn't contain
        # valid data.  We use ... instead of None because None is a valid
        # properties value.
        # """
        cdef vnacal_t *vcp = self.vcp
        self._extend_properties(ci)
        cdef vnaproperty_t *root = NULL
        if self._properties[ci + 1] is ...:
            root = vnacal_property_get_subtree(vcp, ci, ".")
            self._properties[ci + 1] = _property_to_py(root)
        return self._properties[ci + 1]

    def _set_properties(self, int ci, value):
        # """
        # Replace the property tree at ci with value.
        # """
        self._extend_properties(ci)
        self._properties[ci + 1] = value

    def _put_properties(self, int ci):
        # """
        # Put the properties stored at _properties[ci + 1] back into
        # the vnacal_t C structure.  The given slot must exist and must
        # not contain an ellipsis.
        # """
        cdef vnacal_t *vcp = self.vcp
        cdef vnaproperty_t **rootptr
        cdef vnaproperty_t *new_root = NULL
        success = False
        try:
            _py_to_property(self._properties[ci + 1], &new_root)
            rootptr = vnacal_property_set_subtree(vcp, ci, ".")
            rootptr[0] = new_root
            new_root = NULL
            success = True
        finally:
            if not success:
                vnaproperty_delete(&new_root, ".")

    def _put_all_properties(self):
        # """
        # Put the properties for each ci value back into the vnacal_t
        # C structure.
        # """
        for index, value in enumerate(self._properties):
            if value is not ...:
                self._put_properties(index - 1)

    def save(self, filename):
        """
        Save the Calset to a file.

        Parameters:
            filename (str):
                Pathname of the save file.  The recommended file
                extension, ".vnacal", is not added automatically and
                should included in *filename*.
        """
        self._put_all_properties()
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode("UTF-8")
        cdef const unsigned char[:] cfilename = filename
        cdef int rc = vnacal_save(self.vcp, <const char *>&cfilename[0])
        self._handle_error(rc)

    def index(self, name) -> int:
        """
        Convert from calibration name to index:

        Parameters:
            name (str):
                Name of the calibration
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

    @property
    def calibrations(self):
        """
        Vector of Calibration objects in this Calset.  This attribute
        is indexable, iterable, and the elements can be deleted.
        """
        helper = _CalHelper()
        helper.calset = self
        return helper

    @property
    def properties(self):
        """
        A tree of nested dictionaries, lists, scalars and None values
        representing arbitrary user-defined global properties of the
        calibration set, for example, the model and serial number of
        the instrument, test fixture used, and other conditions under
        which the calibration was made.

        Also see the calibration-specific properties,
        Calibration.properties.
        """
        return self._get_properties(-1)

    @properties.setter
    def properties(self, value):
        self._set_properties(-1, value)


    ######################################################################
    # Internal Functions
    ######################################################################

    cdef _handle_error(self, int rc):
        # """
        # Check if _error_fn has saved an exception for this thread or
        # if rc is -1.  In either case, raise an exception or warning.
        #
        # Parameters:
        #    self:  libvna.cal.Calset class reference
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
