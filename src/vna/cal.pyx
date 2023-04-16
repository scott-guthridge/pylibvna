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
cimport numpy as np
from threading import local
import warnings
import vna.data
from vna.data import Data

np.import_array()   # otherwise, PyArray_SimpleNewFromData segfaults

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
                    vnaerr_category_t category):
    # """
    # C callback function for vnaerr
    # """
    self = <CalSet>error_arg
    if self._thread_local.exception is not None:
        return
    umessage = (<bytes>message).decode("utf8")
    if   category == VNAERR_SYSTEM:
        if errno.errorcode == errno.ENOMEM:
            self._thread_local.exception = MemoryError(umessage)
        else:
            self._thread_local.exception = OSError(errno.errorcode, umessage)
    elif category == VNAERR_USAGE:
        self._thread_local.exception = ValueError(umessage)
    elif category == VNAERR_VERSION:
        self._thread_local.exception = ValueError(umessage)
    elif category == VNAERR_SYNTAX:
        self._thread_local.exception = SyntaxError(umessage, None)
    elif category == VNAERR_WARNING:
        if self._thread_local.warning is None:
            self._thread_local.warning = umessage
    elif category == VNAERR_MATH:
        self._thread_local.exception = ArithmeticError(umessage)
    elif category == VNAERR_INTERNAL:
        self._thread_local.exception = AssertionError(umessage)
    else:
        self._thread_local.exception = Exception(umessage)


cdef class Parameter:
    """
    A measured or defined element of an S-parameter matrix.
    """
    cdef CalSet calset
    cdef int ci

    def __cinit__(self, calset, ci: int):
        self.calset = calset
        self.ci = ci

    cdef get_value(self, double frequency):
        cdef vnacal_t *vcp = self.calset.vcp
        cdef complex result
        result = vnacal_get_parameter_value(vcp, self.ci, frequency)
        self.calset._handle_error(0)
        return result


cdef class Solver:
    cdef vnacal_new_t *vnp

    def __cinit__(self, ctype: CalType, rows, columns, frequencies : int):
        self.vnp = NULL

    #ZZ several methods needed here


cdef class Calibration:
    cdef CalSet calset
    cdef int ci

    @property
    def name(self):
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
    def ctype(self):
        """
        Type of the calibration
        """
        cdef vnacal_t *vcp = self.calset.vcp
        cdef int ci = self.ci
        cdef vnacal_type_t ctype = vnacal_get_type(vcp, ci)
        return <CalType>ctype

    @property
    def rows(self):
        """
        Number of rows in the calibration
        """
        cdef vnacal_t *vcp = self.calset.vcp
        cdef int ci = self.ci
        cdef int rows = vnacal_get_rows(vcp, ci)
        assert(rows != -1)
        return rows

    @property
    def columns(self):
        """
        Number of columns in the calibration
        """
        cdef vnacal_t *vcp = self.calset.vcp
        cdef int ci = self.ci
        cdef int columns = vnacal_get_columns(vcp, ci)
        assert(columns != -1)
        return columns

    @property
    def frequencies(self):
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
        Reference frequency all all ports in the calibration
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
        #
        # Get calibration type and dimensions
        #
        cdef vnacal_t *vcp = self.calset.vcp
        cdef int ci = self.ci
        cdef CalType cal_type = <CalType>vnacal_get_type(vcp, ci)
        assert(<int>cal_type != -1)
        cdef int cal_frequencies = vnacal_get_frequencies(vcp, ci)
        assert(cal_frequencies != -1)
        cdef int cal_rows = vnacal_get_rows(vcp, ci)
        assert(cal_rows != -1)
        cdef int cal_columns = vnacal_get_columns(vcp, ci)
        assert(cal_columns != -1)
        cdef int ports = max(cal_rows, cal_columns)

        #
        # Define vnacal_apply parameters
        #
        cdef int frequencies
        cdef const double *frequency_vector
        cdef double complex **a_clfpp = NULL
        cdef int a_rows = 0
        cdef int a_columns = 0
        cdef double complex **b_clfpp = NULL
        cdef vnadata_t *vdp

        #
        # Define temporary variables and views
        #
        cdef int i
        cdef int j
        cdef int k
        cdef int rc
        cdef const double [::1] v1
        cdef double complex [:, :, ::1] v3

        #
        # Process frequency vector.  If f is None, default to the
        # calibration frequencies.
        #
        if f is None:
            frequencies = cal_frequencies
            frequency_vector = vnacal_get_frequency_vector(vcp, ci)
            assert(frequency_vector != NULL)
        else:
            f_array = np.asarray(f, dytpe=np.double, order="C")
            if f_array.ndim != 1:
                raise ValueError("f must be a one-dimensional array or None")
            frequencies = f_array.shape[0]
            v1 = f_array
            frequency_vector = &v1[0]

        #
        # Prepare a_array
        #
        if a is not None:
            if cal_type == E12 or cal_type == UE14:
                a_rows = 1
            else:
                a_rows = ports

            a_columns = ports
            a_array = np.asarray(a, dtype=np.complex128, order="C")
            if a_array.ndim != 3 \
                    or a_array.shape != (a_rows, ports, frequencies):
                raise ValueError(f"a must be a {a_rows} x {ports} x "
                                 f"{frequencies} complex matrix or None")

        #
        # Prepare b_array
        #
        if b is None:
            raise ValueError(f"b must be a {ports} x "
                             f"{ports} x {frequencies} array")
        b_array = np.asarray(b, dtype=np.complex128, order="C")
        if b_array.ndim != 3 or b_array.shape != (ports, ports, frequencies):
            raise ValueError(f"b must have dimension {ports} x "
                             f"{ports} x {frequencies}")

        #
        # Allocate C arrays
        #
        b_clfpp = <double complex **>malloc(
                ports * ports * sizeof(double complex *))
        if b_clfpp == NULL:
            raise MemoryError()

        if a is not None:
            a_clfpp = <double complex **>malloc(
                    a_rows * ports * sizeof(double complex *))
            if a_clfpp == NULL:
                free(b_clfpp)
                raise MemoryError()

        #
        # Avoid leaking the C arrays...
        #
        try:
            #
            # Fill C arrays
            #
            if a is not None:
                v3 = a_array
                k = 0
                for i in range(a_rows):
                    for j in range(ports):
                        a_clfpp[k] = &v3[i, j, 0]
                        k += 1
            v3 = b_array
            k = 0
            for i in range(ports):
                for j in range(ports):
                    b_clfpp[k] = &v3[i, j, 0]
                    k += 1

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
                                  b_clfpp, ports, ports,
                                  vdp)
            else:
                rc = vnacal_apply_m(vcp, ci,
                                  frequency_vector, frequencies,
                                  b_clfpp, ports, ports,
                                  vdp)
            self.calset._handle_error(rc)

            return result

        finally:
            free(b_clfpp)
            free(a_clfpp)

    # ZZ: add getter and setter for properties


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
        for i in range(len(calset._index_to_ci)):
            yield self[i]

    def  __reversed__(self):
        """
        Iterate reversed over the Calibrations
        """
        cdef CalSet calset = self.calset
        for i in reversed(range(len(calset._index_to_ci))):
            yield self[i]


cdef class CalSet:
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
        self._thread_local.exception = None
        self._thread_local.warning = None
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

    def add(self, solver: Solver, name: str) -> int:
        """
        Add a new solved calibration to the CalSet.
        """
        cdef vnacal_t *vcp = self.vcp
        cdef vnacal_new_t *vnp = solver.vnp
        cdef int rc = vnacal_add_calibration(vcp, name, vnp)
        self._handle_error(rc)
        self._index_to_ci.append(rc)

    def index(self, name: str) -> int:
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

    @property
    def calibrations(self):
        """
        Vector of calibrations
        """
        helper = _CalHelper()
        helper.calset = self
        return helper


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
        exception = self._thread_local.exception
        warning = self._thread_local.warning
        self._thread_local.exception = None
        self._thread_local.warning = None
        if rc == -1 and exception is None:
            PyErr_SetFromErrno(OSError)
        if exception is not None:
            raise exception
        if warning is not None:
            warnings.warn(warning)
