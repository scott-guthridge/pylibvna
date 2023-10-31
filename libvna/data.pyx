#cython: language_level=3
#cython: binding=True
# Python Bindings for Vector Network Analyzer Library
# Copyright © 2023 D Scott Guthridge <scott_guthridge@rompromity.net>
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from cpython.exc cimport PyErr_SetFromErrno
from cpython.pycapsule cimport PyCapsule_New
import errno
from libc.stdio cimport FILE, fdopen, fclose
import numpy as np
cimport numpy as cnp
from posix.unistd cimport close, dup, lseek, off_t
import re
from threading import local
import warnings


cnp.import_array()   # otherwise, PyArray_SimpleNewFromData segfaults

cpdef enum PType:
    # Note: current cython version does not allow a docstring here.
    # """
    # Select the network parameter data type.
    # Values:
    #
    # - ANY: unspecified data (matrices may be rectangular)
    # - S:   scattering paramters
    # - T:   scattering-transfer parameters
    # - U:   inverse scattering-transfer parameters
    # - Z:   impedance parameters
    # - Y:   admittance parameters
    # - H:   hybrid parameters
    # - G:   inverse hybrid parameters
    # - A:   ABCD parameters
    # - B:   inverse ABCD parameters
    # - ZIN: input impedances (single row-vector)
    # """
    ANY   =  0
    S     =  1
    T     =  2
    U     =  3
    Z     =  4
    Y     =  5
    H     =  6
    G     =  7
    A     =  8
    B     =  9
    ZIN   = 10

cpdef enum FileType:
    # Note: current cython version does not allow a docstring here.
    # """
    # Select the file type used to save network parameter data.
    # Values:
    #     AUTO:        determine the type based on filename extension (default)
    #     TOUCHSTONE1: use Touchstone version 1 format (.s2p)
    #     TOUCHSTONE2: use Touchstone version 2 format (.ts)
    #     NPD:         use network parameter data format (.npd)
    # """
    AUTO        = 0
    TOUCHSTONE1 = 1
    TOUCHSTONE2 = 2
    NPD         = 3

cdef enum IndexClass:
    # IndexClass:
    #     Bitmasks giving which data_array subscripts are integers
    #     vs. slices.
    IDX_SSS = 0x0
    IDX_ISS = 0x1
    IDX_SIS = 0x2
    IDX_IIS = 0x3
    IDX_SSI = 0x4
    IDX_ISI = 0x5
    IDX_SII = 0x6
    IDX_III = 0x7


cdef void _error_fn(const char *message, void *error_arg,
                    vnaerr_category_t category) noexcept:
    # """
    # C callback function for vnaerr
    # """
    self = <NPData>error_arg
    if self._thread_local._vna_data_exception is not None:
        return
    umessage = (<bytes>message).decode("UTF-8")
    if   category == VNAERR_SYSTEM:
        if errno.errorcode == errno.ENOMEM:
            self._thread_local._vna_data_exception = MemoryError(umessage)
        else:
            self._thread_local._vna_data_exception = OSError(
                    errno.errorcode, umessage)
    elif category == VNAERR_USAGE:
        self._thread_local._vna_data_exception = ValueError(umessage)
    elif category == VNAERR_VERSION:
        self._thread_local._vna_data_exception = ValueError(umessage)
    elif category == VNAERR_SYNTAX:
        m = re.match(r'(.*) \(line (\d+)\)', umessage)
        if m:
            details = (m.group(1), m.group(2), 1, "")
        else:
            details = ("", 1, 1, "")
        # TODO: if >= version 3.10, add (m.group(2), 1) to end
        self._thread_local._vna_data_exception = SyntaxError(umessage, details)
    elif category == VNAERR_WARNING:
        if self._thread_local._vna_data_warning is None:
            self._thread_local._vna_data_warning = umessage
    elif category == VNAERR_MATH:
        self._thread_local._vna_data_exception = ArithmeticError(umessage)
    elif category == VNAERR_INTERNAL:
        self._thread_local._vna_data_exception = AssertionError(umessage)
    else:
        self._thread_local._vna_data_exception = Exception(umessage)


def _convert_indices(shape, indices):
    # """
    # Given the shape of an array and an arbitrary index expression,
    # convert to a list of non-zero integer indices or tuples suitable
    # for use as the arguments to the range operator.
    # Args:
    #     shape (tuple): The shape of the array being indexed
    #     indices: The index expression to __getitem__ or __setitem__, which
    #         can be any of: integer, slice, ellipsis, or tuple of these
    # Returns:
    #     tuple(mask, return_list):
    #     - mask (int):
    #         Bitmask showing the locations of the integer indices with lsb
    #         representing the first (leftmost) index
    #     - return_list (list):
    #         List of items: non-negative integer representing an integer index
    #         or a tuple suitable as the arguments for the range operator,
    #         one for each element in shape
    # Raises:
    #     IndexError: The index expression is invalid or an index is
    #         out of bounds.
    # """

    # Bitmask of positions with integer (as opposed to slice) subscripts
    cdef int mask = 0

    # Make sure indices is a tuple.
    if not isinstance(indices, tuple):
        indices = (indices,)

    # Get shape; fail if too many indices in indices.
    cdef int n_dimensions = len(shape)
    cdef int n_indices = len(indices)
    if n_indices > n_dimensions:
        raise ValueError("too many indices for array")

    # For each index, produce the arguments to range.
    ellipsis_found = False
    cdef int i = 0
    cdef int j = 0
    return_list = []
    while i < n_dimensions:
        # If we're out of indices, assume a trailing ellipsis.
        if j >= n_indices:
            return_list.append((0, shape[i], 1))
            i += 1
            continue
        item = indices[j]

        # Handle integer index.  Index must be in bounds.
        if isinstance(item, int):
            if item < 0:
                if item < -shape[i]:
                    raise IndexError(f"index {item} is out of bounds for axis "
                                     f"{i} with size {shape[i]}")
                item += shape[i]

            elif item >= shape[i]:
                raise IndexError(f"index {item} is out of bounds for axis "
                                 f"{i} with size {shape[i]}")

            return_list.append(item)
            mask = mask | (1 << i)
            i += 1
            j += 1
            continue

        # Handle slice index.
        if isinstance(item, slice):
            return_list.append(item.indices(shape[i]))
            i += 1
            j += 1
            continue

        # Handle ellipsis.
        if item is Ellipsis:
            if ellipsis_found:
                raise IndexError("an index can have only a single "
                                 "ellipsis ('...')")
            ellipsis_found = True
            slots_to_add = n_dimensions - (n_indices - 1)
            for k in range(slots_to_add):
                return_list.append((0, shape[i], 1))
                i += 1
            j += 1
            continue

        # Otherwise, error.  Note that we don't currently accept integer
        # or boolean arrays as indices.
        else:
            raise IndexError("indices must be integers or slices only")

    return (mask, return_list)


def _form_double_array(shape, value):
    # """
    # Form value into the expected size array.
    # Args:
    #    shape: list or tuple giving the required dimensions
    #    value: array-like object
    # Returns:
    #    An array in C order with the given dimensions.
    # """
    array = np.asarray(value, dtype=np.double, order="C")
    if array.shape == shape:
        return array

    # Allow broadcasting, e.g. a scalar, over the array.
    array = np.empty(shape, dtype=np.double, order="C")
    array[:] = value
    return array


def _form_complex_array(shape, value):
    # """
    # Form value into the expected size array.
    # Args:
    #     shape: list or tuple giving the required dimensions
    #     value: array-like object
    # Returns:
    #     A complex array in C order with the given dimensions.
    # """
    array = np.asarray(value, dtype=np.complex128, order="C")
    if array.shape == shape:
        return array

    # Allow broadcasting, e.g. a scalar, over the array.
    array = np.empty(shape, dtype=np.complex128, order="C")
    array[:] = value
    return array


cdef class _FrequencyVectorHelper:
    # """
    # Helper class for frequency_vector
    # """
    cdef NPData npd

    def __len__(self):
        # """
        # Return the number of frequencies
        # """
        return self.npd.frequencies

    def __getitem__(self, indices):
        # """
        # Get subscripted item [].
        # """
        cdef vnadata_t *vdp = self.npd.vdp
        cdef int frequencies = vnadata_get_frequencies(vdp)
        cdef int mask
        cdef int findex
        cdef int i

        # Convert indices to a list of tuple arguments suitable for the
        # range operator.  The returned mask gives the positions with
        # integer indices.
        (mask, indices) = _convert_indices((frequencies,), indices)
        assert len(indices) == 1

        # If given an integer index, return a scalar.
        if mask == 1:
            return vnadata_get_frequency(vdp, indices[0])

        # If given a slice, return a vector.
        f_range = range(*indices[0])
        array = np.empty((len(f_range),), dtype=np.double)
        i = 0
        for findex in f_range:
            array[i] = vnadata_get_frequency(vdp, findex)
            i += 1

        return array

    def __setitem__(self, indices, value):
        # """
        # Assign to subscripted item [].
        # """
        cdef vnadata_t *vdp = self.npd.vdp
        cdef int frequencies = vnadata_get_frequencies(vdp)
        cdef int mask
        cdef int findex
        cdef int i
        cdef int rc

        # Convert indices to a list of tuple arguments suitable for the
        # range operator.  The returned mask gives the positions with
        # integer indices.
        (mask, indices) = _convert_indices((frequencies,), indices)
        assert len(indices) == 1

        # If given an integer index, set one entry.
        if mask == 1:
            rc = vnadata_set_frequency(vdp, indices[0], value)
            self.npd._handle_error(rc)

        # Set the sliced items
        else:
            f_range = range(*indices[0])
            array = _form_double_array((len(f_range),), value)
            i = 0
            for findex in f_range:
                rc = vnadata_set_frequency(vdp, findex, array[i])
                self.npd._handle_error(rc)
                i += 1

    def __array__(self, *args, **kwargs):
        # """
        # Return the vector of frequencies.
        # """
        cdef vnadata_t *vdp = self.npd.vdp
        cdef int frequencies = vnadata_get_frequencies(vdp)
        cdef cnp.npy_intp shape[1]
        shape[0] = <cnp.npy_intp>frequencies
        cdef const double *lfp = vnadata_get_frequency_vector(vdp)
        if lfp == NULL:
            self.npd._handle_error(-1)
        a = cnp.PyArray_SimpleNewFromData(
            1, &shape[0], cnp.NPY_DOUBLE, <void *>lfp)
        return np.asarray(a, *args, **kwargs)

    def __copy__(self):
        # """
        # Copy as array.
        # """
        return np.asarray(self)

    def __iter__(self):
        # """
        # Iterate over the parameter matrices.
        # """
        cdef int frequencies = vnadata_get_frequencies(self.npd.vdp)
        cdef vnadata_t *vdp = self.npd.vdp
        for findex in range(frequencies):
            yield vnadata_get_frequency(vdp, findex)

    def __reversed__(self):
        # """
        # Iterate over the parameter matrices in reverse.
        # """
        cdef int frequencies = vnadata_get_frequencies(self.npd.vdp)
        cdef vnadata_t *vdp = self.npd.vdp
        for findex in range(frequencies)[::-1]:
            yield vnadata_get_frequency(vdp, findex)

    def __str__(self):
        # """
        # Show as string.
        # """
        return str(np.asarray(self))

cdef class _DataArrayHelper:
    # """
    # Helper class for data_array.
    # """
    cdef NPData npd

    def _get_matrix(self, int findex):
        # """
        # Return the parameter matrix at given frequency index.
        # Args:
        #     findex (int): frequency index
        # Returns:
        #    copy of matrix
        # """
        cdef vnadata_t *vdp = self.npd.vdp
        cdef int rows = vnadata_get_rows(vdp)
        cdef int columns = vnadata_get_columns(vdp)
        cdef cnp.npy_intp shape[2]
        shape[0] = <cnp.npy_intp>rows
        shape[1] = <cnp.npy_intp>columns
        cdef double complex *clfp
        clfp = vnadata_get_matrix(vdp, findex)
        if clfp == NULL:
            self.npd._handle_error(-1)
        return cnp.PyArray_SimpleNewFromData(
            2, &shape[0], cnp.NPY_COMPLEX128, <void *>clfp)

    def __len__(self):
        # """
        # Return the number of frequencies.
        # """
        return self.npd.frequencies

    def __getitem__(self, indices):
        # """
        # Return subscripted item [].
        # """
        cdef vnadata_t *vdp = self.npd.vdp
        cdef int frequencies = vnadata_get_frequencies(vdp)
        cdef int rows = vnadata_get_rows(vdp)
        cdef int columns = vnadata_get_columns(vdp)
        cdef int mask
        cdef int findex, row, column
        cdef int i, j, k

        # Convert indices to a list of tuple arguments suitable for
        # the range operator.  The returned mask gives the positions
        # with integer indices.
        (mask, indices) = _convert_indices((frequencies, rows, columns),
                                           indices)
        assert len(indices) == 3

        # Handle each case starting with the common case of no slices.
        # Note that at the time of writing, cython doesn't support
        # match-case.
        if mask == IDX_III:
            return vnadata_get_cell(vdp, indices[0], indices[1], indices[2])

        if mask == IDX_IIS:
            c_range = range(*indices[2])
            a = np.empty((len(c_range),), dtype=np.complex128)
            k = 0
            for column in c_range:
                a[k] = vnadata_get_cell(vdp, indices[0], indices[1], column)
                k += 1
            return a

        if mask == IDX_ISI:
            r_range = range(*indices[1])
            a = np.empty((len(r_range),), dtype=np.complex128)
            j = 0
            for row in r_range:
                a[j] = vnadata_get_cell(vdp, indices[0], row, indices[2])
                j += 1
            return a

        if mask == IDX_ISS:
            # If we're given an index of the form [int, :, :], optimize
            # by using _get_matrix.
            findex = indices[0]
            if indices[1] == (0, rows, 1) and indices[2] == (0, columns, 1):
                return self._get_matrix(findex)

            # Else, iterate over the slices.
            r_range = range(*indices[1])
            c_range = range(*indices[2])
            a = np.empty((len(r_range), len(c_range)), dtype=np.complex128)
            j = 0
            for row in r_range:
                k = 0
                for column in c_range:
                    a[j, k] = vnadata_get_cell(vdp, indices[0], row, column)
                    k += 1
                j += 1
            return a

        if mask == IDX_SII:
            f_range = range(*indices[0])
            a = np.empty((len(f_range),), dtype=np.complex128)
            i = 0
            for findex in f_range:
                a[i] = vnadata_get_cell(vdp, findex, indices[1], indices[2])
                i += 1
            return a

        if mask == IDX_SIS:
            f_range = range(*indices[0])
            c_range = range(*indices[2])
            a = np.empty((len(f_range), len(c_range)), dtype=np.complex128)
            i = 0
            for findex in f_range:
                k = 0
                for column in c_range:
                    a[i, k] = vnadata_get_cell(vdp, findex, indices[1],
                                               column)
                    k += 1
                i += 1
            return a

        if mask == IDX_SSI:
            f_range = range(*indices[0])
            r_range = range(*indices[1])
            a = np.empty((len(f_range), len(r_range)), dtype=np.complex128)
            i = 0
            for findex in f_range:
                j = 0
                for row in r_range:
                    a[i, j] = vnadata_get_cell(vdp, findex, row, indices[2])
                    j += 1
                i += 1
            return a

        if mask == IDX_SSS:
            f_range = range(*indices[0])
            r_range = range(*indices[1])
            c_range = range(*indices[2])
            a = np.empty((len(f_range), len(r_range), len(c_range)),
                         dtype=np.complex128)
            i = 0
            for findex in f_range:
                j = 0
                for row in r_range:
                    k = 0
                    for column in c_range:
                        a[i, j, k] = vnadata_get_cell(vdp, findex, row, column)
                        k += 1
                    j += 1
                i += 1
            return a

        assert(False)

    def __setitem__(self, indices, value):
        # """
        # Assign to subscripted item [].
        # """
        # cdef vnadata_t *vdp = self.npd.vdp
        cdef vnadata_t *vdp = self.npd.vdp
        cdef int frequencies = vnadata_get_frequencies(vdp)
        cdef int rows = vnadata_get_rows(vdp)
        cdef int columns = vnadata_get_columns(vdp)
        cdef int mask
        cdef int findex, row, column
        cdef int i, j, k
        cdef int rc
        cdef double complex [:, ::1] v2
        cdef double complex [:, :, ::1] v3

        # Convert indices to a list of tuple arguments suitable for
        # the range operator.  The returned mask gives the positions
        # with integer indices.
        (mask, indices) = _convert_indices((frequencies, rows, columns),
                                           indices)
        assert len(indices) == 3

        # Handle each case starting with the common case of no slices.
        # Note that at the time of writing, cython doesn't support
        # match-case.
        if mask == IDX_III:
            array = _form_complex_array((1,), value)
            rc = vnadata_set_cell(vdp, indices[0], indices[1], indices[2],
                                  array[0])
            self.npd._handle_error(rc)

        elif mask == IDX_IIS:
            c_range = range(*indices[2])
            c_len = len(c_range)
            array = _form_complex_array((c_len,), value)
            k = 0
            for column in c_range:
                rc = vnadata_set_cell(vdp, indices[0], indices[1], column,
                                      array[k])
                self.npd._handle_error(rc)
                k += 1

        elif mask == IDX_ISI:
            r_range = range(*indices[1])
            r_len = len(r_range)
            array = _form_complex_array((r_len,), value)
            j = 0
            for row in r_range:
                rc = vnadata_set_cell(vdp, indices[0], row, indices[2],
                                      array[j])
                self.npd._handle_error(rc)
                j += 1

        elif mask == IDX_ISS:
            # If we're given an index of the form [int, :, :], optimize
            # by using vnadata_set_matrix.
            findex = indices[0]
            if indices[1] == (0, rows, 1) and indices[2] == (0, columns, 1):
                array = _form_complex_array((rows, columns), value)
                v2 = array
                rc = vnadata_set_matrix(vdp, findex, &v2[0][0])
                self.npd._handle_error(rc)

            # Else, iterate over the slices.
            else:
                r_range = range(*indices[1])
                c_range = range(*indices[2])
                r_len = len(r_range)
                c_len = len(c_range)
                array = _form_complex_array((r_len, c_len), value)
                v2 = array
                j = 0
                for row in r_range:
                    k = 0
                    for column in c_range:
                        rc = vnadata_set_cell(vdp, indices[0], row, column,
                                              v2[j, k])
                        self.npd._handle_error(rc)
                        k += 1
                    j += 1

        elif mask == IDX_SII:
            f_range = range(*indices[0])
            f_len = len(f_range)
            array = _form_complex_array((f_len,), value)
            i = 0
            for findex in f_range:
                rc = vnadata_set_cell(vdp, findex, indices[1], indices[2],
                                      array[i])
                self.npd._handle_error(rc)
                i += 1

        elif mask == IDX_SIS:
            f_range = range(*indices[0])
            c_range = range(*indices[2])
            f_len = len(f_range)
            c_len = len(c_range)
            array = _form_complex_array((f_len, c_len), value)
            v2 = array
            i = 0
            for findex in f_range:
                k = 0
                for column in c_range:
                    rc = vnadata_set_cell(
                        vdp, findex, indices[1], column, v2[i, k])
                    self.npd._handle_error(rc)
                    k += 1
                i += 1

        elif mask == IDX_SSI:
            f_range = range(*indices[0])
            r_range = range(*indices[1])
            f_len = len(f_range)
            r_len = len(r_range)
            array = _form_complex_array((f_len, r_len), value)
            v2 = array
            i = 0
            for findex in f_range:
                j = 0
                for row in r_range:
                    rc = vnadata_set_cell(
                        vdp, findex, row, indices[2], v2[i, j])
                    self.npd._handle_error(rc)
                    j += 1
                i += 1

        elif mask == IDX_SSS:
            # If we're given an index of the form [x:y, :, :], optimize
            # by using vnadata_set_matrix.
            f_range = range(*indices[0])
            f_len = len(f_range)
            if indices[1] == (0, rows, 1) and indices[2] == (0, columns, 1):
                array = _form_complex_array((f_len, rows, columns), value)
                v3 = array
                i = 0
                for findex in f_range:
                    rc = vnadata_set_matrix(vdp, findex, &v3[i, 0, 0])
                    i += 1

            # Else, iterate over the slices.
            else:
                r_range = range(*indices[1])
                c_range = range(*indices[2])
                r_len = len(r_range)
                c_len = len(c_range)
                array = _form_complex_array((f_len, r_len, c_len), value)
                v3 = array
                i = 0
                for findex in f_range:
                    j = 0
                    for row in r_range:
                        k = 0
                        for column in c_range:
                            rc = vnadata_set_cell(
                                vdp, findex, row, column, v3[i, j, k])
                            self.npd._handle_error(rc)
                            k += 1
                        j += 1
                    i += 1

        else:
            assert(False)

    def __array__(self, *args, **kwargs):
        # """
        # Return the vector of parameter matrices as a 3D array.
        # """
        cdef vnadata_t *vdp = self.npd.vdp
        cdef int frequencies = vnadata_get_frequencies(vdp)
        cdef int rows = vnadata_get_rows(vdp)
        cdef int columns = vnadata_get_columns(vdp)
        array = np.empty((frequencies, rows, columns), dtype=np.complex128)
        for findex in range(frequencies):
            for row in range(rows):
                for column in range(columns):
                    array[findex, row, column] = \
                            vnadata_get_cell(vdp, findex, row, column)
        return np.asarray(array, *args, **kwargs)

    def __copy__(self):
        # """
        # Copy as array.
        # """
        return np.asarray(self)

    def __iter__(self):
        # """
        # Iterate over the parameter matrices.
        # """
        cdef int frequencies = vnadata_get_frequencies(self.npd.vdp)
        for findex in range(frequencies):
            yield self._get_matrix(findex)

    def __reversed__(self):
        # """
        # Iterate over the parameter matrices in reverse.
        # """
        cdef int frequencies = vnadata_get_frequencies(self.npd.vdp)
        for findex in range(frequencies)[::-1]:
            yield(self._get_matrix(findex))

    def __str__(self):
        # """
        # Show as string.
        # """
        return str(np.asarray(self))

cdef class _Z0VectorHelper:
    # """
    # Helper class for the z0 vector.
    # """
    cdef NPData npd

    def __len__(self):
        # """
        # Return the number of ports.
        # """
        cdef vnadata_t *vdp = self.npd.vdp
        cdef int rows = vnadata_get_rows(vdp)
        cdef int columns = vnadata_get_columns(vdp)
        cdef int ports = max(rows, columns)
        return ports

    def __getitem__(self, indices):
        # """
        # Get subscripted item [].
        # """
        cdef vnadata_t *vdp = self.npd.vdp
        cdef int rows = vnadata_get_rows(vdp)
        cdef int columns = vnadata_get_columns(vdp)
        cdef int ports = max(rows, columns)
        cdef int mask
        cdef int port
        cdef int i

        # Convert indices to a list of tuple arguments suitable for the
        # range operator.  The returned mask gives the positions with
        # integer indices.
        (mask, indices) = _convert_indices((ports,), indices)
        assert len(indices) == 1

        # If given an integer index, return a scalar.
        if mask == 1:
            return vnadata_get_z0(vdp, indices[0])

        # If given a slice, return a vector.
        p_range = range(*indices[0])
        array = np.empty((len(p_range),), dtype=np.complex128)
        i = 0
        for port in p_range:
            array[i] = vnadata_get_z0(vdp, port)
            i += 1

        return array

    def __setitem__(self, indices, value):
        # """
        # Assign to subscripted item [].
        # """
        cdef vnadata_t *vdp = self.npd.vdp
        cdef int rows = vnadata_get_rows(vdp)
        cdef int columns = vnadata_get_columns(vdp)
        cdef int ports = max(rows, columns)
        cdef int mask
        cdef int port
        cdef int i
        cdef int rc

        # Convert indices to a list of tuple arguments suitable for the
        # range operator.  The returned mask gives the positions with
        # integer indices.
        (mask, indices) = _convert_indices((ports,), indices)
        assert len(indices) == 1

        # If given an integer index, set one entry.
        if mask == 1:
            rc = vnadata_set_z0(vdp, indices[0], value)
            self.npd._handle_error(rc)

        # Set the sliced items.
        else:
            p_range = range(*indices[0])
            array = _form_complex_array((len(p_range),), value)
            i = 0
            for port in p_range:
                rc = vnadata_set_z0(vdp, port, array[i])
                self.npd._handle_error(rc)
                i += 1

    def __array__(self, *args, **kwargs):
        # """
        # Return the vector of frequencies.
        # """
        cdef vnadata_t *vdp = self.npd.vdp
        cdef int rows = vnadata_get_rows(vdp)
        cdef int columns = vnadata_get_columns(vdp)
        cdef int ports = max(rows, columns)
        cdef cnp.npy_intp shape[1]
        shape[0] = <cnp.npy_intp>ports
        cdef const double complex *clfp = vnadata_get_z0_vector(vdp)
        if clfp == NULL:
            self.npd._handle_error(-1)
        a = cnp.PyArray_SimpleNewFromData(
            1, &shape[0], cnp.NPY_COMPLEX128, <void *>clfp)
        return np.asarray(a, *args, **kwargs)

    def __copy__(self):
        # """
        # Copy as array.
        # """
        return np.asarray(self)

    def __iter__(self):
        # """
        # Iterate over the parameter matrices.
        # """
        cdef vnadata_t *vdp = self.npd.vdp
        cdef int rows = vnadata_get_rows(vdp)
        cdef int columns = vnadata_get_columns(vdp)
        cdef int ports = max(rows, columns)
        for port in range(ports):
            yield vnadata_get_z0(vdp, port)

    def __reversed__(self):
        # """
        # Iterate over the parameter matrices in reverse.
        # """
        cdef vnadata_t *vdp = self.npd.vdp
        cdef int rows = vnadata_get_rows(vdp)
        cdef int columns = vnadata_get_columns(vdp)
        cdef int ports = max(rows, columns)
        for port in range(ports)[::-1]:
            yield vnadata_get_z0(vdp, port)

    def __str__(self):
        # """
        # Show as string.
        # """
        return str(np.asarray(self))

cdef class _FZ0ArrayHelper:
    # """
    # Helper class for the fz0 vector.
    # """
    cdef NPData npd

    def __len__(self):
        # """
        # Return the number of ports.
        # """
        cdef vnadata_t *vdp = self.npd.vdp
        return vnadata_get_frequencies(vdp)

    def __getitem__(self, indices):
        # """
        # Get subscripted item [].
        # """
        cdef vnadata_t *vdp = self.npd.vdp
        cdef int frequencies = vnadata_get_frequencies(vdp)
        cdef int rows = vnadata_get_rows(vdp)
        cdef int columns = vnadata_get_columns(vdp)
        cdef int ports = max(rows, columns)
        cdef int mask
        cdef int findex, port
        cdef int i, j

        # Convert indices to a list of tuple arguments suitable for the
        # range operator.  The returned mask gives the positions with
        # integer indices.
        (mask, indices) = _convert_indices((frequencies, ports), indices)
        assert len(indices) == 2

        # II
        if mask == 3:
            return vnadata_get_fz0(vdp, indices[0], indices[1])

        # SI
        elif mask == 2:
            f_range = range(*indices[0])
            array = np.empty((len(f_range),), dtype=np.complex128)
            i = 0
            for findex in f_range:
                array[i] = vnadata_get_fz0(vdp, findex, indices[1])
                i += 1
            return array

        # IS
        elif mask == 1:
            p_range = range(*indices[1])
            array = np.empty((len(p_range),), dtype=np.complex128)
            j = 0
            for port in p_range:
                array[j] = vnadata_get_fz0(vdp, indices[0], port)
                j += 1
            return array

        # SS
        elif mask == 0:
            f_range = range(*indices[0])
            p_range = range(*indices[1])
            array = np.empty((len(f_range), len(p_range)),
                             dtype=np.complex128)
            i = 0
            for findex in f_range:
                j = 0
                for port in p_range:
                    array[i, j] = vnadata_get_fz0(vdp, findex, port)
                    j += 1
                i += 1
            return array

        else:
            assert(False)

    def __setitem__(self, indices, value):
        # """
        # Assign to subscripted item [].
        # """
        cdef vnadata_t *vdp = self.npd.vdp
        cdef int frequencies = vnadata_get_frequencies(vdp)
        cdef int rows = vnadata_get_rows(vdp)
        cdef int columns = vnadata_get_columns(vdp)
        cdef int ports = max(rows, columns)
        cdef int mask
        cdef int findex, port
        cdef int i, j
        cdef double complex [::1] v1
        cdef int rc

        # Convert indices to a list of tuple arguments suitable for the
        # range operator.  The returned mask gives the positions with
        # integer indices.
        (mask, indices) = _convert_indices((frequencies, ports), indices)
        assert len(indices) == 2

        # I=integer index, S=slice
        # II
        if mask == 3:
            rc = vnadata_set_fz0(vdp, indices[0], indices[1], value)
            self.npd._handle_error(rc)

        # SI
        elif mask == 2:
            f_range = range(*indices[0])
            port = indices[1]
            array = _form_complex_array((len(f_range),), value)
            i = 0
            for findex in f_range:
                rc = vnadata_set_fz0(vdp, findex, port, array[i])
                self.npd._handle_error(rc)
                i += 1

        # IS
        elif mask == 1:
            findex = indices[0]
            if indices[1] == (0, ports, 1):
                array = _form_complex_array((ports,), value)
                v1 = array
                rc = vnadata_set_fz0_vector(vdp, findex, &v1[0])
                self.npd._handle_error(rc)

            else:
                p_range = range(*indices[1])
                array = _form_complex_array((len(p_range),), value)
                j = 0
                for port in p_range:
                    rc = vnadata_set_fz0(vdp, findex, port, array[j])
                    self.npd._handle_error(rc)
                    j += 1

        # SS
        elif mask == 0:
            f_range = range(*indices[0])
            p_range = range(*indices[1])
            array = _form_complex_array((len(f_range), len(p_range)), value)
            i = 0
            for findex in f_range:
                j = 0
                for port in p_range:
                    rc = vnadata_set_fz0(vdp, findex, port, array[i, j])
                    j += 1
                i += 1

        else:
            assert(False)

    def __array__(self, *args, **kwargs):
        # """
        # Return the fz0 vector.
        # """
        cdef vnadata_t *vdp = self.npd.vdp
        cdef int frequencies = vnadata_get_frequencies(vdp)
        cdef int rows = vnadata_get_rows(vdp)
        cdef int columns = vnadata_get_columns(vdp)
        cdef int ports = max(rows, columns)
        array = np.empty((frequencies, ports), dtype=np.complex128)
        for findex in range(frequencies):
            for port in range(ports):
                array[findex, port] = vnadata_get_fz0(vdp, findex, port)
        return np.asarray(array, *args, **kwargs)

    def __copy__(self):
        # """
        # Copy as array.
        # """
        return np.asarray(self)

    def __iter__(self):
        # """
        # Iterate over the frequencies.
        # """
        cdef vnadata_t *vdp = self.npd.vdp
        cdef int frequencies = vnadata_get_frequencies(vdp)
        for findex in range(frequencies):
            yield self[findex, :]

    def __reversed__(self):
        # """
        # Iterate over the frequencies reverse.
        # """
        cdef vnadata_t *vdp = self.npd.vdp
        cdef int frequencies = vnadata_get_frequencies(vdp)
        for findex in range(frequencies)[::-1]:
            yield self[findex, :]

    def __str__(self):
        # """
        # Show as string.
        # """
        return str(np.asarray(self))


cdef class NPData:
    """
    In-memory representation of electrical network parameter data.

    Args:
        ptype (PType, optional):
            Type of parameter data to be stored: ANY, S, T, U, Z, Y,
            H, G, A, B, ZIN.  The default is ANY, which stores raw data
            without interpretation.

        rows (int, optional), columns (int, optional):
            Dimensions of the parameter matrix.  Rows and columns must
            be equal for S, Z and Y parameters.  They must both be 2
            for T, U, H, G, A and B parameters.  Rows must be 1 for
            ZIN parameters, which are stored as a row vector.  Rows and
            columns can be any non-negative values for ANY parameters.
            Both dimensions default to zero.  Type and dimensions can be
            changed using the :func:`init`, :func:`resize` or :func:`load`
            functions, or by assignment to the :py:attr:`type`,
            :py:attr:`data_array`, or :py:attr:`fz0_array` attributes.

        frequencies (int, optional): Number of frequency points.
            Number of frequencies defaults to zero.  The number can be
            changed using the :func:`init`, :func:`resize`, :func:`load`,
            or :func:`add_frequency` functions, or by assignment to
            the :py:attr:`frequency_vector` or or :py:attr:`fz0_array`
            attributes.

    Raises:
        ValueError: if rows, columns and ptype are not valid

    This class stores a vector of frequency points
    :py:attr:`frequency_vector`, a frequencies by rows by columns array
    of network parameter data :py:attr:`data_array`, and a per-port array
    of complex reference impedances.  The reference impedances can be
    the same across all frequencies (z0_vector), or can be different
    for each frequency (fz0_array).

    When given non-zero dimensions, the class is constructed with all data
    and frequency values initialized to zero, and z0 values initialized
    to 50 ohms.

    The class provides methods to load and save in Touchstone versions
    1 (.s2p) and 2 (.ts), and in a more general space-separated value
    Network Parameter Data (.npd) format.  It also provides conversion
    between network parameter types.
    """
    cdef vnadata_t *vdp
    cdef object _thread_local

    ######################################################################
    # Internal Functions
    ######################################################################

    cdef _handle_error(self, int rc):
        # """
        # Check if _error_fn has saved an exception for this thread or
        # if rc is -1.  In either case, raise an exception or warning.
        #
        # Args:
        #    rc: Return value from C function
        #
        # Raises:
        #    See exceptions in _error_fn.
        # """
        exception = self._thread_local._vna_data_exception
        warning = self._thread_local._vna_data_warning
        self._thread_local._vna_data_exception = None
        self._thread_local._vna_data_warning = None
        if rc == -1 and exception is None:
            PyErr_SetFromErrno(OSError)
        if exception is not None:
            raise exception
        if warning is not None:
            warnings.warn(warning)

    ######################################################################
    # Allocation and Initialization
    ######################################################################

    def __cinit__(self, ptype: PType = PType.ANY,
                  int rows = 0, int columns = 0, int frequencies = 0):
        self._thread_local = local()
        self._thread_local._vna_data_exception = None
        self._thread_local._vna_data_warning = None
        self.vdp = vnadata_alloc_and_init(
            <vnaerr_error_fn_t *>&_error_fn, <void *>self,
            <vnadata_parameter_type_t>ptype, rows, columns, frequencies)
        self._handle_error(0 if self.vdp else -1)

    def __dealloc__(self):
        """
        Free C resources when a libvna.npd.NPData object is garbage collected.
        """
        vnadata_free(self.vdp)

    def init(self, ptype, rows, columns, frequencies):
        """
        Change the ptype and dimensions as indicated by the arguments,
        and reset all frequency and data elements back to zero and all
        z0 entries to 50 ohms.

        Args:
            ptype (PType): Set the parameter type:
                ANY, S, T, U, Z, Y, H, G, A, B, ZIN.

            rows (int):
            columns (int):
                Dimensions of the parameter matrix.  Rows and columns
                must be equal for S, Z and Y parameters.  They must both
                be 2 for T, U, H, G, A and B parameters.  Rows must
                be 1 for ZIN parameters.  Rows and columns can be any
                non-negative values for ANY parameters.

            frequencies (int):
                Number of frequency points

        """
        cdef int rc
        rc = vnadata_init(self.vdp, ptype, rows, columns, frequencies)
        self._handle_error(rc)

    def resize(self, ptype, rows, columns, frequencies):
        """
        Change the parameter type and dimensions without clearing or
        converting data.  Existing values remain undisturbed when the
        matrix type, number of rows or number of frequencies is changed,
        but shift to other cells when the number of columns is changed.
        Changing the parameter type with this function does not convert
        existing data to the new type.  To convert parameters, use
        :func:`convert` instead.

        Args:
            ptype (PType): Set the parameter type:
                ANY, S, T, U, Z, Y, H, G, A, B, ZIN.

            rows (int):
            columns (int):
                Dimensions of the parameter matrix.  Rows and columns
                must be equal for S, Z and Y parameters.  They must both
                be 2 for T, U, H, G, A and B parameters.  Rows must
                be 1 for ZIN parameters.  Rows and columns can be any
                non-negative values for ANY parameters.

            frequencies (int):
                Number of frequency points

        """
        cdef int rc
        rc = vnadata_resize(self.vdp, ptype, rows, columns, frequencies)
        self._handle_error(rc)

    @property
    def ptype(self):
        """
        The parameter type as a PType enum value.  Valid values are:
            ANY, S, T, U, Z, Y, H, G, A, B, ZIN

        When changing this value, the dimensions must be consistent
        with the new type.  Existing data is *not* converted.  See the
        :func:`resize` and :func:`convert` functions for alternative
        ways to change the type.

        Type: PType
        """
        return vnadata_get_type(self.vdp)

    @ptype.setter
    def ptype(self, value):
        # no docstring for setter
        cdef int rc
        rc = vnadata_set_type(self.vdp, <vnadata_parameter_type_t>value)
        self._handle_error(rc)

    @property
    def ptype_name(self):
        """
        The parameter type as a string (read-only).

        Type: str
        """
        return vnadata_get_type_name(vnadata_get_type(self.vdp))

    @property
    def rows(self):
        """
        The number of rows (read-only).

        Type: int
        """
        return vnadata_get_rows(self.vdp)

    @property
    def columns(self):
        """
        The number of columns (read-only).

        Type: int
        """
        return vnadata_get_columns(self.vdp)

    @property
    def frequencies(self):
        """
        The number of frequencies (read-only).

        Type: int
        """
        return vnadata_get_frequencies(self.vdp)

    def _get_vdp(self):
        # """
        # PyCapsule containing the C vnadata_t pointer
        # """
        return PyCapsule_New(<void *>self.vdp, NULL, NULL)

    ######################################################################
    # The Frequency Vector
    ######################################################################

    @property
    def frequency_vector(self):
        """
        The vector of frequency points (read-write).  Slicing operations
        are supported.

        Type: behaves like array[f_index] of float
        """
        cdef _FrequencyVectorHelper helper = _FrequencyVectorHelper()
        helper.npd = self
        return helper

    @frequency_vector.setter
    def frequency_vector(self, vector):
        # no docstring for setter
        array = np.asarray(vector, dtype=np.double)
        if array.ndim != 1:
            raise ValueError("frequency_vector: value must have exactly "
                             "one dimension")

        cdef int rc
        rc = vnadata_resize(self.vdp,
                            vnadata_get_type(self.vdp),
                            vnadata_get_rows(self.vdp),
                            vnadata_get_columns(self.vdp),
                            array.shape[0])
        self._handle_error(rc)

        # Set the vector.
        cdef double [:] cvector = array
        rc = vnadata_set_frequency_vector(self.vdp, &cvector[0])
        self._handle_error(rc)

    def add_frequency(self, double frequency):
        """
        Add a new frequency entry.  The corresponding new data matrix
        (:py:attr:`data_array`\[-1, :, :]) is initialized to zeros.
        This function is useful when loading data into the object
        when the number of frequency points is not known in advance.

        Args:
            frequency (float): new frequency value to add
        """
        cdef int rc
        rc = vnadata_add_frequency(self.vdp, frequency)
        self._handle_error(rc)

    ######################################################################
    # Data Elements
    ######################################################################

    @property
    def data_array(self):
        """
        The parameter data array (read-write).  Slicing operations
        are supported.

        Type: behaves like array [f_index, row, column] of complex
        """
        cdef _DataArrayHelper helper = _DataArrayHelper()
        helper.npd = self
        return helper

    @data_array.setter
    def data_array(self, vector):
        # no docstring for setter
        array = np.asarray(vector, dtype=np.complex128, order="C")
        if array.ndim != 3:
            raise ValueError("data_array: value must have three dimensions")

        cdef int rc
        rc = vnadata_resize(self.vdp,
                            vnadata_get_type(self.vdp),
                            array.shape[1],
                            array.shape[2],
                            array.shape[0])
        self._handle_error(rc)

        # Fill in the values
        cdef int findex
        cdef double complex [:, :, :] c_array = array
        for findex in range(array.shape[0]):
            rc = vnadata_set_matrix(self.vdp, findex,
                                    &c_array[findex][0][0])
            self._handle_error(rc)

    ######################################################################
    # Ordinary System Impedances
    ######################################################################

    @property
    def z0_vector(self):
        """
        The vector of reference impedances for each port (read-write).
        If the dimensions have been established, setting this attribute
        to a scalar value sets all ports to the same reference impedance.
        The default value is 50 ohms.

        Type: behaves like array[port] of complex
        """
        cdef _Z0VectorHelper helper = _Z0VectorHelper()
        helper.npd = self
        return helper

    @z0_vector.setter
    def z0_vector(self, vector):
        # no docstring for setter
        cdef vnadata_t *vdp = self.vdp
        cdef int rows = vnadata_get_rows(vdp)
        cdef int columns = vnadata_get_columns(vdp)
        cdef int ports = max(rows, columns)
        cdef int rc
        if ports == 0:
            # If rows and columns are both zero, allow automatic resize
            # to square matrix matching number of given elements.
            array = np.asarray(vector)
            if array.ndim == 0:
                raise ValueError("z0_vector: cannot broadcast scalar "
                                 "until dimensions have been set")
            if array.ndim != 1:
                raise ValueError("z0_vector: value must have exactly "
                                 "one dimension")

            ports = array.shape[0]
            rc = vnadata_resize(self.vdp,
                                vnadata_get_type(self.vdp),
                                ports, ports, vnadata_get_frequencies(self.vdp))
            self._handle_error(rc)

        else:
            # Otherwise, size must match or be broadcastable
            array = _form_complex_array((ports,), vector)

        # Set the vector.
        cdef double complex [:] cvector = array
        rc = vnadata_set_z0_vector(vdp, &cvector[0])
        self._handle_error(rc)

    ######################################################################
    # Frequency-Dependent System Impedances
    ######################################################################

    @property
    def has_fz0(self):
        """
        True if frequency-dependent system impedances are in effect
        (read-only).

        Type: bool
        """
        return vnadata_has_fz0(self.vdp)

    @property
    def fz0_array(self):
        """
        Frequency-dependent reference impedances as an array
        (read-write).  Slicing operations are supported.  Assigning to
        or modifying elements in this array automatically establishes
        frequency-dependent reference impedances.

        Type: behaves like array[f_index, port] of complex
        """
        cdef _FZ0ArrayHelper helper = _FZ0ArrayHelper()
        helper.npd = self
        return helper

    @fz0_array.setter
    def fz0_array(self, value):
        # no docstring for setter
        cdef vnadata_t *vdp = self.vdp
        cdef int frequencies = vnadata_get_frequencies(vdp)
        cdef int rows = vnadata_get_rows(vdp)
        cdef int columns = vnadata_get_columns(vdp)
        cdef int ports = max(rows, columns)
        cdef int rc

        # If the vnadata object has all zero dimensions, allow
        # automatic resize to a square data matrix based on the
        # value.
        if ports == 0 and frequencies == 0:
            array = np.asarray(value)
            if array.ndim == 0:
                raise ValueError("fz0_array: cannot broadcast scalar "
                                 "until dimensions have been set")
            if array.ndim != 2:
                raise ValueError("fz0_array: expected (frequencies x ports) "
                                 "array")
            frequencies = array.shape[0]
            ports = array.shape[1]
            rc = vnadata_resize(self.vdp,
                                vnadata_get_type(self.vdp),
                                ports, ports, frequencies)
            self._handle_error(rc)

        # Otherwise, it must be broadcastable to the current size.
        else:
            array = _form_complex_array((frequencies, ports), value)

        # Set the vector.
        cdef double complex [:, :] cvector = array
        for findex in range(frequencies):
            rc = vnadata_set_fz0_vector(vdp, findex, &cvector[findex][0])
            self._handle_error(rc)

    ######################################################################
    # Parameter Conversion
    ######################################################################

    def convert(self, PType new_ptype):
        """
        Return a new Data object with data converted to the new type.

        Args:
            new_type (PType): new parameter type:
                ANY, S, T, U, Z, Y, H, G, A, B, ZIN.

        Returns:
            A new NPData object in the requested type.  The original
            object is left unchanged.

        Raises:
            ValueError: if conversion to new_type is invalid
        """
        result = NPData()
        cdef int rc
        rc = vnadata_convert(self.vdp, result.vdp,
                             <vnadata_parameter_type_t>new_ptype)
        self._handle_error(rc)
        return result

    ######################################################################
    # Load and Save
    ######################################################################

    def load(self, filename):
        """
        Load network parameter data from filename, automatically resetting
        the type and dimensions as needed.

        Args:
            filename (str): pathname to file

        Raises:
            OSError:     if can't open file
            SyntaxError: if the file is badly formed
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode("UTF-8")
        cdef const unsigned char[:] cfilename = filename
        cdef int rc
        rc = vnadata_load(self.vdp, <const char *>&cfilename[0])
        self._handle_error(rc)

    def fload(self, filehandle, filename):
        """
        Load network parameter data from an open file handle,
        automatically resetting the type and dimensions as needed.

        Args:
            filehandle: open file handle
            filename:   name of file (used in error messages only)

        Raises:
            OSError:     if can't read file
            SyntaxError: if the file is badly formed
        """
        # TODO enhancement: catch the UnsupportedOperation exception
        #    around the  call to fileno(), e.g. if passed something like
        #    a StringIO object and emulate the behavior by saving the
        #    stream to a temporary file, rewinding and then handing the
        #    open fd to libvna.
        cdef off_t start_position = filehandle.tell()
        cdef int fd1 = filehandle.fileno()
        lseek(fd1, start_position, 0)
        cdef int fd2 = dup(fd1)
        self._handle_error(fd2)
        cdef FILE *fp = fdopen(fd2, "r")
        if fp is NULL:
            close(fd2)
            self._handle_error(-1)
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode("UTF-8")
        cdef const unsigned char[:] cfilename = filename
        cdef int rc
        rc = vnadata_fload(self.vdp, fp, <const char *>&cfilename[0])
        fclose(fp)
        cdef off_t end_position = lseek(fd1, 0, 1)
        filehandle.seek(end_position)
        self._handle_error(rc)

    def save(self, filename):
        """
        Save network parameter data to filename.

        Args:
            filename (str): pathname to file

        Raises:
            OSError:        if can't open file
            ValueError:     if file type inconsistent with parameter data
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode("UTF-8")
        cdef const unsigned char[:] cfilename = filename
        cdef int rc
        rc = vnadata_save(self.vdp, <const char *>&cfilename[0])
        self._handle_error(rc)

    def fsave(self, filehandle, filename):
        """
        Save network parameter data to an open file handle.

        Args:
            filehandle: open file handle to write
            filename:   name of file used in error messages only

        Raises:
            OSError:        if can't write file
            ValueError:     if file type inconsistent with parameter data
        """
        # TODO enhancement: catch the UnsupportedOperation exception
        #    around the  call to fileno(), e.g. if passed something like
        #    a StringIO object and emulate the behavior by having libvna
        #    save to a temporary file then rewind and read back the file
        #    and write each line to the given stream.
        cdef off_t start_position = filehandle.tell()
        cdef int fd1 = filehandle.fileno()
        lseek(fd1, start_position, 0)
        cdef int fd2 = dup(fd1)
        self._handle_error(fd2)
        cdef FILE *fp = fdopen(fd2, "w")
        if fp is NULL:
            close(fd2)
            self._handle_error(-1)
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode("UTF-8")
        cdef const unsigned char[:] cfilename = filename
        cdef int rc
        rc = vnadata_fsave(self.vdp, fp, <const char *>&cfilename[0])
        fclose(fp)
        cdef off_t end_position = lseek(fd1, 0, 1)
        filehandle.seek(end_position)
        self._handle_error(rc)

    def cksave(self, filename):
        """
        Test if we would be able to save the currently selected format
        in the currently selected filetype without actually saving.

        This function is useful to test if there's a conflict betweent the
        requested parameter format and the save file type before doing
        potentially expensive operations such as VNA measurements only
        to ultimately fail at save time.  The *filename* argument is used
        only to determine the filetype when :py:attr:`filetype` is AUTO.

        Args:
            filename (str): proposed pathname to save file

        Raises:
            ValueError: if file type inconsistent with parameter data
        """
        if isinstance(filename, unicode):
            filename = (<unicode>filename).encode("UTF-8")
        cdef const unsigned char[:] cfilename = filename
        cdef int rc
        rc = vnadata_cksave(self.vdp, <const char *>&cfilename[0])
        self._handle_error(rc)

    @property
    def filetype(self):
        """
        The file type (read-write):

        - AUTO:        Determine type based on filename extension (default)
        - TOUCHSTONE1: Use Touchstone version 1 (.s2p)
        - TOUCHSTONE2: Use Touchstone version 2 (.ts)
        - NPD:         Use network parameter data (.ndp)

        Type: FileType

        When a file is loaded, *filetype* is set to the type of the
        loaded file.  To save in a different format with type based on
        the file extension, set *filetype* back to AUTO before saving.
        """
        cdef int rc = vnadata_get_filetype(self.vdp)
        self._handle_error(rc)
        return <FileType>rc

    @filetype.setter
    def filetype(self, value):
        # no docstring for setter
        cdef rc = vnadata_set_filetype(self.vdp, value)
        self._handle_error(rc)

    @property
    def format(self):
        """
        The file format (str, read-write):

        Valid values:

        - S[ri|ma|dB]:  scattering parameters
        - T[ri|ma|dB]:  scattering-transfer parameters
        - U[ri|ma|dB]:  inverse scattering-transfer parameters
        - Z[ri|ma]:     impedance parameters
        - Y[ri|ma]:     admittance parameters
        - H[ri|ma]:     hybrid parameters
        - G[ri|ma]:     inverse-hybrid parameters
        - A[ri|ma]:     ABCD parameters
        - B[ri|ma]:     inverse ABCD parameters
        - Zin[ri|ma]:   impedance looking into each port
        - PRC:          Zin as parallel resistance and capacitance
        - PRL:          Zin as parallel resistance and inductance
        - SRC:          Zin as series resistance and capacitance
        - SRL:          Zin as series resistance and inducatance
        - IL:           insertion loss (dB)
        - RL:           return loss (dB)
        - VSWR:         voltage standing wave ratio

        where the ri,  ma  or  dB  suffix  is  an  optional  coordinate
        system modifier:

        - ri:  real, imaginary
        - ma:  magnitude, angle
        - dB:  decibels, angle

        Format specifiers are case-insensitive.  Note that not all
        formats can be represented in all file formats.  Touchstone
        formats support only S, Y, Z, H, and G parameters, with version
        1 imposing further restrictions on the maximum number of ports
        (4), and disallowing per-port Z0 values.  The .npd file type
        can store all formats.  With the .npd type only, format may be
        a comma-separated list of formats, e.g. "IL,RL,VSWR".

        Type: str
        """
        cdef const char *fmt = vnadata_get_format(self.vdp)
        if fmt == NULL:
            self._handle_error(-1)
        return (<bytes>fmt).decode("UTF-8")

    @format.setter
    def format(self, value):
        # no docstring for setter
        if isinstance(value, unicode):
            value = (<unicode>value).encode("UTF-8")
        cdef int rc = vnadata_set_format(self.vdp, <const char *>value)
        self._handle_error(rc)

    @property
    def fprecision(self):
        """
        The number of significant figures to use for frequency values
        when saving parameters to a file.  Default is 7.

        Type: int
        """
        cdef int rc = vnadata_get_fprecision(self.vdp)
        self._handle_error(rc)
        return rc

    @fprecision.setter
    def fprecision(self, value):
        # no docstring for setter
        cdef int rc = vnadata_set_fprecision(self.vdp, value)
        self._handle_error(rc)

    @property
    def dprecision(self):
        """
        The number of significant figures to use for data values when
        saving parameters to a file.  Default is 6.

        Type: int
        """
        cdef int rc = vnadata_get_dprecision(self.vdp)
        self._handle_error(rc)
        return rc

    @dprecision.setter
    def dprecision(self, digits):
        # no docstring for setter
        cdef int rc = vnadata_set_dprecision(self.vdp, digits)
        self._handle_error(rc)

    ######################################################################
    # Python Specific Methods
    ######################################################################

    def __len__(self):
        """
        Return the number of frequencies.
        """
        return vnadata_get_frequencies(self.vdp)

    def __copy__(self):
        """
        Make a copy of the NPData object.
        """
        result = NPData()
        rc = vnadata_convert(self.vdp, result.vdp,
                             <vnadata_parameter_type_t>self.ptype)
        self._handle_error(rc)
        return result

    def __iter__(self):
        """
        Return frequency, parameter matrix tuples.
        """
        if self.has_fz0:
            for findex in range(self.frequencies):
                yield((self.frequencies[findex], self._get_matrix(findex)),
                      self.fz0_array[findex])

        else:
            for findex in range(self.frequencies):
                yield((self.frequencies[findex], self._get_matrix(findex)))

    def __reversed__(self):
        """
        Return frequency, parameter matrix tuples in reverse.
        """
        if self.has_fz0:
            for findex in range(self.frequencies)[::-1]:
                yield((self.frequencies[findex], self._get_matrix(findex)),
                      self.fz0_array[findex])

        else:
            for findex in range(self.frequencies)[::-1]:
                yield((self.frequencies[findex], self._get_matrix(findex)))
