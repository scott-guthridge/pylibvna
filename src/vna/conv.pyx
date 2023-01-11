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

# cython: language_level=3
import  numpy as np

#
# fn_2x2: pointer to function taking two double complex [2][2] matrices
#
ctypedef void (*fn_2x2)(const double complex (*in_2x2)[2],
                        double complex (*out_2x2)[2]) nogil

#
# fn_2x2_with_z0: pointer to function taking two double complex [2][2]
#                 matrices and a vector of z0 values
#
ctypedef void (*fn_2x2_with_z0)(const double complex (*in_2x2)[2],
                                double complex (*out_2x2)[2],
                                const double complex *z0) nogil

#
# fn_NxN: pointer to function taking the addresses of the first elements
#       of two double complex matrices
#
ctypedef void (*fn_NxN)(const double complex *in_2x2,
                        double complex *out_2x2, int) nogil

#
# fn_NxN_with_z0: pointer to function taking the addresses of the first
#     elements of two double complex matrices and a double complex vector
#     of z0
#
ctypedef void (*fn_NxN_with_z0)(const double complex *in_NxN,
                                double complex *out_NxN,
                                const double complex *z0, int n) nogil

#
# fn_2x2_to_zin: pointer to function taking a double complex [2][2] matrix,
#             a length 2 vector of z0 values, and a length 2 output vector
#
ctypedef void (*fn_2x2_to_zin)(const double complex (*in_2x2)[2],
                               double complex *zi,
                               const double complex *z0) nogil

#
# fn_NxN_to_zin: pointer to function taking the address of the first
#               element of an NxN doulbe complex matrix, a length N
#               vector of z0 values, and a length
#
ctypedef void (*fn_NxN_to_zin)(const double complex *in_NxN,
                               double complex *zi, const double complex *z0,
                               int n) nogil

#
# canonicalize_z0: return z0 as [m, n] where m is number of outer elements
#
cdef canonicalize_z0(z0, name, outer_shape, n):
    z0 = np.asarray(z0)
    #
    # scalar
    #
    if z0.shape == tuple():
        return np.asarray([[z0 for i in range(n)]],
                          dtype=np.complex128, order="C")

    #
    # array [1]
    #
    if z0.shape == (1,):
        return np.asarray([[z0[0] for i in range(n)]],
                          dtype=np.complex128, order="C")

    #
    # array [n]
    #
    if z0.shape == (n,):
        return np.asarray([z0],
                          dtype=np.complex128, order="C")

    #
    # array[m][1]
    #
    if z0.shape == outer_shape + (n, 1):
        z0 = np.asarray(z0, dtype=np.complex128, order="C")
        z0.shape = (-1, 1)
        return np.asarray([z0[i, 0]
                           for i in range(z0.shape[0]) for j in range(n)],
                          dtype=np.complex128, order="C")

    #
    # array[m][n]
    #
    if z0.shape == outer_shape + (n,):
        z0 = np.asarray(z0, dtype=np.complex128, order="C")
        z0.shape = (-1, n)
        return z0

    #
    # none of the above
    #
    raise ValueError(name + ": invalid dimensions for z0")

#
# helper_2x2: call fn on each 2x2 array in input
#
cdef helper_2x2(fn_2x2 fn, name, array):
    #
    # Convert the input array to C compatible double complex form,
    # validate that the last two dimensions are 2x2, and create the
    # output array with same shape as input.
    #
    input = np.asarray(array, dtype=np.complex128, order="C")
    if input.ndim < 2 or input.shape[-2] != 2 or input.shape[-1] != 2:
        raise ValueError(name + ": expected 2x2 array")
    output = np.empty(input.shape, dtype=np.complex128, order="C")

    #
    # Create views into the input and output arrays that flatten
    # the outer dimensions.
    #
    vin = input.view()
    vin.shape = (-1, 2, 2)
    vout = output.view()
    vout.shape = (-1, 2, 2)
    cdef double complex [:, :, :] cvin = vin
    cdef double complex [:, :, :] cvout = vout

    #
    # Call fn on each 2x2 array...
    #
    cdef int i
    for i in range(vin.shape[0]):
        fn(<const double complex (*)[2]>&cvin[i, 0, 0],
           <double complex (*)[2]>&cvout[i, 0, 0])

    #
    # If we added an extraneous dimension, remove it.
    #
    return output

#
# helper_2x2_with_z0: call fn on each 2x2 array in input
#
cdef helper_2x2_with_z0(fn_2x2_with_z0 fn, name, array, z0):
    #
    # Convert the input array to C compatible double complex form,
    # validate that the last two dimensions are 2x2, and create the
    # output array with same shape as input.
    #
    input = np.asarray(array, dtype=np.complex128, order="C")
    if input.ndim < 2 or input.shape[-2] != 2 or input.shape[-1] != 2:
        raise ValueError(name + ": expected 2x2 array")
    output = np.empty(input.shape, dtype=np.complex128, order="C")

    #
    # Make views into input and output with outer dimensions flattened
    #
    vin = input.view()
    vin.shape = (-1, 2, 2)
    vout = output.view()
    vout.shape = (-1, 2, 2)
    cdef double complex [:, :, :] cvin = vin
    cdef double complex [:, :, :] cvout = vout

    #
    # Convert z0 to C compatible complex m x 2, where m is the
    # number of outer elements.
    #
    z0 = canonicalize_z0(z0, name, input.shape[:-2], 2)
    cdef double complex [:, :] cz0 = z0

    #
    # Call fn on each 2x2 array...
    #
    cdef int i
    if z0.shape[0] == 1:
        for i in range(vin.shape[0]):
            fn(<const double complex (*)[2]>&cvin[i, 0, 0],
               <double complex (*)[2]>&cvout[i, 0, 0],
               &cz0[0][0])

    else:
        for i in range(vin.shape[0]):
            fn(<const double complex (*)[2]>&cvin[i, 0, 0],
               <double complex (*)[2]>&cvout[i, 0, 0],
               &cz0[i][0])

    return output

#
# helper_NxN: call fn on each NxN array in input
#
cdef helper_NxN(fn_NxN fn, name, array):
    #
    # Convert the input array to C compatible double complex form,
    # validate that the last two dimensions are square, and create
    # the output array with same shape as input.
    #
    input = np.asarray(array, dtype=np.complex128, order="C")
    if input.ndim < 2 or input.shape[-2] != input.shape[-1]:
        raise ValueError(name + ": expected NxN array")
    output = np.empty(input.shape, dtype=np.complex128, order="C")
    n = input.shape[-1]

    #
    # Make views into input and output with outer dimensions flattened
    #
    vin = input.view()
    vin.shape = (-1, n, n)
    vout = output.view()
    vout.shape = (-1, n, n)
    cdef double complex [:, :, :] cvin = vin
    cdef double complex [:, :, :] cvout = vout

    #
    # Call fn on each 2x2 array...
    #
    cdef int i
    for i in range(vin.shape[0]):
        fn(&cvin[i, 0, 0], &cvout[i, 0, 0], n)

    #
    # If we added an extraneous dimension, remove it.
    #
    return output

#
# helper_NxN_with_z0: call fn on each NxN array in input
#
cdef helper_NxN_with_z0(fn_NxN_with_z0 fn, name, array, z0):
    #
    # Convert the input array to C compatible double complex form,
    # validate that the last two dimensions are square, and create
    # the output array with same shape as input.
    #
    input = np.asarray(array, dtype=np.complex128, order="C")
    if input.ndim < 2 or input.shape[-2] != input.shape[-1]:
        raise ValueError(name + ": expected NxN array")
    output = np.empty(input.shape, dtype=np.complex128, order="C")
    n = input.shape[-1]

    #
    # Make views into input and output with outer dimensions flattened
    #
    vin = input.view()
    vin.shape = (-1, n, n)
    vout = output.view()
    vout.shape = (-1, n, n)
    cdef double complex [:, :, :] cvin = vin
    cdef double complex [:, :, :] cvout = vout

    #
    # Convert z0 to C compatible complex m x n, where m is the
    # number of outer elements.
    #
    z0 = canonicalize_z0(z0, name, input.shape[:-2], n)
    cdef double complex [:, :] cz0 = z0

    #
    # Call fn on each 2x2 array...
    #
    cdef int i
    if z0.shape[0] == 1:
        for i in range(vin.shape[0]):
            fn(&cvin[i, 0, 0], &cvout[i, 0, 0], &cz0[0][0], n)

    else:
        for i in range(vin.shape[0]):
            fn(&cvin[i, 0, 0], &cvout[i, 0, 0], &cz0[i][0], n)

    return output

#
# helper_2x2_to_zin: call fn on each 2x2 array in input
#
cdef helper_2x2_to_zin(fn_2x2_to_zin fn, name, array, z0):
    #
    # Convert the input array to C compatible double complex form,
    # validate that the last two dimensions are 2x2, and create the
    # output vector with same shape as input with last component removed.
    #
    input = np.asarray(array, dtype=np.complex128, order="C")
    if input.ndim < 2 or input.shape[-2] != 2 or input.shape[-1] != 2:
        raise ValueError(name + ": expected 2x2 array")
    output = np.empty(input.shape[0:-1], dtype=np.complex128, order="C")

    #
    # Make views into input and output with outer dimensions flattened
    #
    vin = input.view()
    vin.shape = (-1, 2, 2)
    vout = output.view()
    vout.shape = (-1, 2)
    cdef double complex [:, :, :] cvin = vin
    cdef double complex [:, :]    cvout = vout

    #
    # Convert z0 to C compatible complex m x 2, where m is the
    # number of outer elements.
    #
    z0 = canonicalize_z0(z0, name, input.shape[:-2], 2)
    cdef double complex [:, :] cz0 = z0

    #
    # Call fn on each 2x2 array...
    #
    cdef int i
    if z0.shape[0] == 1:
        for i in range(vin.shape[0]):
            fn(<const double complex (*)[2]>&cvin[i, 0, 0],
               &cvout[i, 0], &cz0[0][0])

    else:
        for i in range(vin.shape[0]):
            fn(<const double complex (*)[2]>&cvin[i, 0, 0],
               &cvout[i, 0], &cz0[i][0])

    return output

#
# helper_NxN_to_zin: call fn on each NxN array in input
#
cdef helper_NxN_to_zin(fn_NxN_with_z0 fn, name, array, z0):
    #
    # Convert the input array to C compatible double complex form,
    # validate that the last two dimensions are square, and create
    # the output array with same shape as input with last dimension
    # removed.
    #
    input = np.asarray(array, dtype=np.complex128, order="C")
    if input.ndim < 2 or input.shape[-2] != input.shape[-1]:
        raise ValueError(name + ": expected NxN array")
    output = np.empty(input.shape[0:-1], dtype=np.complex128, order="C")
    n = input.shape[-1]

    #
    # Make views into input and output with outer dimensions flattened
    #
    vin = input.view()
    vin.shape = (-1, n, n)
    vout = output.view()
    vout.shape = (-1, n)
    cdef double complex [:, :, :] cvin = vin
    cdef double complex [:, :]    cvout = vout

    #
    # Convert z0 to C compatible complex m x n, where m is the
    # number of outer elements.
    #
    z0 = canonicalize_z0(z0, name, input.shape[:-2], n)
    cdef double complex [:, :] cz0 = z0

    #
    # Call fn on each 2x2 array...
    #
    cdef int i
    if z0.shape[0] == 1:
        for i in range(vin.shape[0]):
            fn(&cvin[i, 0, 0], &cvout[i, 0], &cz0[0][0], n)

    else:
        for i in range(vin.shape[0]):
            fn(&cvin[i, 0, 0], &cvout[i, 0], &cz0[i][0], n)

    return output


###############################################################################
# 2x2 Conversions Without z0
###############################################################################

def atob(array):
    """
    atob(array):
        Convert ABCD parameters to inverse ABCD parameters.
    """
    return helper_2x2(&vnaconv_atob, "atob", array)


def atog(array):
    """
    atog(array):
        Convert ABCD parameters to inverse hybrid parameters.
    """
    return helper_2x2(&vnaconv_atog, "atog", array)


def atoh(array):
    """
    atoh(array):
        Convert ABCD parameters to hybrid parameters.
    """
    return helper_2x2(&vnaconv_atoh, "atoh", array)


def atoy(array):
    """
    atoy(array):
        Convert ABCD parameters to admittance parameters.
    """
    return helper_2x2(&vnaconv_atoy, "atoy", array)


def atoz(array):
    """
    atoz(array):
        Convert ABCD parameters to impedance parameters.
    """
    return helper_2x2(&vnaconv_atoz, "atoz", array)


def btoa(array):
    """
    btoa(array):
        Convert inverse ABCD parameters to ABCD parameters.
    """
    return helper_2x2(&vnaconv_btoa, "btoa", array)


def btog(array):
    """
    btog(array):
        Convert inverse ABCD parameters to inverse hybrid parameters.
    """
    return helper_2x2(&vnaconv_btog, "btog", array)


def btoh(array):
    """
    btoh(array):
        Convert inverse ABCD parameters to hybrid parameters.
    """
    return helper_2x2(&vnaconv_btoh, "btoh", array)


def btoy(array):
    """
    btoy(array):
        Convert inverse ABCD parameters to admittance parameters.
    """
    return helper_2x2(&vnaconv_btoy, "btoy", array)


def btoz(array):
    """
    btoz(array):
        Convert inverse ABCD parameters to impedance parameters.
    """
    return helper_2x2(&vnaconv_btoz, "btoz", array)


def gtoa(array):
    """
    gtoa(array):
        Convert inverse hybrid parameters to ABCD parameters.
    """
    return helper_2x2(&vnaconv_gtoa, "gtoa", array)


def gtob(array):
    """
    gtob(array):
        Convert inverse hybrid parameters to inverse ABCD parameters.
    """
    return helper_2x2(&vnaconv_gtob, "gtob", array)


def gtoh(array):
    """
    gtoh(array):
        Convert inverse hybrid parameters to hybrid parameters.
    """
    return helper_2x2(&vnaconv_gtoh, "gtoh", array)


def gtoy(array):
    """
    gtoy(array):
        Convert inverse hybrid parameters to inverse hybrid parameters.
    """
    return helper_2x2(&vnaconv_gtoy, "gtoy", array)


def gtoz(array):
    """
    gtoz(array):
        Convert inverse hybrid parameters to impedance parameters.
    """
    return helper_2x2(&vnaconv_gtoz, "gtoz", array)


def htoa(array):
    """
    htoa(array):
        Convert hybrid parameters to ABCD parameters.
    """
    return helper_2x2(&vnaconv_htoa, "htoa", array)


def htob(array):
    """
    htob(array):
        Convert hybrid parameters to inverse ABCD parameters.
    """
    return helper_2x2(&vnaconv_htob, "htob", array)


def htog(array):
    """
    htog(array):
        Convert hybrid parameters to inverse hybrid parameters.
    """
    return helper_2x2(&vnaconv_htog, "htog", array)


def htoy(array):
    """
    htoy(array):
        Convert hybrid parameters to admittance parameters.
    """
    return helper_2x2(&vnaconv_htoy, "htoy", array)


def htoz(array):
    """
    htoz(array):
        Convert hybrid parameters to impedance parameters.
    """
    return helper_2x2(&vnaconv_htoz, "htoz", array)


def stot(array):
    """
    stot(array):
        Convert hybrid parameters to scattering-transfer parameters.
    """
    return helper_2x2(&vnaconv_stot, "stot", array)


def stou(array):
    """
    stou(array):
        Convert hybrid parameters to inverse scattering-transfer parameters.
    """
    return helper_2x2(&vnaconv_stou, "stou", array)


def ttos(array):
    """
    ttos(array):
        Convert scattering-transfer parameters to scattering parameters.
    """
    return helper_2x2(&vnaconv_ttos, "ttos", array)


def ttou(array):
    """
    ttou(array):
        Convert scattering-transfer parameters to inverse
        scattering-transfer parameters.
    """
    return helper_2x2(&vnaconv_ttou, "ttou", array)


def utos(array):
    """
    utos(array):
        Convert inverse scattering-transfer parameters to scattering
        parameters.
    """
    return helper_2x2(&vnaconv_utos, "utos", array)


def utot(array):
    """
    utot(array):
        Convert inverse scattering-transfer parameters to
        scattering-transfer parameters.
    """
    return helper_2x2(&vnaconv_utot, "utot", array)


def ytoa(array):
    """
    ytoa(array):
        Convert admittance parameters to ABCD parameters.
    """
    return helper_2x2(&vnaconv_ytoa, "ytoa", array)


def ytob(array):
    """
    ytob(array):
        Convert admittance parameters to inverse ABCD parameters.
    """
    return helper_2x2(&vnaconv_ytob, "ytob", array)


def ytog(array):
    """
    ytog(array):
        Convert admittance parameters to inverse hybrid parameters.
    """
    return helper_2x2(&vnaconv_ytog, "ytog", array)


def ytoh(array):
    """
    ytoh(array):
        Convert admittance parameters to hybrid parameters.
    """
    return helper_2x2(&vnaconv_ytoh, "ytoh", array)


def ytoz(array):
    """
    ytoz(array):
        Convert admittance parameters to impedance parameters.
    """
    return helper_2x2(&vnaconv_ytoz, "ytoz", array)


def ztoa(array):
    """
    ztoa(array):
        Convert impedance parameters to ABCD parameters.
    """
    return helper_2x2(&vnaconv_ztoa, "ztoa", array)


def ztob(array):
    """
    ztob(array):
        Convert impedance parameters to inverse ABCD parameters.
    """
    return helper_2x2(&vnaconv_ztob, "ztob", array)


def ztog(array):
    """
    ztog(array):
        Convert impedance parameters to inverse hybrid parameters..
    """
    return helper_2x2(&vnaconv_ztog, "ztog", array)


def ztoh(array):
    """
    ztoh(array):
        Convert impedance parameters to hybrid parameters..
    """
    return helper_2x2(&vnaconv_ztoh, "ztoh", array)


def ztoy(array):
    """
    ztoy(array):
        Convert impedance parameters to admittance parameters.
    """
    return helper_2x2(&vnaconv_ztoy, "ztoy", array)


###############################################################################
# 2x2 Conversions With z0
###############################################################################

def atos(array, z0=50.0):
    """
    atos(array, z0=50.0):
        Convert ABCD parameters to scattering parameters.
    """
    return helper_2x2_with_z0(&vnaconv_atos, "atos", array, z0)


def atot(array, z0=50.0):
    """
    atot(array, z0=50.0):
        Convert ABCD parameters to scattering-transfer parameters.
    """
    return helper_2x2_with_z0(&vnaconv_atot, "atot", array, z0)


def atou(array, z0=50.0):
    """
    atou(array, z0=50.0):
        Convert ABCD parameters to inverse scattering-transfer parameters.
    """
    return helper_2x2_with_z0(&vnaconv_atou, "atou", array, z0)


def btos(array, z0=50.0):
    """
    btos(array, z0=50.0):
        Convert inverse ABCD parameters to scattering parameters.
    """
    return helper_2x2_with_z0(&vnaconv_btos, "btos", array, z0)


def btot(array, z0=50.0):
    """
    btot(array, z0=50.0):
        Convert inverse ABCD parameters to scattering-transfer parameters.
    """
    return helper_2x2_with_z0(&vnaconv_btot, "btot", array, z0)


def btou(array, z0=50.0):
    """
    btou(array, z0=50.0):
        Convert inverse ABCD parameters to inverse scattering-transfer
        parameters.
    """
    return helper_2x2_with_z0(&vnaconv_btou, "btou", array, z0)


def gtos(array, z0=50.0):
    """
    gtos(array, z0=50.0):
        Convert inverse hybrid parameters to scattering parameters.
    """
    return helper_2x2_with_z0(&vnaconv_gtos, "gtos", array, z0)


def gtot(array, z0=50.0):
    """
    gtot(array, z0=50.0):
        Convert inverse hybrid parameters to scattering-transfer
        parameters.
    """
    return helper_2x2_with_z0(&vnaconv_gtot, "gtot", array, z0)


def gtou(array, z0=50.0):
    """
    gtou(array, z0=50.0):
        Convert inverse hybrid parameters to inverse scattering-transfer
        parameters.
    """
    return helper_2x2_with_z0(&vnaconv_gtou, "gtou", array, z0)


def htos(array, z0=50.0):
    """
    htos(array, z0=50.0):
        Convert hybrid parameters to scattering parameters.
    """
    return helper_2x2_with_z0(&vnaconv_htos, "htos", array, z0)


def htot(array, z0=50.0):
    """
    htot(array, z0=50.0):
        Convert hybrid parameters to scattering-transfer parameters.
    """
    return helper_2x2_with_z0(&vnaconv_htot, "htot", array, z0)


def htou(array, z0=50.0):
    """
    htou(array, z0=50.0):
        Convert hybrid parameters to inverse scattering-transfer
        parameters.
    """
    return helper_2x2_with_z0(&vnaconv_htou, "htou", array, z0)


def stoa(array, z0=50.0):
    """
    stoa(array, z0=50.0):
        Convert scattering parameters to ABCD parameters.
    """
    return helper_2x2_with_z0(&vnaconv_stoa, "stoa", array, z0)


def stob(array, z0=50.0):
    """
    stob(array, z0=50.0):
        Convert scattering parameters to inverse ABCD parameters.
    """
    return helper_2x2_with_z0(&vnaconv_stob, "stob", array, z0)


def stog(array, z0=50.0):
    """
    stog(array, z0=50.0):
        Convert scattering parameters to inverse hybrid parameters.
    """
    return helper_2x2_with_z0(&vnaconv_stog, "stog", array, z0)


def stoh(array, z0=50.0):
    """
    stoh(array, z0=50.0):
        Convert scattering parameters to hybrid parameters.
    """
    return helper_2x2_with_z0(&vnaconv_stoh, "stoh", array, z0)


def stoy(array, z0=50.0):
    """
    stoy(array, z0=50.0):
        Convert scattering parameters to admittance parameters.
    """
    return helper_2x2_with_z0(&vnaconv_stoy, "stoy", array, z0)


def stoz(array, z0=50.0):
    """
    stoz(array, z0=50.0):
        Convert scattering parameters to impedance parameters.
    """
    return helper_2x2_with_z0(&vnaconv_stoz, "stoz", array, z0)


def ttoa(array, z0=50.0):
    """
    ttoa(array, z0=50.0):
        Convert scattering-transfer parameters to ABCD parameters.
    """
    return helper_2x2_with_z0(&vnaconv_ttoa, "ttoa", array, z0)


def ttob(array, z0=50.0):
    """
    ttob(array, z0=50.0):
        Convert scattering-transfer parameters to inverse ABCD parameters.
    """
    return helper_2x2_with_z0(&vnaconv_ttob, "ttob", array, z0)


def ttog(array, z0=50.0):
    """
    ttog(array, z0=50.0):
        Convert scattering-transfer parameters to inverse hybrid
        parameters.
    """
    return helper_2x2_with_z0(&vnaconv_ttog, "ttog", array, z0)


def ttoh(array, z0=50.0):
    """
    ttoh(array, z0=50.0):
        Convert scattering-transfer parameters to hybrid parameters.
    """
    return helper_2x2_with_z0(&vnaconv_ttoh, "ttoh", array, z0)


def ttoy(array, z0=50.0):
    """
    ttoy(array, z0=50.0):
        Convert scattering-transfer parameters to admittance parameters.
    """
    return helper_2x2_with_z0(&vnaconv_ttoy, "ttoy", array, z0)


def ttoz(array, z0=50.0):
    """
    ttoz(array, z0=50.0):
        Convert scattering-transfer parameters to impedance parameters.
    """
    return helper_2x2_with_z0(&vnaconv_ttoz, "ttoz", array, z0)


def utoa(array, z0=50.0):
    """
    utoa(array, z0=50.0):
        Convert inverse scattering-transfer parameters to ABCD parameters.
    """
    return helper_2x2_with_z0(&vnaconv_utoa, "utoa", array, z0)


def utob(array, z0=50.0):
    """
    utob(array, z0=50.0):
        Convert inverse scattering-transfer parameters to inverse ABCD
        parameters.
    """
    return helper_2x2_with_z0(&vnaconv_utob, "utob", array, z0)


def utog(array, z0=50.0):
    """
    utog(array, z0=50.0):
        Convert inverse scattering-transfer parameters to inverse hybrid
        parameters.
    """
    return helper_2x2_with_z0(&vnaconv_utog, "utog", array, z0)


def utoh(array, z0=50.0):
    """
    utoh(array, z0=50.0):
        Convert inverse scattering-transfer parameters to hybrid
        parameters.
    """
    return helper_2x2_with_z0(&vnaconv_utoh, "utoh", array, z0)


def utoy(array, z0=50.0):
    """
    utoy(array, z0=50.0):
        Convert inverse scattering-transfer parameters to admittance
        parameters.
    """
    return helper_2x2_with_z0(&vnaconv_utoy, "utoy", array, z0)


def utoz(array, z0=50.0):
    """
    utoz(array, z0=50.0):
        Convert inverse scattering-transfer parameters to impedance
        parameters.
    """
    return helper_2x2_with_z0(&vnaconv_utoz, "utoz", array, z0)


def ytos(array, z0=50.0):
    """
    ytos(array, z0=50.0):
        Convert admittance parameters to scattering parameters.
    """
    return helper_2x2_with_z0(&vnaconv_ytos, "ytos", array, z0)


def ytot(array, z0=50.0):
    """
    ytot(array, z0=50.0):
        Convert admittance parameters to scattering-transfer parameters.
    """
    return helper_2x2_with_z0(&vnaconv_ytot, "ytot", array, z0)


def ytou(array, z0=50.0):
    """
    ytou(array, z0=50.0):
        Convert admittance parameters to inverse scattering-transfer
        parameters.
    """
    return helper_2x2_with_z0(&vnaconv_ytou, "ytou", array, z0)


def ztos(array, z0=50.0):
    """
    ztos(array, z0=50.0):
        Convert impedance parameters to scattering parameters.
    """
    return helper_2x2_with_z0(&vnaconv_ztos, "ztos", array, z0)


def ztot(array, z0=50.0):
    """
    ztot(array, z0=50.0):
        Convert impedance parameters to scattering-transfer parameters.
    """
    return helper_2x2_with_z0(&vnaconv_ztot, "ztot", array, z0)


def ztou(array, z0=50.0):
    """
    ztou(array, z0=50.0):
        Convert impedance parameters to inverse scattering-transfer
        parameters.
    """
    return helper_2x2_with_z0(&vnaconv_ztou, "ztou", array, z0)


###############################################################################
# NxN Conversions Without z0
###############################################################################

def ytozn(array):
    """
    ytozn(array):
        Convert NxN admittance parameters to impedance parameters.
    """
    return helper_NxN(&vnaconv_ytozn, "ytozn", array)


def ztoyn(array):
    """
    ztoyn(array):
        Convert NxN impedance parameters to admittance parameters.
    """
    return helper_NxN(&vnaconv_ztoyn, "ztoyn", array)


###############################################################################
# NxN Conversions With z0
###############################################################################

def stoyn(array, z0=50.0):
    """
    stoyn(array, z0=50.0):
        Convert NxN scattering parameters to admittance parameters.
    """
    return helper_NxN_with_z0(&vnaconv_stoyn, "stoyn", array, z0)


def stozn(array, z0=50.0):
    """
    stozn(array, z0=50.0):
        Convert NxN scattering parameters to impedance parameters.
    """
    return helper_NxN_with_z0(&vnaconv_stozn, "stozn", array, z0)


def ytosn(array, z0=50.0):
    """
    ytosn(array, z0=50.0):
        Convert NxN admittance parameters to scattering parameters.
    """
    return helper_NxN_with_z0(&vnaconv_ytosn, "ytosn", array, z0)


def ztosn(array, z0=50.0):
    """
    ztosn(array, z0=50.0):
        Convert NxN impedance parameters to scattering parameters.
    """
    return helper_NxN_with_z0(&vnaconv_ztosn, "ztosn", array, z0)


###############################################################################
# 2x2 Conversions To Zin
###############################################################################

def atozi(array, z0=50.0):
    """
    atozi(array, z0=50.0):
        Convert ABCD parameters to impedances into each port.
    """
    return helper_2x2_to_zin(&vnaconv_atozi, "atozi", array, z0)


def btozi(array, z0=50.0):
    """
    btozi(array, z0=50.0):
        Convert inverse ABCD parameters to impedances into each port.
    """
    return helper_2x2_to_zin(&vnaconv_btozi, "btozi", array, z0)


def gtozi(array, z0=50.0):
    """
    gtozi(array, z0=50.0):
        Convert inverse hybrid parameters to impedances into each port.
    """
    return helper_2x2_to_zin(&vnaconv_gtozi, "gtozi", array, z0)


def htozi(array, z0=50.0):
    """
    htozi(array, z0=50.0):
        Convert hybrid parameters to impedances into each port.
    """
    return helper_2x2_to_zin(&vnaconv_htozi, "htozi", array, z0)


def stozi(array, z0=50.0):
    """
    stozi(array, z0=50.0):
        Convert scattering parameters to impedances into each port.
    """
    return helper_2x2_to_zin(&vnaconv_stozi, "stozi", array, z0)


def ttozi(array, z0=50.0):
    """
    ttozi(array, z0=50.0):
        Convert scattering-transfer parameters to impedances into
        each port.
    """
    return helper_2x2_to_zin(&vnaconv_ttozi, "ttozi", array, z0)


def utozi(array, z0=50.0):
    """
    utozi(array, z0=50.0):
        Convert inverse scattering-transfer parameters to impedances
        into each port.
    """
    return helper_2x2_to_zin(&vnaconv_utozi, "utozi", array, z0)


def ytozi(array, z0=50.0):
    """
    ytozi(array, z0=50.0):
        Convert admittance parameters to impedances into each port.
    """
    return helper_2x2_to_zin(&vnaconv_ytozi, "ytozi", array, z0)


def ztozi(array, z0=50.0):
    """
    ztozi(array, z0=50.0):
        Convert impedance parameters to impedances into each port.
    """
    return helper_2x2_to_zin(&vnaconv_ztozi, "ztozi", array, z0)


###############################################################################
# NxN Conversions To Zin
###############################################################################

def stozin(array, z0=50.0):
    """
    stozin(array, z0=50.0):
        Convert NxN scattering parameters to impedances into each port.
    """
    return helper_NxN_to_zin(&vnaconv_stozin, "stozin", array, z0)


def ytozin(array, z0=50.0):
    """
    ytozin(array, z0=50.0):
        Convert NxN admittance parameters to impedances into each port.
    """
    return helper_NxN_to_zin(&vnaconv_ytozin, "ytozin", array, z0)


def ztozin(array, z0=50.0):
    """
    ztozin(array, z0=50.0):
        Convert NxN impedance parameters to impedances into each port.
    """
    return helper_NxN_to_zin(&vnaconv_ztozin, "ztozin", array, z0)
