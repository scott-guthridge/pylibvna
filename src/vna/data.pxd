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

from libc.stdio cimport FILE

cdef extern from "<stdbool.h>":
    ctypedef bint bool

cdef extern from "<vnaerr.h>":
    ctypedef enum vnaerr_category_t:
        VNAERR_SYSTEM      = 0
        VNAERR_USAGE       = 1
        VNAERR_VERSION     = 2
        VNAERR_SYNTAX      = 3
        VNAERR_WARNING     = 4
        VNAERR_MATH        = 5
        VNAERR_INTERNAL    = 6

    ctypedef void vnaerr_error_fn_t(const char *message, void *error_arg,
                                    vnaerr_category_t category)

cdef extern from "<vnadata.h>":
    ctypedef enum vnadata_parameter_type_t:
        pass    # see ParameterType

    ctypedef enum vnadata_filetype_t:
        pass    # see FileType

    ctypedef struct vnadata_t:
        pass    # opaque

    vnadata_t *vnadata_alloc(vnaerr_error_fn_t *error_fn, void *error_arg)
    void vnadata_free(vnadata_t *vdp)
    int vnadata_init(vnadata_t *vdp, vnadata_parameter_type_t ptype,
                     int rows, int columns, int frequencies)
    void _vnadata_bounds_error(const char *function, const vnadata_t *vdp,
                               const char *what, int value)
    vnadata_t *vnadata_alloc_and_init(vnaerr_error_fn_t *error_fn,
                                      void *error_arg,
                                      vnadata_parameter_type_t ptype,
                                      int rows, int columns, int frequencies)
    int vnadata_resize(vnadata_t *vdp, vnadata_parameter_type_t ptype,
                       int rows, int columns, int frequencies)
    int vnadata_get_frequencies(const vnadata_t *vdp)
    int vnadata_get_rows(const vnadata_t *vdp)
    int vnadata_get_columns(const vnadata_t *vdp)
    vnadata_parameter_type_t vnadata_get_type(const vnadata_t *vdp)
    int vnadata_set_type(vnadata_t *vdp, vnadata_parameter_type_t type)
    double vnadata_get_fmin(const vnadata_t *vdp)
    double vnadata_get_fmax(const vnadata_t *vdp)
    double vnadata_get_frequency(const vnadata_t *vdp, int findex)
    int vnadata_set_frequency(vnadata_t *vdp, int findex, double frequency)
    const double *vnadata_get_frequency_vector(const vnadata_t *vdp)
    int vnadata_set_frequency_vector(vnadata_t *vdp,
                                     const double *frequency_vector)
    double complex vnadata_get_cell(const vnadata_t *vdp,
                                    int findex, int row, int column)
    int vnadata_set_cell(vnadata_t *vdp, int findex, int row,
                         int column, double complex value)
    double complex *vnadata_get_matrix(const vnadata_t *vdp,
                                       int findex)
    int vnadata_set_matrix(vnadata_t *vdp, int findex,
                           const double complex *matrix)
    int vnadata_get_to_vector(const vnadata_t *vdp,
                              int row, int column, double complex *vector)
    int vnadata_set_from_vector(vnadata_t *vdp, int row, int column,
                                const double complex *vector)
    double complex vnadata_get_z0(const vnadata_t *vdp, int port)
    int vnadata_set_z0(vnadata_t *vdp, int port, double complex z0)
    int vnadata_set_all_z0(vnadata_t *vdp, double complex z0)
    const double complex *vnadata_get_z0_vector(const vnadata_t *vdp)
    int vnadata_set_z0_vector(vnadata_t *vdp, const double complex *z0_vector)
    bool vnadata_has_fz0(const vnadata_t *vdp)
    complex vnadata_get_fz0(const vnadata_t *vdp, int findex, int port)
    int vnadata_set_fz0(vnadata_t *vdp, int findex, int port,
                        double complex z0)
    const double complex *vnadata_get_fz0_vector(const vnadata_t *vdp,
                                                 int findex)
    int vnadata_set_fz0_vector(vnadata_t *vdp, int findex,
                               const double complex *z0_vector)
    int vnadata_convert(const vnadata_t *vdp_in,
                        vnadata_t *vdp_out,
                        vnadata_parameter_type_t new_parameter)
    int vnadata_add_frequency(vnadata_t *vdp, double frequency)
    const char *vnadata_get_type_name(vnadata_parameter_type_t type)
    vnadata_filetype_t vnadata_get_filetype(const vnadata_t *vdp)
    int vnadata_set_filetype(vnadata_t *vdp, vnadata_filetype_t filetype)
    const char *vnadata_get_format(const vnadata_t *vdp)
    int vnadata_set_format(vnadata_t *vdp, const char *format)
    int vnadata_get_fprecision(const vnadata_t *vdp)
    int vnadata_set_fprecision(vnadata_t *vdp, int precision)
    int vnadata_get_dprecision(const vnadata_t *vdp)
    int vnadata_set_dprecision(vnadata_t *vdp, int precision)
    int vnadata_load(vnadata_t *vdp, const char *filename)
    int vnadata_fload(vnadata_t *vdp, FILE *fp, const char *filename)
    int vnadata_cksave(vnadata_t *vdp, const char *filename)
    int vnadata_save(vnadata_t *vdp, const char *filename)
    int vnadata_fsave(vnadata_t *vdp, FILE *fp, const char *filename)
