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
from libc.stdint cimport uint32_t
from libc.stdio cimport FILE
cimport libvna.data as data

cdef extern from "<stdbool.h>":
    ctypedef bint bool


cdef extern from "<stdarg.h>":
    ctypedef struct va_list:
        pass    # opaque


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


cdef extern from "<vnaproperty.h>":
    ctypedef struct vnaproperty_t:
        pass    # opaque

    int vnaproperty_vtype(const vnaproperty_t *root,
                          const char *format, va_list ap)
    int vnaproperty_vcount(const vnaproperty_t *root,
                           const char *format, va_list ap)
    const char **vnaproperty_vkeys(const vnaproperty_t *root,
                                   const char *format, va_list ap)
    const char *vnaproperty_vget(const vnaproperty_t *root,
                                 const char *format, va_list ap)
    int vnaproperty_vset(vnaproperty_t **rootptr,
                         const char *format, va_list ap)
    int vnaproperty_vdelete(vnaproperty_t **rootptr,
                            const char *format, va_list ap)
    vnaproperty_t *vnaproperty_vget_subtree(const vnaproperty_t *root,
                                            const char *format, va_list ap)
    vnaproperty_t **vnaproperty_vset_subtree(vnaproperty_t **rootptr,
                                             const char *format, va_list ap)
    int vnaproperty_type(const vnaproperty_t *root, const char *format, ...)
    int vnaproperty_count(const vnaproperty_t *root, const char *format, ...)
    const char **vnaproperty_keys(const vnaproperty_t *root,
                                  const char *format, ...)
    const char *vnaproperty_get(const vnaproperty_t *root,
                                const char *format, ...)
    int vnaproperty_set(vnaproperty_t **rootptr, const char *format, ...)
    int vnaproperty_delete(vnaproperty_t **rootptr, const char *format, ...)
    vnaproperty_t *vnaproperty_get_subtree(const vnaproperty_t *root,
                                           const char *format, ...)
    vnaproperty_t **vnaproperty_set_subtree(vnaproperty_t **rootptr,
                                            const char *format, ...)
    int vnaproperty_copy(vnaproperty_t **destination,
                         const vnaproperty_t *source);
    char *vnaproperty_quote_key(const char *key);
    int vnaproperty_import_yaml_from_string(vnaproperty_t **rootptr,
                                            const char *input,
                                            vnaerr_error_fn_t *error_fn,
                                            void *error_arg);
    int vnaproperty_import_yaml_from_file(vnaproperty_t **rootptr, FILE *fp,
                                          const char *filename,
                                          vnaerr_error_fn_t *error_fn,
                                          void *error_arg);
    int vnaproperty_export_yaml_to_file(const vnaproperty_t *root, FILE *fp,
                                        const char *filename,
                                        vnaerr_error_fn_t *error_fn,
                                        void *error_arg);


cdef extern from "<vnadata.h>":
    ctypedef struct vnadata_t:
        pass    # opaque

cdef extern from "<vnacal.h>":
    ctypedef enum vnacal_calkit_type_t:
        VNACAL_CALKIT_SHORT   = 0x636b0000
        VNACAL_CALKIT_OPEN    = 0x636b0001
        VNACAL_CALKIT_LOAD    = 0x636b0002
        VNACAL_CALKIT_THROUGH = 0x636b0003

    enum:
        VNACAL_CKF_TRADITIONAL = 0x0001

    ctypedef struct vnacal_calkit_data_t:
        vnacal_calkit_type_t vcd_type
        uint32_t vcd_flags
        double vcd_offset_delay
        double vcd_offset_loss
        double vcd_offset_z0
        double vcd_fmin
        double vcd_fmax
        # Note: these are really in a union, but Cython doesn't make
        # any assumptions about memory layout in these structures.
        double vcd_l_coefficients[4]
        double vcd_c_coefficients[4]
        double complex vcd_zl

    ctypedef enum vnacal_type_t:
        VNACAL_NOTYPE       = -1
        VNACAL_T8           =  0
        VNACAL_U8           =  1
        VNACAL_TE10         =  2
        VNACAL_UE10         =  3
        VNACAL_T16          =  4
        VNACAL_U16          =  5
        VNACAL_UE14         =  6
        VNACAL_E12          =  7

    enum:
        VNACAL_MATCH        =  0
        VNACAL_OPEN         =  1
        VNACAL_SHORT        = -1
        VNACAL_ZERO         =  0
        VNACAL_ONE          =  1

    ctypedef struct vnacal_t:
        pass    # opaque

    ctypedef struct vnacal_new_t:
        pass    # opaque

    ctypedef enum vnacal_z0_type_t:
        VNACAL_Z0_SCALAR = 0
        VNACAL_Z0_VECTOR = 1
        VNACAL_Z0_MATRIX = 2

    vnacal_type_t vnacal_name_to_type(const char *name)
    const char *vnacal_type_to_name(vnacal_type_t type)
    vnacal_new_t *vnacal_new_alloc(vnacal_t *vcp, vnacal_type_t type,
                                   int rows, int columns, int frequencies)
    int vnacal_new_set_frequency_vector(vnacal_new_t *vnp,
                                        const double *frequency_vector)
    int vnacal_new_set_z0(vnacal_new_t *vnp, double complex z0)
    int vnacal_new_set_z0_vector(vnacal_new_t *vnp,
        const double complex *z0_vector, int length)
    int vnacal_new_set_m_error(vnacal_new_t *vnp,
                               const double *frequency_vector,
                               int frequencies,
                               const double *sigma_nf_vector,
                               const double *sigma_tr_vector)
    int vnacal_new_set_p_tolerance(vnacal_new_t *vnp, double tolerance)
    int vnacal_new_set_et_tolerance(vnacal_new_t *vnp, double tolerance)
    int vnacal_new_set_iteration_limit(vnacal_new_t *vnp, int iterations)
    int vnacal_new_set_pvalue_limit(vnacal_new_t *vnp, double significance)
    int vnacal_new_add_single_reflect(vnacal_new_t *vnp,
                                      double complex *const *a,
                                      int a_rows, int a_columns,
                                      double complex *const *b,
                                      int b_rows, int b_columns,
                                      int s11, int port)
    int vnacal_new_add_single_reflect_m(vnacal_new_t *vnp,
                                        double complex *const *m,
                                        int m_rows, int m_columns,
                                        int s11, int port)
    int vnacal_new_add_double_reflect(vnacal_new_t *vnp,
                                      double complex *const *a,
                                      int a_rows, int a_columns,
                                      double complex *const *b,
                                      int b_rows, int b_columns,
                                      int s11, int s22,
                                      int port1, int port2)
    int vnacal_new_add_double_reflect_m(vnacal_new_t *vnp,
                                        double complex *const *m,
                                        int m_rows, int m_columns,
                                        int s11, int s22,
                                        int port1, int port2)
    int vnacal_new_add_line(vnacal_new_t *vnp,
                            double complex *const *a,
                            int a_rows, int a_columns,
                            double complex *const *b,
                            int b_rows, int b_columns,
                            const int *s_2x2, int port1, int port2)
    int vnacal_new_add_line_m(vnacal_new_t *vnp,
                              double complex *const *m,
                              int m_rows, int m_columns,
                              const int *s_2x2,
                              int port1, int port2)
    int vnacal_new_add_through(vnacal_new_t *vnp,
                               double complex *const *a,
                               int a_rows, int a_columns,
                               double complex *const *b,
                               int b_rows, int b_columns,
                               int port1, int port2)
    int vnacal_new_add_through_m(vnacal_new_t *vnp,
                                 double complex *const *m,
                                 int m_rows, int m_columns,
                                 int port1, int port2)
    int vnacal_new_add_mapped_matrix(vnacal_new_t *vnp,
                                     double complex *const *a,
                                     int a_rows, int a_columns,
                                     double complex *const *b,
                                     int b_rows, int b_columns,
                                     const int *s,
                                     int s_rows, int s_columns,
                                     const int *port_map)
    int vnacal_new_add_mapped_matrix_m(vnacal_new_t *vnp,
                                       double complex *const *m,
                                       int m_rows, int m_columns,
                                       const int *s,
                                       int s_rows, int s_columns,
                                       const int *port_map)
    int vnacal_new_solve(vnacal_new_t *vnp)
    void vnacal_new_free(vnacal_new_t *vnp)
    int vnacal_make_scalar_parameter(vnacal_t *vcp, double complex coefficient)
    int vnacal_make_vector_parameter(vnacal_t *vcp,
                                     const double *frequency_vector,
                                     int frequencies,
                                     const double complex *coefficient_vector)
    int vnacal_make_calkit_parameter(vnacal_t *vcp,
                                     const vnacal_calkit_data_t *vcdp)
    int vnacal_make_calkit_parameter_matrix(vnacal_t *vcp,
                                            const vnacal_calkit_data_t *vcdp,
                                            int *parameter_matrix,
                                            size_t parameter_matrix_size)
    int vnacal_make_data_parameter(vnacal_t *vcp, const vnadata_t *vdp)
    int vnacal_make_data_parameter_matrix(vnacal_t *vcp,
                                          const vnadata_t *vcdp,
                                          int *parameter_matrix,
                                          size_t parameter_matrix_size)
    int vnacal_load_data_parameter_matrix(vnacal_t *vcp,
        const char *filename, int *parameter_matrix,
        size_t parameter_matrix_size)
    int vnacal_embed_parameter(vnacal_t *vcp, int dut,
            const int (*fixture_matrix)[2])
    int vnacal_embed_parameter_matrix(vnacal_t *vcp,
        const int *dut_matrix, int dut_ports,
        const int *fixture_matrix,
        int *result_matrix, size_t result_matrix_size)
    int vnacal_deembed_parameter(vnacal_t *vcp, int embedded,
        const int (*fixture_matrix)[2])
    int vnacal_deembed_parameter_matrix(vnacal_t *vcp,
        const int *embedded_matrix, int embedded_ports,
        const int *fixture_matrix,
        int *result_matrix, size_t result_matrix_size)
    int vnacal_embed(vnacal_t *vcp, vnadata_t *vdp_in, vnadata_t *vdp_out,
        const int *fixture_matrix, int fixture_ports)
    int vnacal_deembed(vnacal_t *vcp, vnadata_t *vdp_in, vnadata_t *vdp_out,
        const int *fixture_matrix, int fixture_ports)
    int vnacal_make_unknown_parameter(vnacal_t *vcp, int initial_guess)
    int vnacal_make_correlated_parameter(vnacal_t *vcp, int other,
                                         const double *sigma_frequency_vector,
                                         int sigma_frequencies,
                                         const double *sigma_vector)
    double complex vnacal_get_parameter_value(vnacal_t *vcp,
                                              int parameter,
                                              double frequency)
    double complex vnacal_eval_parameter(vnacal_t *vcp, int parameter,
                                         double frequency, double complex z0)
    int vnacal_eval_parameter_matrix(vnacal_t *vcp,
                                     const int *parameter_matrix,
                                     int rows, int columns, double frequency,
                                     const double complex *z0_vector,
                                     double complex *result_matrix)
    int vnacal_parameter_matrix_to_data(vnacal_t *vcp,
                                        const int *parameter_matrix,
                                        int rows, int columns, vnadata_t *vdp)
    int vnacal_delete_parameter(vnacal_t *vcp, int parameter)
    int vnacal_delete_parameter_matrix(vnacal_t *vcp,
                                       const int *parameter_matrix,
                                       int rows, int columns)
    vnacal_t *vnacal_create(vnaerr_error_fn_t *error_fn, void *error_arg)
    vnacal_t *vnacal_load(const char *pathname,
                          vnaerr_error_fn_t *error_fn, void *error_arg)
    int vnacal_save(vnacal_t *vcp, const char *pathname)
    const char *vnacal_get_filename(const vnacal_t *vcp)
    int vnacal_add_calibration(vnacal_t *vcp, const char *name,
                               vnacal_new_t *vnp)
    int vnacal_find_calibration(const vnacal_t *vcp, const char *name)
    int vnacal_delete_calibration(vnacal_t *vcp, int ci)
    int vnacal_get_calibration_end(const vnacal_t *vcp)
    const char *vnacal_get_name(const vnacal_t *vcp, int ci)
    vnacal_type_t vnacal_get_type(const vnacal_t *vcp, int ci)
    int vnacal_get_rows(const vnacal_t *vcp, int ci)
    int vnacal_get_columns(const vnacal_t *vcp, int ci)
    int vnacal_get_frequencies(const vnacal_t *vcp, int ci)
    double vnacal_get_fmin(const vnacal_t *vcp, int ci)
    double vnacal_get_fmax(const vnacal_t *vcp, int ci)
    const double *vnacal_get_frequency_vector(const vnacal_t *vcp, int ci)
    vnacal_z0_type_t vnacal_get_z0_type(const vnacal_t *vcp, int ci)
    double complex vnacal_get_z0(const vnacal_t *vcp, int ci)
    int vnacal_get_z0_vector(const vnacal_t *vcp, int ci, double f,
                             double complex *vector, int max_entries)
    int vnacal_set_fprecision(vnacal_t *vcp, int precision)
    int vnacal_set_dprecision(vnacal_t *vcp, int precision)
    int vnacal_property_type(vnacal_t *vcp, int ci, const char *format, ...)
    int vnacal_property_count(vnacal_t *vcp, int ci, const char *format, ...)
    const char **vnacal_property_keys(vnacal_t *vcp, int ci,
                                      const char *format, ...)
    const char *vnacal_property_get(vnacal_t *vcp, int ci,
                                    const char *format, ...)
    int vnacal_property_set(vnacal_t *vcp, int ci, const char *format, ...)
    int vnacal_property_delete(vnacal_t *vcp, int ci, const char *format, ...)
    vnaproperty_t *vnacal_property_get_subtree(vnacal_t *vcp, int ci,
                                               const char *format, ...)
    vnaproperty_t **vnacal_property_set_subtree(vnacal_t *vcp, int ci,
                                                const char *format, ...)
    void vnacal_free(vnacal_t *vcp)
    int vnacal_apply(vnacal_t *vcp, int ci,
                     const double *frequency_vector, int frequencies,
                     double complex *const *a, int a_rows, int a_columns,
                     double complex *const *b, int b_rows, int b_columns,
                     data.vnadata_t *s_parameters)
    int vnacal_apply_m(vnacal_t *vcp, int ci,
                       const double *frequency_vector, int frequencies,
                       double complex *const *m, int m_rows, int m_columns,
                       data.vnadata_t *s_parameters)
