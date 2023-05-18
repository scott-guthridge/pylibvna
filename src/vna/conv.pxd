#cython: language_level=3
#cython: binding=True
#
# Python Bindings for Vector Network Analyzer Library
# Copyright © 2022 D Scott Guthridge <scott_guthridge@rompromity.net>
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

cdef extern from "<vnaconv.h>":
    #
    # 2x2 Conversions
    #
    void vnaconv_atob(const double complex (*a)[2], double complex (*b)[2])
    void vnaconv_atog(const double complex (*a)[2], double complex (*g)[2])
    void vnaconv_atoh(const double complex (*a)[2], double complex (*h)[2])
    void vnaconv_atos(const double complex (*a)[2], double complex (*s)[2],
                      const double complex *z0)
    void vnaconv_atot(const double complex (*a)[2], double complex (*t)[2],
                      const double complex *z0)
    void vnaconv_atou(const double complex (*a)[2], double complex (*u)[2],
                      const double complex *z0)
    void vnaconv_atoy(const double complex (*a)[2], double complex (*y)[2])
    void vnaconv_atoz(const double complex (*a)[2], double complex (*z)[2])
    void vnaconv_atozi(const double complex (*a)[2], double complex *zi,
                       const double complex *z0)
    void vnaconv_btoa(const double complex (*b)[2], double complex (*a)[2])
    void vnaconv_btog(const double complex (*b)[2], double complex (*g)[2])
    void vnaconv_btoh(const double complex (*b)[2], double complex (*h)[2])
    void vnaconv_btos(const double complex (*b)[2], double complex (*s)[2],
                      const double complex *z0)
    void vnaconv_btot(const double complex (*b)[2], double complex (*t)[2],
                      const double complex *z0)
    void vnaconv_btou(const double complex (*b)[2], double complex (*u)[2],
                      const double complex *z0)
    void vnaconv_btoy(const double complex (*b)[2], double complex (*y)[2])
    void vnaconv_btoz(const double complex (*b)[2], double complex (*z)[2])
    void vnaconv_btozi(const double complex (*b)[2], double complex *zi,
                       const double complex *z0)
    void vnaconv_gtoa(const double complex (*g)[2], double complex (*a)[2])
    void vnaconv_gtob(const double complex (*g)[2], double complex (*b)[2])
    void vnaconv_gtoh(const double complex (*g)[2], double complex (*h)[2])
    void vnaconv_gtos(const double complex (*g)[2], double complex (*s)[2],
                      const double complex *z0)
    void vnaconv_gtot(const double complex (*g)[2], double complex (*t)[2],
                      const double complex *z0)
    void vnaconv_gtou(const double complex (*g)[2], double complex (*u)[2],
                      const double complex *z0)
    void vnaconv_gtoy(const double complex (*g)[2], double complex (*y)[2])
    void vnaconv_gtoz(const double complex (*g)[2], double complex (*z)[2])
    void vnaconv_gtozi(const double complex (*g)[2], double complex *zi,
                       const double complex *z0)
    void vnaconv_htoa(const double complex (*h)[2], double complex (*a)[2])
    void vnaconv_htob(const double complex (*h)[2], double complex (*b)[2])
    void vnaconv_htog(const double complex (*h)[2], double complex (*g)[2])
    void vnaconv_htos(const double complex (*h)[2], double complex (*s)[2],
                      const double complex *z0)
    void vnaconv_htot(const double complex (*h)[2], double complex (*t)[2],
                      const double complex *z0)
    void vnaconv_htou(const double complex (*h)[2], double complex (*u)[2],
                      const double complex *z0)
    void vnaconv_htoy(const double complex (*h)[2], double complex (*y)[2])
    void vnaconv_htoz(const double complex (*h)[2], double complex (*z)[2])
    void vnaconv_htozi(const double complex (*h)[2], double complex *zi,
                       const double complex *z0)
    void vnaconv_stoa(const double complex (*s)[2], double complex (*a)[2],
                      const double complex *z0)
    void vnaconv_stob(const double complex (*s)[2], double complex (*b)[2],
                      const double complex *z0)
    void vnaconv_stog(const double complex (*s)[2], double complex (*g)[2],
                      const double complex *z0)
    void vnaconv_stoh(const double complex (*s)[2], double complex (*h)[2],
                      const double complex *z0)
    void vnaconv_stot(const double complex (*s)[2], double complex (*t)[2])
    void vnaconv_stou(const double complex (*s)[2], double complex (*u)[2])
    void vnaconv_stoy(const double complex (*s)[2], double complex (*y)[2],
                      const double complex *z0)
    void vnaconv_stoz(const double complex (*s)[2], double complex (*z)[2],
                      const double complex *z0)
    void vnaconv_stozi(const double complex (*s)[2], double complex *zi,
                       const double complex *z0)
    void vnaconv_ttoa(const double complex (*t)[2], double complex (*a)[2],
                      const double complex *z0)
    void vnaconv_ttob(const double complex (*t)[2], double complex (*b)[2],
                      const double complex *z0)
    void vnaconv_ttog(const double complex (*t)[2], double complex (*g)[2],
                      const double complex *z0)
    void vnaconv_ttoh(const double complex (*t)[2], double complex (*h)[2],
                      const double complex *z0)
    void vnaconv_ttos(const double complex (*t)[2], double complex (*s)[2])
    void vnaconv_ttou(const double complex (*t)[2], double complex (*u)[2])
    void vnaconv_ttoy(const double complex (*t)[2], double complex (*y)[2],
                      const double complex *z0)
    void vnaconv_ttoz(const double complex (*t)[2], double complex (*z)[2],
                      const double complex *z0)
    void vnaconv_ttozi(const double complex (*t)[2], double complex *zi,
                       const double complex *z0)
    void vnaconv_utoa(const double complex (*u)[2], double complex (*a)[2],
                      const double complex *z0)
    void vnaconv_utob(const double complex (*u)[2], double complex (*b)[2],
                      const double complex *z0)
    void vnaconv_utog(const double complex (*u)[2], double complex (*g)[2],
                      const double complex *z0)
    void vnaconv_utoh(const double complex (*u)[2], double complex (*h)[2],
                      const double complex *z0)
    void vnaconv_utos(const double complex (*u)[2], double complex (*s)[2])
    void vnaconv_utot(const double complex (*u)[2], double complex (*t)[2])
    void vnaconv_utoy(const double complex (*u)[2], double complex (*y)[2],
                      const double complex *z0)
    void vnaconv_utoz(const double complex (*u)[2], double complex (*z)[2],
                      const double complex *z0)
    void vnaconv_utozi(const double complex (*u)[2], double complex *zi,
                       const double complex *z0)
    void vnaconv_ytoa(const double complex (*y)[2], double complex (*a)[2])
    void vnaconv_ytob(const double complex (*y)[2], double complex (*b)[2])
    void vnaconv_ytog(const double complex (*y)[2], double complex (*g)[2])
    void vnaconv_ytoh(const double complex (*y)[2], double complex (*h)[2])
    void vnaconv_ytos(const double complex (*y)[2], double complex (*s)[2],
                      const double complex *z0)
    void vnaconv_ytot(const double complex (*y)[2], double complex (*t)[2],
                      const double complex *z0)
    void vnaconv_ytou(const double complex (*y)[2], double complex (*u)[2],
                      const double complex *z0)
    void vnaconv_ytoz(const double complex (*y)[2], double complex (*z)[2])
    void vnaconv_ytozi(const double complex (*y)[2], double complex *zi,
                       const double complex *z0)
    void vnaconv_ztoa(const double complex (*z)[2], double complex (*a)[2])
    void vnaconv_ztob(const double complex (*z)[2], double complex (*b)[2])
    void vnaconv_ztog(const double complex (*z)[2], double complex (*g)[2])
    void vnaconv_ztoh(const double complex (*z)[2], double complex (*h)[2])
    void vnaconv_ztos(const double complex (*z)[2], double complex (*s)[2],
                      const double complex *z0)
    void vnaconv_ztot(const double complex (*z)[2], double complex (*t)[2],
                      const double complex *z0)
    void vnaconv_ztou(const double complex (*z)[2], double complex (*u)[2],
                      const double complex *z0)
    void vnaconv_ztoy(const double complex (*z)[2], double complex (*y)[2])
    void vnaconv_ztozi(const double complex (*z)[2], double complex *zi,
                       const double complex *z0)

    #
    # NxN conversions
    #
    void vnaconv_stoyn(const double complex *s, double complex *y,
                       const double complex *z0, int n)
    void vnaconv_stozin(const double complex *s, double complex *zi,
                        const double complex *z0, int n)
    void vnaconv_stozn(const double complex *s, double complex *z,
                       const double complex *z0, int n)
    void vnaconv_ytosn(const double complex *y, double complex *s,
                       const double complex *z0, int n)
    void vnaconv_ytozin(const double complex *y, double complex *zi,
                        const double complex *z0, int n)
    void vnaconv_ytozn(const double complex *y, double complex *s, int n)
    void vnaconv_ztosn(const double complex *z, double complex *s,
                       const double complex *z0, int n)
    void vnaconv_ztoyn(const double complex *z, double complex *y, int n)
    void vnaconv_ztozin(const double complex *z, double complex *zi,
                        const double complex *z0, int n)
