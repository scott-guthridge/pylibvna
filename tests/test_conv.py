#!/usr/bin/python3
"""
Test the libvna.conv module.
"""

import unittest
import numpy as np
import libvna.conv as vc

TRIALS = 20

M_SQRT1_2 = 1.41421356237309504880


def crandn(size=None):
    """
    Return a gaussian complex random number
    """
    return M_SQRT1_2 * (np.random.normal(size=size) +
                        1j * np.random.normal(size=size))


class TestModule(unittest.TestCase):
    def test_2x2(self):
        """
        Test all the 2x2 conversion cases
        """
        for _ in range(TRIALS):
            #
            # Make random system impedances and their conjugates.
            #
            Z1 = crandn()
            Z2 = crandn()
            Z1c = np.conjugate(Z1)
            Z2c = np.conjugate(Z2)

            #
            # Make the scaling factors to put a and b into units of
            # of sqrt(Watt).
            #
            k1i = np.sqrt(abs(np.real(Z1)))
            k2i = np.sqrt(abs(np.real(Z2)))

            #
            # Make random incident power, a random S matrix, and
            # from these, calculate the reflected power.
            #
            a1 = crandn()
            a2 = crandn()
            s = np.array([[crandn(), crandn()], [crandn(), crandn()]])
            b1 = s[0, 0] * a1 + s[0, 1] * a2
            b2 = s[1, 1] * a2 + s[1, 0] * a1

            #
            # Calculate voltage at and current into each DUT port.
            #
            v1 = k1i * (Z1c * a1 + Z1 * b1) / np.real(Z1)
            v2 = k2i * (Z2c * a2 + Z2 * b2) / np.real(Z2)
            i1 = k1i * (a1 - b1) / np.real(Z1)
            i2 = k2i * (a2 - b2) / np.real(Z2)

            #
            # Calculate input impedance looking into each port when
            # the other ports are terminated in the system impendances.
            #
            zi = np.array([(s[0, 0] * Z1 + Z1c) / (1.0 - s[0, 0]),
                           (s[1, 1] * Z2 + Z2c) / (1.0 - s[1, 1])])

            #
            # Convert s to t and verify against the defition of t.
            #
            t = vc.stot(s)
            c = np.array([[a2], [b2]])
            d = np.array([[b1], [a1]])
            self.assertTrue(np.allclose(np.matmul(t, c), d))

            #
            # Convert s to u and verify against the defition of u.
            #
            u = vc.stou(s)
            c = np.array([[b1], [a1]])
            d = np.array([[a2], [b2]])
            self.assertTrue(np.allclose(np.matmul(u, c), d))

            #
            # Convert s to z and verify against the defition of z.
            #
            z = vc.stoz(s, [Z1, Z2])
            c = np.array([[i1], [i2]])
            d = np.array([[v1], [v2]])
            self.assertTrue(np.allclose(np.matmul(z, c), d))

            #
            # Convert s to y and verify against the defition of y.
            #
            y = vc.stoy(s, [Z1, Z2])
            c = np.array([[v1], [v2]])
            d = np.array([[i1], [i2]])
            self.assertTrue(np.allclose(np.matmul(y, c), d))

            #
            # Convert s to h and verify against the defition of h.
            #
            h = vc.stoh(s, [Z1, Z2])
            c = np.array([[i1], [v2]])
            d = np.array([[v1], [i2]])
            self.assertTrue(np.allclose(np.matmul(h, c), d))

            #
            # Convert s to g and verify against the defition of g.
            #
            g = vc.stog(s, [Z1, Z2])
            c = np.array([[v1], [i2]])
            d = np.array([[i1], [v2]])
            self.assertTrue(np.allclose(np.matmul(g, c), d))

            #
            # Convert s to a and verify against the defition of a.
            #
            a = vc.stoa(s, [Z1, Z2])
            c = np.array([[v2], [-i2]])
            d = np.array([[v1], [i1]])
            self.assertTrue(np.allclose(np.matmul(a, c), d))

            #
            # Convert s to b and verify against the defition of b.
            #
            b = vc.stob(s, [Z1, Z2])
            c = np.array([[v1], [i1]])
            d = np.array([[v2], [-i2]])
            self.assertTrue(np.allclose(np.matmul(b, c), d))

            #
            # Convert s to Zin and verify.
            #
            x = vc.stozi(s, [Z1, Z2])
            self.assertTrue(np.allclose(x, zi))

            #
            # Convert t to each other parameter type and verify.
            #
            x = vc.ttos(t)
            self.assertTrue(np.allclose(x, s))
            x = vc.ttou(t)
            self.assertTrue(np.allclose(x, u))
            x = vc.ttoz(t, [Z1, Z2])
            self.assertTrue(np.allclose(x, z))
            x = vc.ttoy(t, [Z1, Z2])
            self.assertTrue(np.allclose(x, y))
            x = vc.ttoh(t, [Z1, Z2])
            self.assertTrue(np.allclose(x, h))
            x = vc.ttog(t, [Z1, Z2])
            self.assertTrue(np.allclose(x, g))
            x = vc.ttoa(t, [Z1, Z2])
            self.assertTrue(np.allclose(x, a))
            x = vc.ttob(t, [Z1, Z2])
            self.assertTrue(np.allclose(x, b))
            x = vc.ttozi(t, [Z1, Z2])
            self.assertTrue(np.allclose(x, zi))

            #
            # Convert u to each other parameter type and verify.
            #
            x = vc.utos(u)
            self.assertTrue(np.allclose(x, s))
            x = vc.utot(u)
            self.assertTrue(np.allclose(x, t))
            x = vc.utoz(u, [Z1, Z2])
            self.assertTrue(np.allclose(x, z))
            x = vc.utoy(u, [Z1, Z2])
            self.assertTrue(np.allclose(x, y))
            x = vc.utoh(u, [Z1, Z2])
            self.assertTrue(np.allclose(x, h))
            x = vc.utog(u, [Z1, Z2])
            self.assertTrue(np.allclose(x, g))
            x = vc.utoa(u, [Z1, Z2])
            self.assertTrue(np.allclose(x, a))
            x = vc.utob(u, [Z1, Z2])
            self.assertTrue(np.allclose(x, b))
            x = vc.utozi(u, [Z1, Z2])
            self.assertTrue(np.allclose(x, zi))

            #
            # Convert z to each other parmeter type and verify.
            #
            x = vc.ztos(z, [Z1, Z2])
            self.assertTrue(np.allclose(x, s))
            x = vc.ztot(z, [Z1, Z2])
            self.assertTrue(np.allclose(x, t))
            x = vc.ztou(z, [Z1, Z2])
            self.assertTrue(np.allclose(x, u))
            x = vc.ztoy(z)
            self.assertTrue(np.allclose(x, y))
            x = vc.ztoh(z)
            self.assertTrue(np.allclose(x, h))
            x = vc.ztog(z)
            self.assertTrue(np.allclose(x, g))
            x = vc.ztoa(z)
            self.assertTrue(np.allclose(x, a))
            x = vc.ztob(z)
            self.assertTrue(np.allclose(x, b))
            x = vc.ztozi(z, [Z1, Z2])
            self.assertTrue(np.allclose(x, zi))

            #
            # Convert y to each other parameter type and verify.
            #
            x = vc.ytos(y, [Z1, Z2])
            self.assertTrue(np.allclose(x, s))
            x = vc.ytot(y, [Z1, Z2])
            self.assertTrue(np.allclose(x, t))
            x = vc.ytou(y, [Z1, Z2])
            self.assertTrue(np.allclose(x, u))
            x = vc.ytoz(y)
            self.assertTrue(np.allclose(x, z))
            x = vc.ytoh(y)
            self.assertTrue(np.allclose(x, h))
            x = vc.ytog(y)
            self.assertTrue(np.allclose(x, g))
            x = vc.ytoa(y)
            self.assertTrue(np.allclose(x, a))
            x = vc.ytob(y)
            self.assertTrue(np.allclose(x, b))
            x = vc.ytozi(y, [Z1, Z2])
            self.assertTrue(np.allclose(x, zi))

            #
            # Convert h to each other parameter type and verify.
            #
            x = vc.htos(h, [Z1, Z2])
            self.assertTrue(np.allclose(x, s))
            x = vc.htot(h, [Z1, Z2])
            self.assertTrue(np.allclose(x, t))
            x = vc.htou(h, [Z1, Z2])
            self.assertTrue(np.allclose(x, u))
            x = vc.htoz(h)
            self.assertTrue(np.allclose(x, z))
            x = vc.htoy(h)
            self.assertTrue(np.allclose(x, y))
            x = vc.htog(h)
            self.assertTrue(np.allclose(x, g))
            x = vc.htoa(h)
            self.assertTrue(np.allclose(x, a))
            x = vc.htob(h)
            self.assertTrue(np.allclose(x, b))
            x = vc.htozi(h, [Z1, Z2])
            self.assertTrue(np.allclose(x, zi))

            #
            # Convert g to each other paramter type and verify.
            #
            x = vc.gtos(g, [Z1, Z2])
            self.assertTrue(np.allclose(x, s))
            x = vc.gtot(g, [Z1, Z2])
            self.assertTrue(np.allclose(x, t))
            x = vc.gtou(g, [Z1, Z2])
            self.assertTrue(np.allclose(x, u))
            x = vc.gtoz(g)
            self.assertTrue(np.allclose(x, z))
            x = vc.gtoy(g)
            self.assertTrue(np.allclose(x, y))
            x = vc.gtoh(g)
            self.assertTrue(np.allclose(x, h))
            x = vc.gtoa(g)
            self.assertTrue(np.allclose(x, a))
            x = vc.gtob(g)
            self.assertTrue(np.allclose(x, b))
            x = vc.gtozi(g, [Z1, Z2])
            self.assertTrue(np.allclose(x, zi))

            #
            # Convert a to each other parameter type and verify.
            #
            x = vc.atos(a, [Z1, Z2])
            self.assertTrue(np.allclose(x, s))
            x = vc.atot(a, [Z1, Z2])
            self.assertTrue(np.allclose(x, t))
            x = vc.atou(a, [Z1, Z2])
            self.assertTrue(np.allclose(x, u))
            x = vc.atoz(a)
            self.assertTrue(np.allclose(x, z))
            x = vc.atoy(a)
            self.assertTrue(np.allclose(x, y))
            x = vc.atoh(a)
            self.assertTrue(np.allclose(x, h))
            x = vc.atog(a)
            self.assertTrue(np.allclose(x, g))
            x = vc.atob(a)
            self.assertTrue(np.allclose(x, b))
            x = vc.atozi(a, [Z1, Z2])
            self.assertTrue(np.allclose(x, zi))

            #
            # Convert b to each other parameter type and verify.
            #
            x = vc.btos(b, [Z1, Z2])
            self.assertTrue(np.allclose(x, s))
            x = vc.btot(b, [Z1, Z2])
            self.assertTrue(np.allclose(x, t))
            x = vc.btou(b, [Z1, Z2])
            self.assertTrue(np.allclose(x, u))
            x = vc.btoz(b)
            self.assertTrue(np.allclose(x, z))
            x = vc.btoy(b)
            self.assertTrue(np.allclose(x, y))
            x = vc.btoh(b)
            self.assertTrue(np.allclose(x, h))
            x = vc.btog(b)
            self.assertTrue(np.allclose(x, g))
            x = vc.btoa(b)
            self.assertTrue(np.allclose(x, a))
            x = vc.btozi(b, [Z1, Z2])
            self.assertTrue(np.allclose(x, zi))

    def test_3x3(self):
        """
        Test all the NxN conversion cases using 3 ports.
        """
        for _ in range(TRIALS):
            #
            # Make random system impedances and their conjugates.
            #
            Z1 = crandn()
            Z2 = crandn()
            Z3 = crandn()
            Z1c = np.conjugate(Z1)
            Z2c = np.conjugate(Z2)
            Z3c = np.conjugate(Z3)

            #
            # Make the scaling factors to put a and b into units of
            # of sqrt(Watt).
            #
            k1i = np.sqrt(abs(np.real(Z1)))
            k2i = np.sqrt(abs(np.real(Z2)))
            k3i = np.sqrt(abs(np.real(Z3)))

            #
            # Make random incident power, a random S matrix, and
            # from these, calculate the reflected power.
            #
            a1 = crandn()
            a2 = crandn()
            a3 = crandn()
            s = np.array([[crandn(), crandn(), crandn()],
                          [crandn(), crandn(), crandn()],
                          [crandn(), crandn(), crandn()]])
            b1 = s[0, 0] * a1 + s[0, 1] * a2 + s[0, 2] * a3
            b2 = s[1, 0] * a1 + s[1, 1] * a2 + s[1, 2] * a3
            b3 = s[2, 0] * a1 + s[2, 1] * a2 + s[2, 2] * a3

            #
            # Calculate voltage at and current into each DUT port.
            #
            v1 = k1i * (Z1c * a1 + Z1 * b1) / np.real(Z1)
            v2 = k2i * (Z2c * a2 + Z2 * b2) / np.real(Z2)
            v3 = k3i * (Z3c * a3 + Z3 * b3) / np.real(Z3)
            i1 = k1i * (a1 - b1) / np.real(Z1)
            i2 = k2i * (a2 - b2) / np.real(Z2)
            i3 = k3i * (a3 - b3) / np.real(Z3)

            #
            # Calculate input impedance looking into each port when
            # the other ports are terminated in the system impendances.
            #
            zi = np.array([(s[0, 0] * Z1 + Z1c) / (1.0 - s[0, 0]),
                           (s[1, 1] * Z2 + Z2c) / (1.0 - s[1, 1]),
                           (s[2, 2] * Z3 + Z3c) / (1.0 - s[2, 2])])

            #
            # Convert s to z and verify against the defition of z.
            #
            z = vc.stoz(s, [Z1, Z2, Z3])
            c = np.array([[i1], [i2], [i3]])
            d = np.array([[v1], [v2], [v3]])
            self.assertTrue(np.allclose(np.matmul(z, c), d))

            #
            # Convert s to y and verify against the defition of y.
            #
            y = vc.stoy(s, [Z1, Z2, Z3])
            c = np.array([[v1], [v2], [v3]])
            d = np.array([[i1], [i2], [i3]])
            self.assertTrue(np.allclose(np.matmul(y, c), d))

            #
            # Convert s to Zin and verify.
            #
            x = vc.stozi(s, [Z1, Z2, Z3])
            self.assertTrue(np.allclose(x, zi))

            #
            # Convert z to each other parmeter type and verify.
            #
            x = vc.ztos(z, [Z1, Z2, Z3])
            self.assertTrue(np.allclose(x, s))
            x = vc.ztoy(z)
            self.assertTrue(np.allclose(x, y))
            x = vc.ztozi(z, [Z1, Z2, Z3])
            self.assertTrue(np.allclose(x, zi))

            #
            # Convert y to each other parameter type and verify.
            #
            x = vc.ytos(y, [Z1, Z2, Z3])
            self.assertTrue(np.allclose(x, s))
            x = vc.ytoz(y)
            self.assertTrue(np.allclose(x, z))
            x = vc.ytozi(y, [Z1, Z2, Z3])
            self.assertTrue(np.allclose(x, zi))

    def test_z0_dimensions(self):
        # Test 75-50 ohm L-pad at 50 ohms.
        z = [[129.9038105676658, 86.60254037844386],
             [86.60254037844386, 86.60254037844386]],
        s = vc.ztos(z)
        e = [[0.2000000000000000, 0.5071796769724491],
             [0.5071796769724491, -0.05358983848622454]]
        self.assertTrue(np.allclose(s, e))

        # Test 75-50 ohm L-pad at 75, 50 ohms.
        s = vc.ztos(z, [75, 50])
        e = [[0, 0.5176380902050415],
             [0.5176380902050415, 0]]
        self.assertTrue(np.allclose(s, e))

        # Expect errors if z vector has wrong length or dimensions.
        with self.assertRaises(ValueError):
            vc.ztos(z, [50])
        with self.assertRaises(ValueError):
            vc.ztos(z, [75, 50, 110])
        with self.assertRaises(ValueError):
            vc.ztos(z, [[75, 50], [50, 75]])

        # Test with extra dimensions.
        z = [
                [
                    # 75-50 LPAD
                    [[129.9038105676658, 86.60254037844386],
                     [86.60254037844386, 86.60254037844386]],

                    # 50-75 LPAD
                    [[-70.71067811865475j, -106.0660171779821j],
                     [-106.0660171779821j, -106.0660171779821j]],

                    # s 1, 2j, 3, -4j @ 50 ohms
                    [[-116.66666666666667+16.66666666666667j,
                      -33.33333333333333], [50.0j, -50.0]]
                ],
                [
                    # 1/2 attentutor @ 50
                    [[83.33333333, 66.66666667],
                     [66.66666667, 83.33333333]],

                    # 1/10 attenuator @ 50
                    [[51.0101, 10.101],
                     [10.101, 51.0101]],

                    # 50, 75 double reflect
                    [[50, 0],
                     [0, 75]]
                ]
            ]

        # expected S at 50 ohms
        e = [
                [
                    [[0.2000000000000000, 0.5071796769724491],
                     [0.5071796769724491, -0.05358983848622454]],
                    [[-0.0666666666666667+0.1885618083164127j,
                      0.8000000000000000-0.5656854249492380j],
                     [0.8000000000000000-0.5656854249492380j,
                      0.2000000000000000]],
                    [[1.000000000000000, 2.000000000000000j],
                     [3.000000000000000, -4.000000000000000j]]
                ],
                [
                    [[0, 0.5000000000000000],
                     [0.5000000000000000, 0]],
                    [[0, 0.1000000000000000],
                     [0.1000000000000000, 0]],
                    [[0, 0],
                     [0, 0.2000000000000000]]
                ]
            ]
        s = vc.ztos(z)
        self.assertTrue(np.allclose(s, e))

        # expected S at at 75, 50 ohms
        e = [
                [
                    [[0, 0.5176380902050415],
                     [0.5176380902050415, 0]],
                    [[-0.2697095435684647+0.1760431820381446j,
                      0.7927809126020244-0.5174591624272165j],
                     [0.7927809126020244-0.5174591624272165j,
                      0.2697095435684647-0.1760431820381446j]],
                    [[1.000000000000000, 2.449489742783178j],
                     [3.674234614174767, -2.500000000000000j]]
                ], [
                    [[-0.2000000000000000, 0.4898979485566356],
                     [0.4898979485566356, 0.05000000000000000]],
                    [[-0.2000000000000000, 0.09797958971132712],
                     [0.09797958971132712, 0.002000000000000000]],
                    [[-0.2000000000000000, 0],
                     [0, 0.2000000000000000]]
                ]
            ]
        s = vc.ztos(z, [75, 50])
        self.assertTrue(np.allclose(s, e))

        # expected S at various z0
        e = [
                [
                    # 75-50 LPAD at 75, 50
                    [[0, 0.5176380902050415],
                     [0.5176380902050415, 0]],
                    # 50-75 LPAD at 50, 75
                    [[0, 0.8164965809277260-0.5773502691896258j],
                     [0.8164965809277260-0.5773502691896258j, 0]],
                    # s 1, 2j, 3, -4j @ 50 ohms
                    [[1.000000000000000, 2.000000000000000j],
                     [3.000000000000000, -4.000000000000000j]]
                ],
                [
                    # 1/2 attentutor @ 50
                    [[0, 0.5000000000000000],
                     [0.5000000000000000, 0]],
                    # 1/10 attenuator @ 50
                    [[0, 0.1000000000000000],
                     [0.1000000000000000, 0]],
                    # 50, 75 double reflect at 50, 75
                    [[0, 0],
                     [0, 0]]
                ]
            ]
        s = vc.ztos(z, [[[75, 50], [50, 75], [50, 50]],
                        [[50, 50], [50, 50], [50, 75]]])
        self.assertTrue(np.allclose(s, e))

        # Expect errors if z vector has wrong dimensions.
        with self.assertRaises(ValueError):
            vc.ztos(z, [[[75, 50], [50, 75]],
                        [[50, 50], [50, 75]]])

    def renormalize_helper(self, ports):
        """
        Test renormalizing conversions for a given set of ports.
        """
        z01v = crandn(size=(ports,))
        z02v = crandn(size=(ports,))
        z01  = np.diag(z01v)
        z02  = np.diag(z02v)
        z01c = np.diag(np.conj(z01v))
        z02c = np.diag(np.conj(z02v))
        zr1v = np.real(z01v)
        zr2v = np.real(z02v)
        k1ri = np.diag(np.sqrt(abs(zr1v)) / zr1v)
        k2h  = np.diag(0.5 * np.sqrt(abs(zr2v))**(-1))
        a1 = crandn(size=(ports, ports))
        s1 = crandn(size=(ports, ports))
        b1 = s1 @ a1
        v  = k1ri @ (z01c @ a1 + z01 @ b1)
        i  = k1ri @ (a1 - b1)
        a2 = k2h @ (v + z02  @ i)
        b2 = k2h @ (v - z02c @ i)
        s2 = np.linalg.solve(a2.T, b2.T).T

        # check s1 -> s2
        x = vc.stos(s1, z01v, z02v)
        self.assertTrue(np.allclose(x, s2))

        # all other conversions are 2x2 only
        if ports != 2:
            return

        # check s1 -> t1
        t1 = vc.stot(s1)
        c = np.vstack((a1[1], b1[1]))
        d = np.vstack((b1[0], a1[0]))
        self.assertTrue(np.allclose(np.matmul(t1, c), d))

        # check s1 -> u1
        u1 = vc.stou(s1)
        c = np.vstack((b1[0], a1[0]))
        d = np.vstack((a1[1], b1[1]))
        self.assertTrue(np.allclose(np.matmul(u1, c), d))

        # check s1 -> t2
        t2 = vc.stot(s1, z01v, z02v)
        c = np.vstack((a2[1], b2[1]))
        d = np.vstack((b2[0], a2[0]))
        self.assertTrue(np.allclose(np.matmul(t2, c), d))

        # check s1 -> u2
        u2 = vc.stou(s1, z01v, z02v)
        c = np.vstack((b2[0], a2[0]))
        d = np.vstack((a2[1], b2[1]))
        self.assertTrue(np.allclose(np.matmul(u2, c), d))

        # check t1 -> s2
        x = vc.ttos(t1, z01v, z02v)
        self.assertTrue(np.allclose(x, s2))

        # check t1 -> t2
        x = vc.ttot(t1, z01v, z02v)
        self.assertTrue(np.allclose(x, t2))

        # check t1 -> u2
        x = vc.ttou(t1, z01v, z02v)
        self.assertTrue(np.allclose(x, u2))

        # check u1 -> s2
        x = vc.utos(u1, z01v, z02v)
        self.assertTrue(np.allclose(x, s2))

        # check u1 -> t2
        x = vc.utot(u1, z01v, z02v)
        self.assertTrue(np.allclose(x, t2))

        # check u1 -> u2
        x = vc.utou(u1, z01v, z02v)
        self.assertTrue(np.allclose(x, u2))

    def test_try(self):
        """
        Test renormalizing conversions.
        """
        for _ in range(TRIALS):
            for ports in range(1, 6):
                self.renormalize_helper(ports)


if __name__ == "__main__":
    unittest.main()
