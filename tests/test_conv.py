#!/usr/bin/python3
"""
Test the libvna.conv module.
"""

import unittest
import numpy as np
import libvna.conv as vc

TRIALS = 1

M_SQRT1_2 = 1.41421356237309504880

def crandn():
    """
    Return a gaussian complex random number
    """
    return M_SQRT1_2 * complex(np.random.normal(), np.random.normal())

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

if __name__ == '__main__':
    unittest.main()
