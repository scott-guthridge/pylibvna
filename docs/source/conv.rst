libvna.conv: Network Parameter Conversion
=========================================

.. toctree::
   :maxdepth: 2

The functions in this module convert between different representations
of electrical n-port network parameters.  Supported parameter
types include **s** (scattering), **t** (scattering-transfer), **u**
(inverse scattering-transfer), **z** (impedance), **y** (admittance),
**h** (hybrid), **g** (inverse hybrid), **a** (ABCD), and **b**
(inverse ABCD).  While **s**, **z** and **y** parameters are defined for
any number of ports (NxN matrices), **t**, **u**, **h**, **g**, **a** and
**b** parameters are defined for two ports networks (2x2 matrices) only.

The **s**, **t** and **u** parameters are defined in relation to given
reference impedances, **z0**.  When converting between these types and
the others, the reference impedances must be known.  The *z0* parameter
can be a scalar, or a vector with length equal to the number of ports.
It can be real or complex.  If not specified, it defaults to 50 ohms.

All functions take an array-like input parameter, *array*.
This parameter must have at least two dimensions, with the final two
dimensions square and consistent with the parameter type.  If the array
has additional dimensions on the left, the matrix conversion functions
operate on each matrix of the larger array and return a result in the
same shape as the input.

In addition to the matrix conversions, the module provides functions that
convert from each parameter type to a vector of input impedances giving
the impedance looking into each port if the other ports are terminated
in the reference impedances.


Two-Port Matrix Conversions
---------------------------

.. autofunction:: libvna.conv.atob
.. autofunction:: libvna.conv.atog
.. autofunction:: libvna.conv.atoh
.. autofunction:: libvna.conv.atos
.. autofunction:: libvna.conv.atot
.. autofunction:: libvna.conv.atou
.. autofunction:: libvna.conv.atoy
.. autofunction:: libvna.conv.atoz
.. autofunction:: libvna.conv.btoa
.. autofunction:: libvna.conv.btog
.. autofunction:: libvna.conv.btoh
.. autofunction:: libvna.conv.btos
.. autofunction:: libvna.conv.btot
.. autofunction:: libvna.conv.btou
.. autofunction:: libvna.conv.btoy
.. autofunction:: libvna.conv.btoz
.. autofunction:: libvna.conv.gtoa
.. autofunction:: libvna.conv.gtob
.. autofunction:: libvna.conv.gtoh
.. autofunction:: libvna.conv.gtos
.. autofunction:: libvna.conv.gtot
.. autofunction:: libvna.conv.gtou
.. autofunction:: libvna.conv.gtoy
.. autofunction:: libvna.conv.gtoz
.. autofunction:: libvna.conv.htoa
.. autofunction:: libvna.conv.htob
.. autofunction:: libvna.conv.htog
.. autofunction:: libvna.conv.htos
.. autofunction:: libvna.conv.htot
.. autofunction:: libvna.conv.htou
.. autofunction:: libvna.conv.htoy
.. autofunction:: libvna.conv.htoz
.. autofunction:: libvna.conv.stoa
.. autofunction:: libvna.conv.stob
.. autofunction:: libvna.conv.stog
.. autofunction:: libvna.conv.stoh
.. autofunction:: libvna.conv.stot
.. autofunction:: libvna.conv.stou
.. autofunction:: libvna.conv.ttoa
.. autofunction:: libvna.conv.ttob
.. autofunction:: libvna.conv.ttog
.. autofunction:: libvna.conv.ttoh
.. autofunction:: libvna.conv.ttos
.. autofunction:: libvna.conv.ttou
.. autofunction:: libvna.conv.ttoy
.. autofunction:: libvna.conv.ttoz
.. autofunction:: libvna.conv.utoa
.. autofunction:: libvna.conv.utob
.. autofunction:: libvna.conv.utog
.. autofunction:: libvna.conv.utoh
.. autofunction:: libvna.conv.utos
.. autofunction:: libvna.conv.utot
.. autofunction:: libvna.conv.utoy
.. autofunction:: libvna.conv.utoz
.. autofunction:: libvna.conv.ytoa
.. autofunction:: libvna.conv.ytob
.. autofunction:: libvna.conv.ytog
.. autofunction:: libvna.conv.ytoh
.. autofunction:: libvna.conv.ytot
.. autofunction:: libvna.conv.ytou
.. autofunction:: libvna.conv.ztoa
.. autofunction:: libvna.conv.ztob
.. autofunction:: libvna.conv.ztog
.. autofunction:: libvna.conv.ztoh
.. autofunction:: libvna.conv.ztot
.. autofunction:: libvna.conv.ztou

N-Port Matrix Conversions
-------------------------

.. autofunction:: libvna.conv.stoy
.. autofunction:: libvna.conv.stoz
.. autofunction:: libvna.conv.ytos
.. autofunction:: libvna.conv.ytoz
.. autofunction:: libvna.conv.ztos
.. autofunction:: libvna.conv.ztoy

Two-Port Matrix to Input Impedance
----------------------------------

.. autofunction:: libvna.conv.atozi
.. autofunction:: libvna.conv.btozi
.. autofunction:: libvna.conv.gtozi
.. autofunction:: libvna.conv.htozi
.. autofunction:: libvna.conv.ttozi
.. autofunction:: libvna.conv.utozi

N-Port Matrix to Input Impedance
--------------------------------

.. autofunction:: libvna.conv.stozi
.. autofunction:: libvna.conv.ytozi
.. autofunction:: libvna.conv.ztozi

Mathematical Model
------------------

Let:

* :math:`a_1` and :math:`a_2` be the incident voltages [#]_ into ports 1 and 2,
* :math:`b_1` and :math:`b_2` be the reflected voltages out of ports 1 and 2,
* :math:`v_1` and :math:`v_2` be the voltages at ports 1 and 2,
* :math:`i_1` and :math:`i_2` be the currents into ports 1 and 2, and
* :math:`Z_1` and :math:`Z_2` be the system impedances the device sees looking out of its ports

The relationships between the above for port 1 are:

.. math::

    \begin{aligned}
    a_1 &= \frac{1}{2} K_1 (v_1 + Z_1   i_1) \\
    b_1 &= \frac{1}{2} K_1 (v_1 - Z_1^* i_1)
    \end{aligned}
    \quad
    \begin{aligned}
    v_1 &= \frac{Z_1^* a_1 + Z_1 b_1}{K_1 \, \Re(Z_1)} \\
    i_1 &= \frac{a_1 - b_1}{K_1 \, \Re(Z_1)} \\
    \end{aligned}

and for port 2:

.. math::

    \begin{aligned}
    a_2 &= \frac{1}{2} K_2 (v_2 + Z_2   i_2) \\
    b_2 &= \frac{1}{2} K_2 (v_2 - Z_2^* i_2)
    \end{aligned}
    \quad
    \begin{aligned}
    v_2 &= \frac{Z_2^* a_2 + Z_2 b_2}{K_2 \, \Re(Z_2)} \\
    i_2 &= \frac{a_2 - b_2}{K_2 \, \Re(Z_2)}
    \end{aligned}

where :math:`K_1` and :math:`K_2` are unit scaling constants:

.. math::

    \begin{aligned}
    K_1 &= \frac{1}{\sqrt{\left|\Re(Z_1)\right|}} \\
    K_2 &= \frac{1}{\sqrt{\left|\Re(Z_2)\right|}}
    \end{aligned}

and :math:`*` is the conjugate operator.  The same pattern applies for
additional ports.

We can now show the definitions of each network parameter representation.
The **s** (scattering) parameters are defined as:

.. math::

   \begin{bmatrix}
   b_1 \\
   b_2
   \end{bmatrix} =
   \begin{bmatrix}
   s_{11} & s_{12} \\
   s_{21} & s_{22}
   \end{bmatrix}
   \begin{bmatrix}
   a_1 \\
   a_2
   \end{bmatrix}

The **t** (scattering-transfer) parameters are defined as:

.. math::

   \begin{bmatrix}
   b_1 \\
   a_1
   \end{bmatrix} =
   \begin{bmatrix}
   t_{11} & t_{12} \\
   t_{21} & t_{22}
   \end{bmatrix}
   \begin{bmatrix}
   a_2 \\
   b_2
   \end{bmatrix}

The **t** parameters for a cacade of two-port networks is the
left-to-right matrix product of the **t** parameters of each successive
stage.

The **u** (inverse scattering-transfer) parameters are defined as:

.. math::

   \begin{bmatrix}
   a_2 \\
   b_2
   \end{bmatrix} =
   \begin{bmatrix}
   u_{11} & u_{12} \\
   u_{21} & u_{22}
   \end{bmatrix}
   \begin{bmatrix}
   b_1 \\
   a_1
   \end{bmatrix}

The **u** parameters for a cacade of two-port networks is the
right-to-left matrix product of the **u** parameters of each successive
stage.

The **z** (impedance) parameters are defined as:

.. math::

   \begin{bmatrix}
   v_1 \\
   v_2
   \end{bmatrix} =
   \begin{bmatrix}
   z_{11} & z_{12} \\
   z_{21} & z_{22}
   \end{bmatrix}
   \begin{bmatrix}
   i_1 \\
   i_2
   \end{bmatrix}

The **y** (admittance) parameters are defined as:

.. math::

   \begin{bmatrix}
   i_1 \\
   i_2
   \end{bmatrix} =
   \begin{bmatrix}
   y_{11} & y_{12} \\
   y_{21} & y_{22}
   \end{bmatrix}
   \begin{bmatrix}
   v_1 \\
   v_2
   \end{bmatrix}

The **h** (hybrid) parameters are defined as:

.. math::

   \begin{bmatrix}
   v_1 \\
   i_2
   \end{bmatrix} =
   \begin{bmatrix}
   h_{11} & h_{12} \\
   h_{21} & h_{22}
   \end{bmatrix}
   \begin{bmatrix}
   i_1 \\
   v_2
   \end{bmatrix}

The **g** (inverse hybrid) parameters are defined as:

.. math::

   \begin{bmatrix}
   i_1 \\
   v_2
   \end{bmatrix} =
   \begin{bmatrix}
   g_{11} & g_{12} \\
   g_{21} & g_{22}
   \end{bmatrix}
   \begin{bmatrix}
   v_1 \\
   i_2
   \end{bmatrix}

The **a** (ABCD) parameters are defined as:

.. math::

   \begin{bmatrix}
   v_1 \\
   i_1
   \end{bmatrix} =
   \begin{bmatrix}
   a_{11} & a_{12} \\
   a_{21} & a_{22}
   \end{bmatrix}
   \begin{bmatrix}
   v_2 \\
   -i_2
   \end{bmatrix}

The **a** parameters for a cascade of two-port networks is the
left-to-right matrix product of the **a** parameters for each
successive stage.  Don't confuse the **a** matrix with the :math:`a_1`
and :math:`a_2` variables above.

The **b** (inverse ABCD) parameters are defined as:

.. math::

   \begin{bmatrix}
   v_2 \\
   -i_2
   \end{bmatrix} =
   \begin{bmatrix}
   b_{11} & b_{12} \\
   b_{21} & b_{22}
   \end{bmatrix}
   \begin{bmatrix}
   v_1 \\
   i_1
   \end{bmatrix}

The **b** parameters for a cascade of two-port networks is the
right-to-left matrix product of the **b** parameters for each
successive stage.  Don't confuse the **b** matrix with the :math:`b_1`
and :math:`b_2` variables above.

.. rubric:: Footnotes

.. [#] Note that for :math:`a1, a2, b1` and :math:`b2`, we've used
       "voltage" imprecisely.  These variables are really defined as
       root power in units of :math:`\text{Watt}^{1/2}`.  In most cases,
       this distinction is unimportant because the units divide out.

Examples
--------

Use the conversion functions to analyze the 75 ohm to 50 ohm impedance
matching L-pad.

.. code-block:: python

    import libvna.conv as vc
    from math import sqrt

    # System Impedances
    Z1 = 75
    Z2 = 50

    # Resistor values for the impedance-matching L-pad
    R1 = sqrt(Z1) * sqrt(Z1 - Z2)
    R2 = sqrt(Z1) * Z2 / sqrt(Z1 - Z2)

    # Z-parameters of the L-pad
    Z = [[R1+R2, R2],
         [R2,    R2]]

    # Convert from Z to S parameters
    S = vc.ztos(Z, [Z1, Z2])
    print("S =\n"
          f"  {S[0, 0].real:8.5f}{S[0, 0].imag:+8.5f}j"
          f"  {S[0, 1].real:8.5f}{S[0, 1].imag:+8.5f}j\n"
          f"  {S[1, 0].real:8.5f}{S[1, 0].imag:+8.5f}j"
          f"  {S[1, 1].real:8.5f}{S[1, 1].imag:+8.5f}j")

    # Convert from S parameters to input impedances at each port
    Zi = vc.stozi(S, [Z1, Z2])
    print("Zi =\n"
          f"  {Zi[0].real:8.5f}{Zi[0].imag:+8.5f}j"
          f"  {Zi[1].real:8.5f}{Zi[1].imag:+8.5f}j")

Result:

.. code-block::

    S =
      -0.00000+0.00000j   0.51764+0.00000j
       0.51764+0.00000j   0.00000+0.00000j
    Zi =
      75.00000+0.00000j  50.00000+0.00000j

From the S-parameters, we can see that L-pad attenuates the signal in
both directions by 0.51764 (5.719 dB) and has no reflection.  When we
convert to input impedances, we can see that the L-pad has the expected
impedances at each port.
