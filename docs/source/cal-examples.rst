Calibration Examples
====================

One-Port Reflect Only
---------------------

Create the calibration from measurements of short, open and load standards.

.. literalinclude:: ../../examples/1x1-calibrate.py
   :language: python

Apply the calibration to a device under test.

.. literalinclude:: ../../examples/1x1-apply.py
   :language: python


One-Port Reflect Only with Measured Standards
---------------------------------------------

Create the calibration from measurements of imperfect standards that have
been measured on a more trusted instrument, loading the S-parameters of
each standard from a file.

.. literalinclude:: ../../examples/1x1m-calibrate.py
   :language: python

Apply the calibration to a device under test.

.. literalinclude:: ../../examples/1x1m-apply.py
   :language: python


2x1 SOLT
--------

Example of SOLT calibration for a VNA that measures only :math:`S_{11}`
and :math:`S_{21}`

.. literalinclude:: ../../examples/2x1-calibrate.py
   :language: python

Example of applying the calibration to a device under test, where we
first measure the raw values of :math:`S_{11}` and :math:`S_{21}`,
then exchange the probes and measure the raw values of :math:`S_{22}`
and :math:`S_{12}`.  With all four values and the calibration, we can
solve for the full S-parameters of the device.

.. literalinclude:: ../../examples/2x1-apply.py
   :language: python


2x2 with A & B Measurements
---------------------------

In this example, our VNA measures full S parameters.  In addition, it
measures the transmitted power (`a` matrix) as well as the reflected power
(`b` matrix).  Having the `a` matrix makes it possible to compensate for
errors in the RF switch(es) without requiring separate calibration error
terms for each switch setting.  We use the `TE10` calibration type which
consists of 8-term T error terms and two internal leakage terms.  For the
calibration, we use three standards: short-open, short-match and through.

.. literalinclude:: ../../examples/2x2ab-calibrate.py
   :language: python

Apply the calibration to a device under test.

.. literalinclude:: ../../examples/2x2ab-apply.py
   :language: python


Two Port Reflect Only
---------------------

In this example, we calibrate a two-port VNA for reflection measurements
only, useful in cases where we don't need to make through measurements.
This example demonstrates using multiple solvers simultaneously and
saving more than one calibration in the same calibration file.

.. literalinclude:: ../../examples/2PR-calibrate.py
   :language: python

Apply the calibration to a device under test.

.. literalinclude:: ../../examples/2PR-apply.py
   :language: python


TRL
---

Example of through, reflect, line (TRL) calibration.  We need to know
our reflect and line standards only approximately -- the calibration
process solves for the actual parameters of the standards as well as
the error terms.

.. literalinclude:: ../../examples/TRL-calibrate.py
   :language: python

Apply the calibration to a device under test.

.. literalinclude:: ../../examples/TRL-apply.py
   :language: python


Unknown Through
---------------

Example of unknown through calibration.  For this calibration, we need
three reflect standards on each port and the unknown through between
them.  To reduce calibration steps, we measure two reflect standards
at a time.  We have arbitrarily selected them as short-open, open-match,
and match-short.

.. literalinclude:: ../../examples/UT-calibrate.py
   :language: python

Apply the calibration to a device under test.

.. literalinclude:: ../../examples/UT-apply.py
   :language: python
