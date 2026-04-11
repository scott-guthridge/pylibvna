Calibration Examples
====================

One-Port Reflect Only
---------------------

Calibrate the VNA from measurements of short, open and load standards
using Agilent calibration kit models.

.. literalinclude:: ../../examples/1x1-calibrate.py
   :language: python

Apply the calibration to a device under test.

.. literalinclude:: ../../examples/1x1-apply.py
   :language: python


One-Port Reflect Only with Measured Standards
---------------------------------------------

Calibrate the VNA from measurements of short, open and load standards
using measured values of the standards.

.. literalinclude:: ../../examples/1x1m-calibrate.py
   :language: python

Apply the calibration to a device under test.

.. literalinclude:: ../../examples/1x1m-apply.py
   :language: python


2x1 SOLT
--------

Example of SOLT calibration for a VNA that measures only :math:`S_{11}`
and :math:`S_{21}`.  In this example, we use perfect standards: -1 for
short, 1 for open, 0 for load, and solver.add_through() for a perfect
through.

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
errors in the RF switches without requiring separate calibration error
terms for each switch setting.  We use the `TE10` calibration type,
which consists of 8-term T error terms and two internal leakage terms.
To reduce the number of calibration steps, we measure one-port standards
two at a time: short-open, short-load and through.

.. literalinclude:: ../../examples/2x2ab-calibrate.py
   :language: python

Apply the calibration to a device under test.

.. literalinclude:: ../../examples/2x2ab-apply.py
   :language: python


Two-Port Reflect Only
---------------------

In this example, we calibrate a two-port VNA for reflection measurements
only, useful in cases where we don't need to make through measurements.
This example demonstrates using multiple solvers simultaneously and
saving more than one calibration in the same calibration file.  To keep
the example simple, we're using perfect standards.

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
at a time.  We have arbitrarily selected them as short-open, open-load,
and load-short.

.. literalinclude:: ../../examples/UT-calibrate.py
   :language: python

Apply the calibration to a device under test.

.. literalinclude:: ../../examples/UT-apply.py
   :language: python


Test Fixture Embedding and De-Embedding
---------------------------------------

Calibrate the VNA, embedding a male-to-male coaxial adapter between VNA
and standards.

.. literalinclude:: ../../examples/embed-calibrate.py

Measure the DUT with the same coaxial adapter and de-embed the adapter
from the result.

.. literalinclude:: ../../examples/embed-apply.py


Advanced 16 Term with Measurement Error Modeling
------------------------------------------------

In this example we calibrate a two-port VNA using 16-term
T parameters, taking into account both measurement noise and
connection non-repeatability.  We introduce four main new elements:
the `solver.set_m_error()` method, the `solver.et_tolerance` and
`solver.p_tolerance` attributes, and the `CorrelatedParameter` class.

The `solver.set_m_error()` method takes a frequency vector, a noise floor
vector, and a tracking noise error.  The noise floor vector describes the
complex standard deviation in VNA measurements at each frequency due to
noise in the VNA's detectors when no signal from the DUT is applied.
The optional tracking noise vector describes the complex standard
deviation in VNA measurements at each frequency, proportional to the
amplitude of the received signal due to noise in the VNA's transmitter.
It's assumed that the two noise sources are Gaussian and independent.

The `solver.et_tolerance` and `solver.p_tolerance` attributes set
the change in RMS value of the error terms and unknown parameters,
respectively, sufficiently low to stop iteration.  We set them about
10x smaller than smallest resolution we can realistically expect to
attain given the limits of our calibration accuracy.  Not shown is the
`solver.iteration_limit` attribute that limits the number of iterations
allowed before reporting convergence failure.

In the example, **all** parameters of the calibration standards
are unknown.  The transmission term of the through standard is an
`UnknownParameter` as in the unknown-through example; all others are
`CorrelatedParameters`, known only to be statistically related to
other parameters, in this case, to known perfect short, open, and
load standards.  Calibration must not only solve for the 16 error terms
(really 15, because one is a free variable), but also for the 12 unknown
non-repeatable connection parameters of the calibration standards.

We start the calibration by connecting port 1 to a load standard and
port 2 to a short standard.  From there, we change one connection at
a time, providing two measurements of most parameters.  Finally, we
connect the unknown through, not assuming perfect impedance matches.

.. literalinclude:: ../../examples/T16-EM-calibrate.py
   :language: python

Apply the calibration to a device under test.

.. literalinclude:: ../../examples/T16-EM-apply.py
   :language: python
