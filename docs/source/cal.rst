libvna.cal: VNA Calibration
===========================

.. toctree::
   :maxdepth: 2

.. automodule:: libvna.cal

.. autoclass:: libvna.cal.Calset
   :members: calibrations, index, properties, save

.. autoclass:: libvna.cal.Calibration
   :members: name, ctype, rows, columns, frequencies, frequency_vector,
        z0, properties, apply

.. autoclass:: libvna.cal.Parameter
   :members:

.. autoclass:: libvna.cal.ScalarParameter
   :members:
   :show-inheritance:

.. autoclass:: libvna.cal.VectorParameter
   :members:
   :show-inheritance:

.. autoclass:: libvna.cal.UnknownParameter
   :members:
   :show-inheritance:

.. autoclass:: libvna.cal.CorrelatedParameter
   :members:
   :show-inheritance:

.. autoclass:: libvna.cal.Solver
   :members: add_single_reflect, add_double_reflect, add_through, add_line,
        add_mapped_matrix, set_m_error, et_tolerance, p_tolerance,
        iteration_limit, pvalue_limit, solve, add_to_calset

SOLT Calibration Example
------------------------

Example of SOLT calibration for a VNA that measures only :math:`S_{11}`
and :math:`S_{21}`

.. code:: python

    from libvna.cal import Calset, CalType, Solver

    f_vector = **get_calibration_frequencies**

    calset = Calset()
    solver = Solver(calset, CalType.E12, rows=2, columns=1,
                    frequency_vector=f_vector)

    measurement = **measure_short_standard**
    solver.add_single_reflect(a=None, b=measurement, s11=-1.0)

    measurement = **measure_open_standard**
    solver.add_single_reflect(a=None, b=measurement, s11=1.0)

    measurement = **measure_load_standard**
    solver.add_single_reflect(a=None, b=measurement, s11=0.0)

    measurement = **measure_through_standard**
    solver.add_through(a=None, b=measurement, port1=1, port2=2)

    solver.solve()
    solver.add_to_calset("cal_2x1")
    calset.save("MySOLT.vnacal")

Example of applying the calibration to a device under test, switching
the probes to get :math:`S_{11}`, :math:`S_{12}`, :math:`S_{21}` and
:math:`S_{22}`, and saving the result to a Touchstone v1 file:

.. code:: python

    from libvna.cal import Calset, Calibration
    from libvna.data import Data
    from numpy import flipud, hstack

    calset = Calset("MySOLT.vnacal")
    calibration = calset.calibrations["cal_2x1"]

    f_vector = **get_measurement_frequencies**

    m1 = **make 2x1 forward measurement**
    m2 = **make 2x1 reverse measurement**
    measurement = hstack((m1, flipud(m2)))

    corrected = calibration.apply(f_vector, None, measurement)
    corrected.format = "SdB"
    corrected.save("MyMeasurement.s2p")
