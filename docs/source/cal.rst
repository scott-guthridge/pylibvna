libvna.cal Module
=================

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
