libvna.cal Module
=================

.. automodule:: libvna.cal

Calset
------

.. autoclass:: libvna.cal.Calset(filename=None) -> Calset
   :no-members:

**Quick Reference** - Calset methods by category:

- Calibration access: :attr:`~libvna.cal.Calset.calibrations`
- Standards: :meth:`~libvna.cal.Calset.short_standard`, :meth:`~libvna.cal.Calset.open_standard`, :meth:`~libvna.cal.Calset.load_standard`, :meth:`~libvna.cal.Calset.through_standard`, :meth:`~libvna.cal.Calset.data_standard`
- Parameters: :meth:`~libvna.cal.Calset.scalar_parameter`, :meth:`~libvna.cal.Calset.vector_parameter`, :meth:`~libvna.cal.Calset.unknown_parameter`, :meth:`~libvna.cal.Calset.correlated_parameter`, :meth:`~libvna.cal.Calset.parameter`, :meth:`~libvna.cal.Calset.parameter_matrix`
- Solving: :meth:`~libvna.cal.Calset.solver`
- Persistence: :attr:`~libvna.cal.Calset.properties`, :meth:`~libvna.cal.Calset.save`

.. autoattribute:: libvna.cal.Calset.calibrations

Calibration Objects
^^^^^^^^^^^^^^^^^^^

.. autoclass:: libvna.cal.Calibration()
   :members: name, ctype, rows, columns, frequencies, frequency_vector,
        z0, properties, apply

Calibration Standards
^^^^^^^^^^^^^^^^^^^^^

The following methods create calibration kit and data-based standards that
can be used with the Solver class below.

.. py:method:: libvna.cal.Calset.short_standard(offset_delay: float = 0.0, offset_loss: float = 0.0, offset_z0: float = 50.0, fmin: float = 0.0, fmax: float = INFINITY, traditional: bool = False, L: Optional[List[float]] = None) -> ShortStandard

.. py:method:: libvna.cal.Calset.open_standard(offset_delay: float = 0.0, offset_loss: float = 0.0, offset_z0: float = 50.0, fmin: float = 0.0, fmax: float = INFINITY, traditional: bool = False, C: Optional[List[float]] = None) -> OpenStandard

.. py:method:: libvna.cal.Calset.load_standard(offset_delay: float = 0.0, offset_loss: float = 0.0, offset_z0: float = 50.0, fmin: float = 0.0, fmax: float = INFINITY, traditional: bool = False, Zl: complex = 50.0) -> LoadStandard

.. py:method:: libvna.cal.Calset.through_standard(offset_delay: float = 0.0, offset_loss: float = 0.0, offset_z0: float = 50.0, fmin: float = 0.0, fmax: float = INFINITY, traditional: bool = False) -> ThroughStandard

.. py:method:: libvna.cal.Calset.data_standard(npdata: NPData) -> DataStandard

These standards have defined characteristic impedances based on the
*offset_z0* parameter or reference impedance given with the data standard.
If used with a VNA port with a different reference impedance, they
automatically re-normlize the reference impedance.

Parameters
^^^^^^^^^^

The library also provides a lower-level way to specify standards,
where each element of the S-parameter matrix of the standard is given
as a Parameter object.  Parameter objects include scalar type that have
a constant value for all frequencies, such as -1.0 for short, 1.0 for
open, and 0.0 for match, vector type given as a vector of frequency
points and a vector of values, unknown type that are only approximately
known such as the R and L parameters in TRL calibration, and correlated
type that are known only to be statistically related to other parameters.
The library solves for unknown and correlated parameters at the same time
it solves for the error terms.  These Parameter objects do not have a
defined reference impedance and its the caller's responsibilty to ensure
that parameters are meaningful for the VNA ports on which they are used.

.. automethod:: libvna.cal.Calset.scalar_parameter
.. automethod:: libvna.cal.Calset.vector_parameter
.. automethod:: libvna.cal.Calset.unknown_parameter
.. automethod:: libvna.cal.Calset.correlated_parameter
.. automethod:: libvna.cal.Calset.parameter
.. automethod:: libvna.cal.Calset.parameter_matrix

Standard and Parameter Classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: libvna.cal.Parameter
   :members: eval, embed, deembed

.. autoclass:: libvna.cal.ParameterMatrix()
   :show-inheritance:
   :members: eval, to_npdata, embed, deembed

.. py:class:: libvna.cal.ScalarParameter(calset: Calset, value) -> ScalarParameter

    Bases: :class:`Parameter`

.. py:class:: libvna.cal.VectorParameter(calset: Calset, frequency_vector, value_vector) -> VectorParameter

    Bases: :class:`Parameter`

.. py:class:: libvna.cal.UnknownParameter(calset: Calset, initial_guess) -> UnknownParameter

    Bases: :class:`Parameter`

.. py:class:: libvna.cal.CorrelatedParameter(calset: Calset, other, frequency_vector, sigma_vector) -> CorrelatedParameter

    Bases: :class:`Parameter`

.. py:class:: libvna.cal.ShortStandard(calset: Calset, offset_delay: float = 0.0, offset_loss: float = 0.0, offset_z0: float = 50.0, fmin: float = 0.0, fmax: float = INFINITY, traditional: bool = False, L: Optional[List[float]] =None) -> ShortStandard

    Bases: :class:`Parameter`

.. py:class:: libvna.cal.OpenStandard(calset, offset_delay=0.0, offset_loss=0.0, offset_z0=50.0, fmin=0.0, fmax=INFINITY, traditional=False, C=None) -> OpenStandard

    Bases: :class:`Parameter`

.. py:class:: libvna.cal.LoadStandard(calset, offset_delay=0.0, offset_loss=0.0, offset_z0=50.0, fmin=0.0, fmax=INFINITY, traditional=False, Zl=50.0) -> LoadStandard

    Bases: :class:`Parameter`

.. py:class:: libvna.cal.ThroughStandard(calset, offset_delay=0.0, offset_loss=0.0, offset_z0=50.0, fmin=0.0, fmax=INFINITY, traditional=False) -> ThroughStandard

    Bases: :class:`ParameterMatrix`

.. py:class:: libvna.cal.DataStandard(calset: Calset, npdata: NPData) -> DataStandard

    Bases: :class:`ParameterMatrix`

Solving a Calibration
^^^^^^^^^^^^^^^^^^^^^

The :class:`Solver` class takes measurements of calibration standards,
solves for the error terms and adds the new calibration to the Calset.

.. automethod:: libvna.cal.Calset.solver

.. autoclass:: libvna.cal.Solver
   :members: add_single_reflect, add_double_reflect, add_through, add_line,
        add_mapped_matrix, set_m_error, et_tolerance, p_tolerance,
        iteration_limit, pvalue_limit, solve, add_to_calset

Embedding and De-embedding Network Parameter Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automethod:: libvna.cal.Calset.embed_npdata
.. automethod:: libvna.cal.Calset.deembed_npdata

Metadata and Persistence
^^^^^^^^^^^^^^^^^^^^^^^^

Calset.properties and Calibration.properties can be assigned arbitrary
trees of nested lists, dictionaries, scalars and None values making
it possible to save use-defined metadata with the calibrations.
The :func:`save` function saves the Calset object to a file.

.. autoattribute:: libvna.cal.Calset.properties
.. automethod:: libvna.cal.Calset.save
