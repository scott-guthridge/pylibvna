libvna.data: Network Parameter Data
===================================

.. toctree::
   :maxdepth: 2

.. autoclass:: libvna.data.NPData

Managing Type and Dimensions
----------------------------
.. automethod:: libvna.data.NPData.init
.. automethod:: libvna.data.NPData.resize
.. automethod:: libvna.data.NPData.convert
.. autoattribute:: libvna.data.NPData.ptype
.. autoattribute:: libvna.data.NPData.ptype_name
.. autoattribute:: libvna.data.NPData.rows
.. autoattribute:: libvna.data.NPData.columns
.. autoattribute:: libvna.data.NPData.frequencies

The Frequency Vector
--------------------
.. autoattribute:: libvna.data.NPData.frequency_vector
.. automethod:: libvna.data.NPData.add_frequency

The Data Array
--------------
.. autoattribute:: libvna.data.NPData.data_array

Reference Impedances
--------------------
.. autoattribute:: libvna.data.NPData.z0_vector
.. autoattribute:: libvna.data.NPData.fz0_array
.. autoattribute:: libvna.data.NPData.has_fz0

Loading and Saving
------------------
.. automethod:: libvna.data.NPData.load
.. automethod:: libvna.data.NPData.save
.. automethod:: libvna.data.NPData.fload
.. automethod:: libvna.data.NPData.fsave
.. automethod:: libvna.data.NPData.cksave
.. autoattribute:: libvna.data.NPData.format
.. autoattribute:: libvna.data.NPData.filetype
.. autoattribute:: libvna.data.NPData.fprecision
.. autoattribute:: libvna.data.NPData.dprecision
