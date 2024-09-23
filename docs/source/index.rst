.. libvna documentation master file, created by
   sphinx-quickstart on Sat May  6 19:22:39 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for the libvna Package
====================================
**libvna** is a collection of modules that perform advanced calibration,
parameter conversion, and management of network parameter data for vector
network analyzers (VNAs).

Libvna uses a general calibration solver that handles the classic
SOLT, TRD, etc. techniques, the 8-term T parameter techniques: TRL,
LRL, TRM, LRM, TXYZ, LXYZ, TXYX, LXYX, LRRM, UXYZ (a.k.a. unknown
through), etc., and the 16-term T parameter techniques for any number
of ports. Calibration can model both measurement and connection
non-repeatability errors. The library includes 72 network parameter
conversions plus conversion to impedance looking into each DUT
port. It supports Touchstone 1 and 2 file formats and a more general
space-separated-field network parameter data (.npd) file format.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   cal-menu
   data
   conv


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
