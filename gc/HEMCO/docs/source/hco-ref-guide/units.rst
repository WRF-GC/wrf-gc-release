.. _hco-units:

##############
Units in HEMCO
##############

.. _hco-units-overview:

========
Overview
========

.. attention::

   We recommend that you provide explicit scale factors for unit
   conversions in :ref:`the HEMCO configuration file <hco-cfg>`.  This
   will avoid some :ref:`known issues <hco-known-bugs>` with unit
   conversions that were recently discovered.

HEMCO classifies all data fields as fluxes, concentrations, or unitless
data. Data are internally stored in HEMCO standard units of
:literal:`[kg emitted species/m2/s]` for fluxes, and :literal:`[kg
emitted species/m3]` for concentrations. No unit conversion is
performed on unitless data.

The classification of a data field depends on the units attribute in the
netCDF file, the :option:`SrcUnit` attribute in :ref:`the HEMCO
configuration file <hco-cfg>`, and the unit tolerance setting in the
HEMCO configuration file (see below). In general, the original units
of the input data is determined based on the units attribute on the
netCDF file, and data is converted to HEMCO units accordingly. The
mass conversion factor is determined based on the species assigned to
the given field throuh attribute :option:`Species` in the HEMCO
configuration file. It depends on the species molecular weight (MW),
the MW of the emitted species, and the molecular ratio (molecules of
emitted species per molecules of species). If the input data is found
to be in non-standard units (e.g. :literal:`L` instead of
:literal:`m3`, :literal:`g` instead of :literal:`kg`, etc.), HEMCO
will attempt to convert to standard units.
This feature is not fully tested yet, and it is recommended to provide
input data in standard units wherever possible.

.. _hco-units-srcunit:

=================
SrcUnit attribute
=================

The :option:`SrcUnit` attribute in :ref:`the HEMCO configuration file
<hco-cfg>` gives the user some control on unit conversion.

If :option:`SrcUnit` is set to :literal:`1`, data are treated as
unitless irrespective of the units attribute on the file. This option
works on all fields only if unit tolerance is relaxed to :literal:`2`
(for unit tolerance of :literal:`1`, the input data must be in one of
the units recognized by  HEMCO as :literal:`unitless`).

If :option:`SrcUnit` is set to :literal:`count`, the input data is
assumed to represent index-based scalar fields (e.g. land types). No
unit conversion is performed on these data and regridding will
preserve the absolute values.

Special attention needs to be paid to species that are emitted in
quantities other than mass of species, e.g. :literal:`kg C`. For these
species, the species MW differs from the emitted species MW, and the
molecular ratio determines how many molecules are being emitted per
molecule species. By default, HEMCO will attempt to convert all input
data to kg emitted species. If a species is emitted as
:literal:`kgC/m2/s` and the input data is in kg/m2/s, the mass will be
adjusted based on the emitted MW, species MW, and the ratio
emitted MW / species MW. Only input data that is already in
:literal:`kgC/m2/s` will not be converted. This behavior can be
avoided by explicitly set the :option:`SrcUnit` to the same unit as on
the input file. In this case, HEMCO will not convert between species MW
and emitted MW. This is useful for cases where the input data does not
contain data of the actual species, e.g. if VOC emissions are calculated
by scaling CO emissions (see examples below).

.. _hco-units-unit-tolerance:

======================
Unit tolerance setting
======================

The unit tolerance setting (see the :ref:`Settings <hco-cfg-settings>`
section of :ref:`the HEMCO configuration file <hco-cfg>` indicates the
tolerance of HEMCO if discrepancies are found between the units found in
the input file and attribute :option:`SrcUnit` of the configuration
file.

- If the unit tolerance is set to :literal:`0`, HEMCO stops with an
  error if the :option:`SrcUnit` attribute does not exactly match with the units
  attribute found in the input data.

- Unit tolerance of :literal:`1` enables the default behavior.

- Unit tolerance of :literal:`2` will take the :option:`SrcUnit`
  attribute as the data input unit, regardless netCDF units attribute.

.. _hco-units-unitless:

=============
Unitless data
=============

The following units are currently recognized as 'unitless' by HEMCO

- :literal:`1`
- :literal:`count`
- :literal:`unitless`
- :literal:`fraction`
- :literal:`factor`
- :literal:`scale`
- :literal:`hours`
- :literal:`v/v`
- :literal:`v/v/s`
- :literal:`s-1`
- :literal:`m2/m2`
- :literal:`kg/kg`
- :literal:`K`
- :literal:`W/m2`
- :literal:`pptv`
- :literal:`ppt`
- :literal:`ppbv`
- :literal:`ppb`
- :literal:`ppmv`
- :literal:`ppm`
- :literal:`ms-1`
- :literal:`m`
- :literal:`cm2cm-2`
- :literal:`dobsons`
- :literal:`dobsons/day`
- :literal:`hPa`
- :literal:`Pa`

.. _hco-unit-example:

===================
Examples with units
===================

.. attention::

   We recommend that you provide explicit scale factors for unit
   conversions in :ref:`the HEMCO configuration file <hco-cfg>`.  This
   will avoid some :ref:`known issues <hco-known-bugs>` with unit
   conversions that were recently discovered.

File :file:`file1.nc` contains field :literal:`DATA` in units of
:literal:`kg/m2/s`. It shall be applied to species acetone
(:literal:`ACET`), which is emitted as :literal:`kg C`. The species
molecular weight of ACET  is :literal:`58`, the emitted molecular
weight is :literal:`12` (i.e. that of carbon), and the molecular ratio
is :literal:`3` (3 molecules of carbon per molecule of acetone).

The following entry in the HEMCO configuration file will interpret the
input data as :literal:`kg acetone/m2/s`, and convert it to
:literal:`kg C/m2/s` using a scale factor of :literal:`0.62` (= 12/58*3):

.. code-block:: kconfig

  #--> data is converted from kg acetone/m2/s to kgC/m2/s
  0 ACET  /path/to/file1.nc  DATA 2000/1/1/0 C xy kgC/m2/s ACET - 1 1

The following entry will avoid the unit conversion from kg to kgC:

.. code-block:: kconfig

   #--> data is kept in kg species/m2/s
   0 ACET  /path/to/file1.nc  DATA 2000/1/1/0 C xy kg/m2/s ACET - 1 1

Note that the opposite does not work: If :file:`file2.nc` contains
data in units of :file:`kgC/m2/s`, it is not possible to convert to kg
species/m2/s and the following two entries have the same effect:

.. code-block:: kconfig

    #--> data is converted from kgC/m2/s to kg emitted species/m2/s,
    #    which is also kgC/m2/s``
   0 ACET  /path/to/file2.nc  DATA 2000/1/1/0 C xy kg/m2/s  ACET - 1 1

   #--> data is kept in kgC/m2/s
   0 ACET  /path/to/file2.nc  DATA 2000/1/1/0 C xy kgC/m2/s ACET - 1 1

However, if one wants to use file2 for a species not emitted as kg
carbon, say CO, the source unit attribute matters!

.. code-block:: kconfig

    #--> data is converted from kgC/m2/s to kg CO/m2/s
   0 ACETasCO  /path/to/file2.nc  DATA 2000/1/1/0 C xy kg/m2/s  CO - 1 1

   #--> data is kept in kgC/m2/s
   0 ACETasCO  /path/to/file2.nc  DATA 2000/1/1/0 C xy kgC/m2/s CO - 1 1

.. _hco-unit-tips:

================
Tips for testing
================

The unit factor applied by HEMCO is written into the HEMCO log file if
:option:`Verbose` is set to 2 or higher.
