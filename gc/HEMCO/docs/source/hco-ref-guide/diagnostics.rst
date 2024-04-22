.. _hco-diag:

#################
HEMCO diagnostics
#################

.. _hco-diag-overview:

========
Overview
========

HEMCO diagnostics are organized in **collections**, with each
collection consisting of a dynamic number of diagnostic fields (aka
**diagnostic containers**). Each collection has a fixed output
frequency (:option:`DiagnFreq`) assigned to it.  All fields within a
collection are written out at the same interval: :option:`Hourly`,
:option:`Daily`, etc.

The contents of a collection (i.e. the diagnostics containers) are
defined at the beginning of a simulation and become continuously updated
and written out during the simulation. A number of attributes attached
to each diagnostic define the properties of a given field and how to
perform field operations such as time averaging, unit conversion, etc.
These attributes include the **field name** (this will also be the netCDF
variable name), the designated field **output units**, the **averaging
method**, and an explicit **unit conversion factor**. The latter three
determine how data is internally stored and returned for output. The
data returned for output is not necessarily in the same units as it is
internally stored.

By default, HEMCO assumes the passed fields are in kg/m2/s, stores
them in kg/m2, and returns the average flux over the designated output
interval in the units assigned to this field (default is
kg/m2/s). This behavior can be avoided by explicitly setting the
averaging method.

**TODO: Find out where these get defined**

Currently supported averaging methods are:

.. option:: instantaneous

   Instantaneous values (recommended method).

.. option:: mean

   Arithmetic mean over the diagnostic interval.

.. option:: sum

   Total sum over the diagnostic interval.

.. option:: cumulsum

   Cumulative sum since simulation start.

Explicitly setting the averaging method will disable automatic unit
conversion and the fields passed to this diagnostic will be stored as
is. The optional unit conversion factor can be set to perform a unit
conversion before returning the field for output.

.. note::

   It is highly recommended to explicitly set the averaging method for
   all fields that are not in kg/m2/s.

.. _hco-diag-builtin:

===============================
Built-in diagnostic collections
===============================

HEMCO has three built-in diagnostic collections (:ref:`Default
<hco-diag-default>`, :ref:`Restart <hco-diag-restart>`, and
:ref:`Manual <hco-diag-manual>`) that are automatically created
on every HEMCO run. These collections are used by HEMCO for internal
data exchange and to write out restart variables. These collections
are 'open', i.e. the user can add additional diagnostic fields to them
if needed. The user can also define new collections (see below).

.. _hco-diag-default:

The Default collection
----------------------

The **Default** collection contains emission diagnostics intended to
be written to disk, e.g. for analysis purposes. All fields of the
default collection are written out at the frequency provided in
setting :option:`DiagnFreq` in the settings section of the HEMCO
configuration file. The name of the corresponding diagnostics files
can be specified via the :code:`DiagnPrefix` setting. The simulation
date at the time of output will be appended to the diagnostics prefix,
e.g. the diagnostics for Aug 1, 2008 will be written as
:file:`HEMCO_Diagnostics.200808010000.nc`. The  datetime can denote
the beginning, middle, or end (default) of the time interval, as
specified by setting :option:`DiagnTimeStamp` (see below).

Several :ref:`options for the default diagnostic collection
<hco-cfg-settings-diagnostics>` can be specified at the top of the
:ref:`HEMCO configuration file <hco-cfg>` file.  Commonly-used options
are :option:`DiagnFile`, :option:`DiagnFreq`, and
:option:`DiagnPrefix`.

.. _hco-diag-configfile:

Configuration file for the Default collection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Adding the following entries to the diagnostic configuration file
(i.e. the same file specified by :option:`DiagnFreq`, commonly called
:file:`HEMCO_Diagn.rc`) will make HEMCO write out total NO and CO
emissions, as well as GFED biomass burning CO emissions (e.g. only
emissions from extension 111):

   .. code-block:: console

      # Name         Spec ExtNr  Cat Hier Dim Unit     LongName
      EmisNO_Total   NO   -1     -1  -1   2   kg/m2/s  NO_emission_flux_from_all_sectors
      EmisCO_Total   CO   -1     -1  -1   2   kg/m2/s  CO_emission_flux_from_all_sectors
      EmisCO_GFED    CO   111    -1  -1   2   kg/m2/s  CO_emission_flux_from_biomass_burning

If you want to just diagnose regional emissions, then you need to
set the diagnostics extension number, category and hierarchy
accordingly. For example, if you want EPA16 emissions for CO over
the USA, then add this line:

   .. code-block:: console

      #Name          Spec ExtNr  Cat Hier Dim Unit     Longname
      EmisCO_EPA16   CO   0      1   50   2   kg/m2/s  CO_emission_flux_from_EPA16_inventory

It is important that you define valid values for all attributes up
to the hierarchy. As soon as you set an attribute to default
(:literal:`-1`),  HEMCO will take the sum up to this attribute. For
example, the following diagnostics would simply return total base
emissions:

   .. code-block:: console

     #Name           Spec ExtNr  Cat Hier Dim Unit     Longname
     EmisCO_EPA16    CO   0      -1  50   2   kg/m2/s  CO_emission_flux_from_EPA16_inventory

.. _hco-diag-restart:

Restart
-------

The output frequency of the **Restart** collection is :literal:`End`,
meaning that its content is only written out at the end of a
simulation. The HEMCO Restart collection primarily consists of a suite
of fields needed by some of the HEMCO extensions for a "warm" HEMCO
restart (e.g. the 10-day running mean temperature, etc.). These fields
are automatically added to the HEMCO restart collection and filled
within the respective extensions. Once archived, fields can be made
available to an extension via the HEMCO configuration file.

.. _hco-diag-manual:

Manual
------

Fields in the **Manual** collection do not become written out to
disk. Rather, they provide a tool to exchange data files within and
outside of HEMCO, e.g. to pass sector-specific emission fluxes from
HEMCO to the atmospheric model.

Some HEMCO extensions automatically create and fill a number of manual
diagnostics. For example, the PARANOX extension (used in `GEOS-Chem
<https://geos-chem.readthedocs.io>`_) stores the O3 and HNO3 loss
fluxes in the manual diagnostics :literal:`PARANOX_O3_DEPOSITION_FLUX`
and :literal:`PARANOX_HNO3_DEPOSITION_FLUX`, respectively.

.. _hco-diag-importing:

===================================================
Importing diagnostic content into an external model
===================================================

The content of the :ref:`Default collection <hco-diag-default>` can
be specified through the HEMCO diagnostics definitions file (specified
by the :option:`DiagnFile` option).

The content of the :ref:`Manual <hco-diag-manual>` and
:ref:`Restart <hco-diag-restart>` collections currently need to
be defined within the model code (e.g. it is hard-coded). This should
be done in high-level routines (at the HEMCO-to-model interface
level).

Module :file:`hco_diagn_mod.F90` (found in :file:`HEMCO/src/Core/`)
provides a suite of routines to define, fill, obtain, etc. diagnostic
fields. Similarly, :file:`hco_restart_mod.F90` (also found in
:file:`HEMCO/src/Core/`) provides routines for managing HEMCO restart
variables.
