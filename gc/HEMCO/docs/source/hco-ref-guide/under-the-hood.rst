.. |br| raw:: html

   <br />

.. _hco-hood:

####################
HEMCO under the hood
####################

.. _hco-hood-overview:

This section provides a short description of the main principles of
HEMCO. More details are provided in the source code, and references to
the corresponding modules is given where appropriate.

========
Overview
========

The HEMCO code can be broken up into three parts: :ref:`core code
<hco-hood-core>`, :ref:`extensions <hco-hood-ext>` and :ref:`interfaces
<hco-hood-int>`.

- The :ref:`core code <hco-hood-core>` consists of all core modules
  that are essential for every HEMCO simulation. |br|
  |br|

- The :ref:`extensions <hco-hood-ext>` are a collection of emission
  parameterizations that can be optionally selected (e.g. dust
  emissions, air-sea exchange, etc.). Most of the extensions require
  meteorological variables (2D or 3D fields) passed from an
  atmospheric model or an external input file to HEMCO.  (See the
  :ref:`hco-ext` section for more information.) |br|
  |br|

- The :ref:`interfaces <hco-hood-int>` are top-level routines that are
  only required in a given model environment (e.g. in stand-alone
  mode or under an ESMF framework). The HEMCO-model interface
  routines are located outside of the HEMCO code structure, calling
  down to the HEMCO driver routines for both the HEMCO core and
  extensions.

HEMCO stores all emission data (:ref:`base emissions <hco-cfg-base>`,
:ref:`scale factors <hco-cfg-scalefac>`, :ref:`masks <hco-cfg-masks>`)
in a generic data structure (a **HEMCO data container**). Input data
read from disk is translated into this data structure by the HEMCO
input/output module (:file:`src/Core/hcoio_dataread_mod.F90`). This
step includes unit conversion and regridding.

.. _hco-hood-data-objects:

==================
HEMCO data objects
==================

All emission data (:ref:`hco-cfg-base`, :ref:`hco-cfg-scalefac`,
:ref:`hco-cfg-masks`) are internally stored in a **data
container**. For each data element of :ref:`the HEMCO configuration
file <hco-cfg>`, a separate data container object is created when
reading the configuration file at the beginning of the simulation. The
data container object is a Fortran derived type that holds information
of one entry of the configuration file. All file data information such
as filename, file variable, time slice options, etc. are stored in a
:code:`FileData` derived type object (defined in
:file:`src/Core/hco_filedata_mod.F90`). This object also holds a
pointer to the data itself. All data is stored as 2 or 3 dimensional
data arrays. HEMCO can keep multiple time slices in memory
simultaneously, e.g. for diurnal scale factors, in which case a vector
of data arrays is created. Data arrays are defined in module
:file:`/src/core/hco_arr_mod.F90`.

Data containers (and as such, emissions data) are accessed through three
different linked lists: :code:`ConfigList`, :code:`ReadList`, and
:code:`EmisList`. These lists all point to the same content (i.e. the
same containers) but ordered in a manner that is most efficient for
the intended purpose:

- For example, :code:`ReadList` contains sub-lists of all containers
  that need to be updated annually, monthly, daily, hourly, or
  never. Thus, if a new month is entered, only a few lists (monthly,
  daily and hourly) have to be scanned and updated instead of going
  through the whole list of data containers. |br|
  |br|

- Similarly, :code:`EmisList` sorts the data containers by model
  species, emission category (:option:`Cat`) and hierarchy
  (:option:`Hier`) . This allows an efficient emission calculation
  since the EmisList has to be scanned only once.

List containers and generic linked list routines are defined  in
:file:`src/Core/hco_datacont_mod.F90`. Specific routines for
:code:`ConfigList`, :code:`ReadList` and :code:`EmisList` are defined
in :file:`src/Core/hco_config_mod.F90`,
:file:`src/Core/hco_readlist_mod.F90`, and
and :file:`src/Core/hco_emislist_mod.F90` respectively.

.. _hco-hood-core:

=========
Core code
=========

HEMCO core consists of all routines and variables required to read,
store, and update data used for emissions calculation. The driver
routines to execute (initialize, run and finalize) a HEMCO core
simulation are (see hco_driver_mod.F90: :code:`HCO_INIT`, :code:`HCO_RUN`,
:code:`HCO_FINAL`). These are also the routines that are called at the
interface level (see :ref:`the HEMCO-to-model interface
<hco-hood-int-to-model>` section).

Each HEMCO simulation is defined by its state object :code:`HcoState`,
which is a derived type that holds all simulation information,
including a list of the defined HEMCO species,  emission grid
information, configuration file name, and additional run options. More
details on the HEMCO state object can be found in
:file:`src/Core/hco_state_mod.F90`. :code:`HcoState` is defined at the
interface level and then passed down to all HEMCO routines

.. _hco-hood-init:

Initialize: HCO_INIT
--------------------

Before running HEMCO, all variables and objects have to be initialized
properly. The initialization of HEMCO occurs in three steps:

#. Read the HEMCO configuration file (subroutine :code:`Config_ReadFile` in
   :file:`src/Core/hco_config_mod.F90`). This writes the content of
   the entire configuration file into buffer, and creates a data
   container for each data item (:ref:`base emission <hco-cfg-base>`
   :ref:`scale factor <hco-cfg-scalefac>`, :ref:`mask <hco-cfg-masks>`)
   in :code:`ConfigList`. |br|
   |br|

#. Initialize :code:`HcoState`. |br|
   |br|

#. Call :code:`HCO_INIT`, passing :code:`HcoState` to it. This
   initializes the HEMCO clock object (see
   :file:`src/Core/hco_clock_mod.F90`) and creates the
   :file:`ReadList` (:file:`src/Core/hco_readlist_mod.F90`). The
   :file:`ReadList` links to the data containers in
   :code:`ConfigList`, but sorted by data update frequency. Data that
   is not used at all (e.g. scale factors that are not used by any
   base emission, or regional emissions that are outside of the
   emission grid). The :code:`EmisList` linked list is only created in
   the run call.

Note that steps 1 and 2 occur at the :ref:`the HEMCO-to-model
interface level <hco-hood-int-to-model>`.

.. _hco-hood-run:

Run: HCO_RUN
------------

This is the main function to run HEMCO. It can be repeated as often as
necessary. Before calling this routine, the internal clock object has to
be updated to the current simulation time (subroutine :code:`HcoClock_Set`
in :file:`src/Core/hco_clock_mod.F90`). :code:`HCO_RUN` performs the
following steps:

#. Updates the time slice index pointers. This is to make sure
   that the correct time slices are used for every data container. For
   example, hourly scale factors can be stored in a data container
   holding 24 individual 2D fields. Module
   :file:`src/Core/hco_tidx_mod.F90` organizes how to properly access
   these fields. |br|
   |br|

#. Read/update the content of the data containers (:code:`ReadList_Read`).
   Checks if there are any fields that need to be read/updated (e.g. if
   this is a new month compared to the previous time step) and updates
   these fields if so by calling the data interface (see
   :ref:`hco-hood-int`). |br|
   |br|

#. Creates/updates the :code:`EmisList` object. Similar to
   :code:`ReadList`, :code:`EmisList` points to the data containers in
   :code:`ConfigList`, but sorted according to species, emission
   hierarchy, emissions category. To optimize emission calculations,
   EmisList already combines :base emission fields that share the same
   species, category, hierarchy, scale factors, and field name
   (without the field name tag, see :ref:`Base Emissions
   <hco-cfg-base>`). |br|
   |br|

#. Calculate core emissions for the current simulation time. This is
   performed by subroutine :code:`hco_calcemis` in
   :file:`src/Core/hco_calc_mod.F90`. This routine walks through
   :code:`EmisList`  and calculates the emissions for every base
   emission field by applying the assigned scale factors to it. The
   (up to 10) container IDs of all scale factors connected to the
   given base emission field (as set in :ref:`the HEMCO configuration
   file <hco-cfg>`) are stored in the data container variable
   :code:`ScalIDs`. A container ID index list is used to efficiently
   retrieve a pointer to  each of those containers (see
   :code:`cIDList` in :file:`src/Core/hco_datacont_mod.F90`).

.. _hco-hood-final:

Finalize: HCO_FINAL
-------------------

This routine cleans up all internal lists, variables, and objects. This
does not clean up the HEMCO state object, which is removed at the
interface level.

.. _hco-hood-ext:

==========
Extensions
==========

:ref:`hco-ext` are used to calculate emissions based on
meteorological input variables and/or non-linear
parameterizations. Each extension is provided in a separate Fortran
module. Each module must contain a public subroutine to initialize,
run and finalize the extension. Emissions calculated in the extensions
are added to the HEMCO emission array using subroutine
:code:`HCO_Emis_Add` in :file:`src/Core/HCO_FluxArr_mod.F90`.

Meteorological input data is passed to the individual extension
routines through the extension state object ExtState, which provides a
pointer slot for all met fields used by any of the extension (see
:file:`src/Extensions/hcox_state_mod.F90`). These pointers must be
assigned at the interface level (see :ref:`the HEMCO-model interface
section <hco-hood-int>`).

In analogy to the core module, the three main routines for the
extensions are (in :file:`src/Extensions/hcox_driver_mod.F90`):

-  :code:`HCOX_Init`
-  :code:`HCOX_Run`
-  :code:`HCOX_Final`

These subroutines invoke the corresponding calls of all (enabled)
:ref:`extensions <hco-ext>` and must be called at the interface level
(after the core routines).

Extension settings (as specified in the configuration file, see also
:ref:`hco-cfg-ext-switches`) areautomatically read by HEMCO. For any
given extension, routines :code:`GetExtNr` and :code:`GetExtOpt` can
be used to obtain the extension number (:option:`ExtNr`) and desired
setting value, respectively (see
:file:`src/Core/HCO_ExtList_Mod.F90`). Routine :code:`HCO_GetExtHcoID`
should be used to extract the HEMCO species IDs of all species
registered for this extension.

Gridded data associated to an extension (i.e. listed in section
extension data of the configuration file) is automatically added to
the :code:`EmisList`, but ignored by the HEMCO core module during emissions
calculation. Pointers to these data arrays can be obtained through
routine EmisList_GetDataArr in HCO_EmisList_Mod.F90. Note that
this routine identifies the array based on its container name. It is
therefore important that the container name set in the configuration
file matches the names used by this routine!

.. _hco-hood-int:

==========
Interfaces
==========

.. _hco-hood-int-to-model:

HEMCO-to-model interface
-------------------------

.. note::

   For additional information about coupling HEMCO to other models,
   please see our :ref:`hemco-coupling` chapter.

The interface provides the link between HEMCO and the model environment.
This may be a sophisticated Earth System model or a simple environment
that allows the user to run HEMCO in standalone mode. The standalone
interface is provided along with the HEMCO distribution
(:file:`src/Interfaces/hcoi_standalone_mod.F90`). The
HEMCO-to-GEOS-Chem model interface is included in the GEOS-Chem source
code (:file:`GeosCore/hcoi_gc_main_mod.F90`). HEMCO has also been
successfully employed as a stand-alone gridded component within an
ESMF environment. Please contact Christoph Keller for more information
on the ESMF implementation.

The interface routines provide HEMCO with all the necessary information
to perform the emission calculation. This includes the following tasks:

.. _hco-hood-int-to-model-init:

Initialization:
~~~~~~~~~~~~~~~

-  Read the :ref:`configuration file <hco-cfg>` (:code:`Config_ReadFile` in
   :file:`src/Core/hco_config_mod.F90`). |br|
   |br|

-  Initialize :code:`HcoState object` (:code:`HcoState_Init` in
   :file:`src/Core/hco_state_mod.F90`). |br|
   |br|

-  Define the emission grid. Grid definitions are stored in
   :code:`HcoState%Grid`. The emission grid is defined by its
   horizontal mid points and edges (all 2D fields), the hybrid sigma
   coordinate edges (3D), the grid box areas (2D), and the grid box
   heights. The latter is only used by some extensions
   (:option:`DustDead`, :option:`LightNOx`') and may be left undefined
   if those are not used. |br|
   |br|

-  Define emission species. Species definitions are stored in vector
   :code:`HcoState%Spc(:)` (one entry per species). For each species, the
   following parameter are required:

   #. HEMCO species ID: unique integer index for species identification.
      For internal use only.
   #. Model species ID: the integer index assigned to this species by
      the employed model.
   #. Species name
   #. Species molecular weight in g/mol.
   #. Emitted species molecular weight in g/mol. This value can be
      different to the species molecular weight if species are emitted
      on a molecular basis, e.g. in mass carbon (in which case the
      emitted molecular weight becomes 12 g/mol).
   #. Molecular ratio: molecules of emitted species per molecules of
      species. For example, if C3H8 is emitted as kg C, the molecular
      ratio becomes 3.
   #. K0: Liquid over gas Henry constant in M/atm.
   #. CR: Temperature dependency of K0 in K.
   #. pKa: The species pKa, used for correction of the Henry constant.

   The molecular weight - together with the molecular ratio - determine the
   mass scaling factors used for unit conversion in hco_unit_mod.F90. The
   Henry coefficients are only used by the air-sea exchange extension (and
   only for the specified species) and may be left undefined for other
   species and/or if the extension is not used.

-  Define simulation time steps. The emission, chemical and dynamic time
   steps can be defined separately. |br|
   |br|

-  Initialize HEMCO core (:code:`HCO_Init` in
   :file:`src/Core/hco_driver_mod.F90`) |br|
   |br|

-  Initialize HEMCO extensions (code:`HCOX_Init` in
   :file:`src/Core/hcox_driver_mod.F90`)


.. _hco-hood-int-to-model-run:

Run:
~~~~

-  Set current time (:code:`HcoClock_Set` in
   :file:`src/Core/hco_clock_mod.F90`) |br|
   |br|

-  Reset all emission and deposition values (:code:`HCO_FluxArrReset` in
   :file:`src/Core/hco_fluxarr_mod.F90`) |br|
   |br|

-  Run HEMCO core to calculate emissions (:code:`HCO_Run` in
   :file:`src/Core/hco_driver_mod.F90`) |br|
   |br|

-  Link the used meteorology field objects of :code:`ExtState` to
   desired data arrays (this step may also be done during
   initialization) |br|
   |br|

-  Run :ref:`hco-ext` to add extensions emissions (:code:`HCOX_Run` in
   :file:`src/Core/hcox_driver_mod.F90`) |br|
   |br|

-  Export HEMCO emissions into desired environment

.. _hco-hood-int-to-model-final:

Finalization:
~~~~~~~~~~~~~

-  Finalize HEMCO extensions and extension state object ExtState
   (HCOX_Final in hcox_driver_mod.F90).
-  Finalize HEMCO core (HCO_Final in hco_driver_mod.F90).
-  Clean up HEMCO state object HcoState (HcoState_Final in
   hco_state_mod.F90).

.. _hco-hood-int-data:

Data interface (reading and regridding)
---------------------------------------

The data interface (in :file:`src/Core/hcoi_dataread_mod.F90`)
organizes reading, unit conversion, and remapping of data from source
files. Its public routine :file:`HCOI_DataRead` is only called by subroutine
:code:`ReadList_Fill` in :file:`src/Core/hco_readlist_mod.F90`. Data
processing is performed in three steps:

#. Read data from file using the source file information (file name,
   source variable, desired time stamp) provided in the configuration
   file. |br|
   |br|

#. Convert unit to HEMCO units based on the unit attribute read from
   disk and the srcUnit attribute set in the configuration file. See
   :ref:`hco-filefmt` for more information. |br|
   |br|

#. Remap original data onto the HEMCO emission grid. The grid dimensions
   of the input field are determined from the source file. If only
   horizontal regridding is required, e.g. for 2D data or if the number
   of vertical levels of the input data is equal to the number of
   vertical levels of the HEMCO grid, the horizontal interpolation
   routine used by GEOS-Chem is invoked. If vertical regridding is
   required or to interpolate index-based values (e.g. discrete integer
   values), the NcRegrid tool described in `Joeckel
   (2006) <#References>`__ is used.

.. _hco-hood-multi:

===============================
Run multiple instances of HEMCO
===============================

It is possible to run multiple instances of HEMCO at the same
time. These instances can operate on different grids, use different
configuration files, etc. This is made possible by wrapping all
information of a HEMCO simulation into a :code:`HCO_State` derived
type object (defined in
:file:`src/Core/hco_state_mod.F90`). Similarly, all emission extension
information is included in an :code:`Ext_State` derived type (in
:file:`src/Extensions/hcox_state_mod.F90`). These two objects together
fully define the HEMCO setup and are being passed to the top level
HEMCO routines (INIT/RUN/FINALIZE), e.g.:

.. code-block:: fortran

   CALL HCO_Run( am_I_Root, HcoState, Phase, RC )
   # ...etc ...
   CALL HCOX_Run( am_I_Root, HcoState, ExtState, RC )

To run more than one HEMCO instance in parallel, one need to define
multiple HcoState instances and then call each of these separately,
e.g.:

.. code-block:: fortran

   CALL HCO_Run( am_I_Root, HcoStateA, Phase, RC )
   CALL HCO_Run( am_I_Root, HcoStateB, Phase, RC )
   # ... etc ...

The HEMCO state objects also carry the 3D emission arrays, and when
using multiple instances one needs to ensure that these arrays are
properly connected to the 'emission end user', e.g. PBL mixing routine,
etc. In the GEOS-Chem implementation of HEMCO, the module
hco_interface_mod.F90 (in GeosCore) provides the interface between
HEMCO and GEOS-Chem: it is the owner of the HcoState and ExtState
object, and contains a number of wrapper routines to exchange
information between HEMCO and GEOS-Chem. In the GEOS model, the
standalone HEMCO component uses a linked list that can carry a dynamic
number of HEMCO instances, and then loops over the linked list to
perform all model operations (init,run,finalize) on all members of the
linked list.

.. important::

   Several HEMCO extensions still use global arrays and currently
   cannot be used in multi-instance simulations. As of 8/29/2018, the
   following extensions are likely to cause problems in multi-instance
   simulations: Ginoux dust emissions, FINN biomass burning, GFED
   biomass burning, Iodine emissions, PARANOx ship emissions, sea flux
   emissions, sea salt emissions.
