.. _hemco-coupling:

##############################
Coupling HEMCO to other models
##############################

This page details technical information useful for developers who wish
to couple :program:`HEMCO` (the "Harmonized" Emissions Component)
emissions component to other models.

The description of :program:`HEMCO` coupling to other models is
available in :cite:`Lin_et_al._2021`, which describes coupling to
`GEOS-Chem Classic <https://geos-chem.readthedocs.io>`_,
`GCHP <https://gchp.readthedocs.io>`_,
`WRF-GC <http://wrf.geos-chem.org>`_,
:program:`CESM2-GC`, and future NOAA models.

========
Overview
========

This work is made possible by a restructuring of :program:`HEMCO`, named HEMCO
3.0. HEMCO 3.0 separates model-specific components such as I/O,
Regridding and the model speciation interface, into modular
components, and isolate the HEMCO emissions Core.

This work is currently being actively worked on by the GEOS-Chem
Support Team and Haipeng Lin (Harvard) as part of coupling GEOS-Chem
with the CESM model.

================
Useful resources
================
- HEMCO Repository: `geoschem/HEMCO <https://github.com/geoschem/HEMCO geoschem/HEMCO>`_ on GitHub.
- Original description paper: :cite:`Keller_et_al._2014`.
- Coupling and HEMCO 3.0 description paper: :cite:`Lin_et_al._2021`.
- `The HEMCO User's Guide <http://wiki.seas.harvard.edu/geos-chem/index.php/The_HEMCO_User%27s_Guide>`_
- `HEMCO versions <http://wiki.seas.harvard.edu/geos-chem/index.php/HEMCO_versions>`_

===========
Terminology
===========

As part of the :program:`HEMCO 3.0` restructuring, "HEMCO" is now divided into
three pieces depending on their function:

- **The HEMCO Core.** Emissions calculations logic, containers, data types, etc.
- **Data Input Layer.** I/O (previously
  :file:`HCOIO_Read/Write_*_Mod`), Regridding
  (:file:`HCO_MESSY_REGRID`, :file:`HCO_INTERP_MOD`), ... This will be
  rearranged into :file:`Regrid/` and :file:`IO/` folders in a future
  version. Right now due to dependencies, some of these files still
  live in the :file:`Core/` folder.
- **Model Interface Layer.** Code that couples :program:`HEMCO` with other
  models. There are common utilities available at
  :file:`Interfaces/HCO_Interface_Common.F90`.

.. note::

   Note that not all code pertinent to model coupling actually lives
   inside of :program:`HEMCO`; this is by design, as data types that
   are external to :program:`HEMCO` (i.e. GEOS-Chem types such as
   ``State_Met``, CESM types such as ``physics_state``, WRF types such
   as ``domain``) must be maintained with the model and not inside
   HEMCO. Some code lives in :file:`Interfaces/`, and some will live
   inside the model.

==================================
Technical Notes (Data Input Layer)
==================================

TBD

=======================================
Technical Notes (Model Interface Layer)
=======================================

HEMCO 3.0 Model Interface Layer Overview
-----------------------------------------

In order to interface :program:`HEMCO` with the target model, there are a few
primary tasks that need to be performed as outlined below.

Data/code that needs to be provided to :program:`HEMCO` based on the
target model's data structures include:

- The clock and time-step of the target model
- List of species and physical properties (molecular weight required;
  other properties such as Henry's law constants are optional, only
  for extensions such as SeaFlux)
- Grid information (``I``, ``J``, ``L`` atmospheric '0-D box'
  dimensions required; if using HEMCO built-in regrid, then specifics
  are needed. See below)

Data/code that needs to be **retrieved from HEMCO** into the target
model's data structures (i.e. state object for constituent
flux/concentrations) include:

- Emissions fluxes (kg/m2/s format) retrieved from HEMCO, aggregated
  per species ID, for current time step
- Other data retrieved from HEMCO (using :code:`HCO_GetPtr` or
  :code:`HCO_EvalFld`)

.. important::

   Avoid calling HEMCO functions directly from outside of a specific
   module designed to interface HEMCO with the model. This is so the
   interface can be updated more easily if subroutines within HEMCO
   such as :code:`HCO_GetPtr` change, and the HEMCO state
   (:code`HcoState`) doesn't need to be passed to everywhere in your
   model that needs to retrieve data from HEMCO. **It is also useful
   so regridding to/from HEMCO can be performed in a centralized
   location, if so needed by the model.** For example, GEOS-Chem wraps
   :code:`HCO_GetPtr` and  :code:`HCO_EvalFld` into its own interface,
   :code:`HCO_GC_GetPtr`, :code:`HCO_GC_EvalFld`, which will
   auto-magically add the :code:`HcoState` argument, in addition to
   handling regridding if necessary.

Things that come out-of-the-box and generally do not require
customization to a specific model:

- Reading configuration file (:file:`HEMCO_Config.rc`), although the
  path needs to be specified
- HEMCO "driver" (run) routines
- Managing HEMCO memory (initializing HEMCO state in ``HcoState``,
  extensions state in ``ExtState``, etc.)

Reading the HEMCO configuration file and defining species list
---------------------------------------------------------------

This is a three-step process. First initialize the configuration
object (:code:`HcoConfig`):

.. code-block:: Fortran

   call ConfigInit(HcoConfig, HMRC, nModelSpecies=nSpc)

You have to register the species first in addition to some other
HcoConfig properties:

.. code-block:: fortran

   HcoConfig%amIRoot   = masterproc
   HcoConfig%MetField  = 'MERRA2'
   HcoConfig%GridRes   = ''
   HcoConfig%nModelSpc = nHcoSpc
   HcoConfig%nModelAdv = nHcoSpc            ! # of adv spc?

   do N = 1, nHcoSpc
      HcoConfig%ModelSpc(N)%ModID   = N ! model id
      HcoConfig%ModelSpc(N)%SpcName = trim(solsym(N))
   enddo

Then open the configuration file in two phases; after phase 1,
initialize the log file on the MPI root process:

.. code-block:: Fortran

   call Config_ReadFile(HcoConfig%amIRoot, HcoConfig, HcoConfigFile, 1, HMRC, IsDryRun=.false.)

   ! Open the log file
   if(masterproc) then
      call HCO_LOGFILE_OPEN(HcoConfig%Err, RC=HMRC)
   endif

   call Config_ReadFile(HcoConfig%amIRoot, HcoConfig, HcoConfigFile, 2, HMRC, IsDryRun=.false.)

.. warning::

   **Note that the species count has to be populated three times.**
   Once above at :code:`ConfigInit`, and twice inside the *initialized
   HEMCO Config object*.

Some species physical properties need to be defined for :program:`HEMCO`
extensions, such as molecular weight and henry's law constants:

.. code-block:: fortran

	!-----------------------------------------------------------------------
	! Register HEMCO species information (HEMCO state object)
	!-----------------------------------------------------------------------
	do N = 1, nHcoSpc
	    HcoState%Spc(N)%ModID         = N               ! model id
	    HcoState%Spc(N)%SpcName       = trim(solsym(N)) ! species name
	    HcoState%Spc(N)%MW_g          = adv_mass(N)     ! mol. weight [g/mol]

	    ! HcoState%Spc(N)%HenryK0 ! [M/atm]
	    ! HcoState%Spc(N)%HenryCR ! [K]
	    ! HcoState%Spc(N)%HenryPKA ! [1]
	enddo

.. note::
	If you are not using HEMCO extensions, only ``ModID``, ``SpcName`` and ``MW_g`` need to be defined.

Defining Grid
-------------

Define atmospheric column numbers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: fortran

   HcoState%NX = my_IM
   HcoState%NY = my_JM
   HcoState%NZ = LM

Define the vertical grid
~~~~~~~~~~~~~~~~~~~~~~~~

There are many ways of defining the vertical discretization. Check
:code:`HCO_VertGrid_Define`.

.. code-block:: fortran

	! Pass Ap, Bp values, units [Pa], [unitless]
	call HCO_VertGrid_Define(HcoState%Config,                &
	                         zGrid = HcoState%Grid%zGrid,    &
	                         nz    = HcoState%NZ,            &
	                         Ap    = Ap,                     &
	                         Bp    = Bp,                     &
	                         RC    = HMRC)

Define horizontal grid parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

   HEMCO **requires HORIZONTAL grid information only if it is using
   internal regridding routines**, i.e. :code:`MAP_A2A` or
   MESSy. Otherwise, this can be filled with dummy information.

.. warning::

   If :program:`HEMCO` internal regridding (:code:`MAP_A2A`) regridding
   routines are used, **only rectilinear grids are supported.**

   This is because :code:`XMid`, :code:`YMid`, ... arrays are
   **1-dimensional** and thus curvilinear coordinates cannot be
   stored. The underlying :code:`MAP_A2A` algorithm **can** handle
   curvilinear; it is just due to the data structure. This will be
   fixed in a future HEMCO version.

.. code-block:: fortran

   ! Point to grid variables
   HcoState%Grid%XMID%Val         => XMid   (my_IS:my_IE  , my_JS:my_JE  )
   HcoState%Grid%YMID%Val         => YMid   (my_IS:my_IE  , my_JS:my_JE  )
   HcoState%Grid%XEdge%Val        => XEdge  (my_IS:my_IE+1, my_JS:my_JE  )
   HcoState%Grid%YEdge%Val        => YEdge  (my_IS:my_IE  , my_JS:my_JE+1)
   HcoState%Grid%YSin%Val         => YSin   (my_IS:my_IE  , my_JS:my_JE+1)
   HcoState%Grid%AREA_M2%Val      => AREA_M2(my_IS:my_IE  , my_JS:my_JE  )

Here we point :program:`HEMCO`'s variables to structures we have
created in the model. Examples in how to create these structures are
available `in the HEMCO-CESM interface
<https://github.com/jimmielin/HEMCO_CESM/blob/development/hco_esmf_grid.F90>`_.

Defining Met Fields for HEMCO Extensions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An example to translate and define meteorological quantities such as
temperature, humidity, etc. is available in the HEMCO-CESM interface.

Running HEMCO
--------------

Prerequisites:

.. code-block:: fortran

   ! HEMCO
   use HCO_Interface_Common,   only: GetHcoVal, GetHcoDiagn
   use HCO_Clock_Mod,          only: HcoClock_Set, HcoClock_Get
   use HCO_Clock_Mod,          only: HcoClock_EmissionsDone
   use HCO_Diagn_Mod,          only: HcoDiagn_AutoUpdate
   use HCO_Driver_Mod,         only: HCO_Run
   use HCO_EmisList_Mod,       only: Hco_GetPtr
   use HCO_FluxArr_Mod,        only: HCO_FluxArrReset
   use HCO_GeoTools_Mod,       only: HCO_CalcVertGrid, HCO_SetPBLm

Update the HEMCO clock
~~~~~~~~~~~~~~~~~~~~~~

Also make sure the time steps are set correctly.
Use from the common utilities:

.. code-block:: fortran

   call HCOClock_Set(HcoState, year, month, day,  &
                     hour, minute, second, IsEmisTime=.true., RC=HMRC)


Reset fluxes for new timestep
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: fortran

   call HCO_FluxArrReset(HcoState, HMRC)

Update vertical grid parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:program:`HEMCO` needs an updated vertical grid at each time step. Data passed
into :code:`HCO_CalcVertGrid` can vary and the definition can be checked
for acceptable parameters.

.. code-block:: fortran

   call HCO_CalcVertGrid(HcoState, PSFC, ZSFC, TK, BXHEIGHT, PEDGE, HMRC)

   call HCO_SetPBLm(HcoState, PBLM=State_HCO_PBLH, &
                    DefVal=1000.0_hp, & ! default value
                    RC=HMRC)

Some dummy setup (advanced)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
To document.

.. code-block:: fortran

   ! Range of species and emission categories.
   ! Set Extension number ExtNr to 0, indicating that the core
   ! module shall be executed.
   HcoState%Options%SpcMin = 1
   HcoState%Options%SpcMax = -1
   HcoState%Options%CatMin = 1
   HcoState%Options%CatMax = -1
   HcoState%Options%ExtNr  = 0

   ! Use temporary array?
   HcoState%Options%FillBuffer = .FALSE.

Run HEMCO driver
~~~~~~~~~~~~~~~~

.. code-block:: fortran

   call HCO_Run( HcoState, 1, HMRC, IsEndStep=.false. )
   call HCO_Run( HcoState, 2, HMRC, IsEndStep=.false. )

Run HEMCO extensions driver
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Necessary only if you are using :program:`HEMCO` extensions.

.. code-block:: fortran

   call HCOX_Run(HcoState, ExtState, HMRC)

Close timestep
~~~~~~~~~~~~~~

.. code-block:: fortran

   !-----------------------------------------------------------------------
   ! Update "autofill" diagnostics.
   ! Update all 'AutoFill' diagnostics. This makes sure that all
   ! diagnostics fields with the 'AutoFill' flag are up-to-date. The
   ! AutoFill flag is specified when creating a diagnostics container
   ! (Diagn_Create).
   !-----------------------------------------------------------------------
   call HcoDiagn_AutoUpdate(HcoState, HMRC)

   !-----------------------------------------------------------------------
   ! Tell HEMCO we are done for this timestep...
   !-----------------------------------------------------------------------
   call HcoClock_EmissionsDone(HcoState%Clock, HMRC)

Retrieving emissions data from HEMCO
--------------------------------------
You can either use the common utilities, where data is retrieved using
:code:`GetHcoValEmis`, or tap into the arrays directly.

For generic data containers, pass the container name like so:

.. code-block:: fortran

   ! For grabbing data from HEMCO Ptrs (uses HEMCO single-precision)
   real(sp), pointer                     :: Ptr2D(:,:)
   real(sp), pointer                     :: Ptr3D(:,:,:)

   logical                               :: FND

   call HCO_GetPtr(HcoState, 'CONTAINER_NAME', Ptr2D, HMRC, FOUND=FND)

Retrieving deposition velocities (depv) from HEMCO
---------------------------------------------------

.. warning::

   **Important:** Note that deposition (sink terms) fluxes are handled
   separately from emissions in HEMCO. This is particularly important
   if you use HEMCO to calculate deposition terms, e.g. the sink term
   in :code:`SeaFlux` (sea-air exchange). The standard in HEMCO is that
   the sink terms are stored as deposition velocities (:code:`depv`,
   unit :code:`1/s`) so HEMCO generally does not need to be aware of
   concentrations.

A thorough discussion of this is in `the HEMCO GitHub issue tracker
<https://github.com/geoschem/HEMCO/issues/72#issuecomment-789409266>`_. The
code to handle deposition velocities from HEMCO is generally as
follows:

.. code-block:: fortran

   !------------------------------------------------------------------
   ! Also add drydep frequencies calculated by HEMCO (e.g. from the
   ! air-sea exchange module) to DFLX.  These values are stored
   ! in 1/s.  They are added in the same manner as the drydep freq values
   ! from drydep_mod.F90.  DFLX will be converted to kg/m2/s later.
   ! (ckeller, 04/01/2014)
   !------------------------------------------------------------------
   CALL GetHcoValDep( NA, I, J, L, found, dep )
   IF ( found ) THEN
      dflx(I,J,NA) = dflx(I,J,NA)                                     &
                   + ( dep * spc(I,J,NA) / (AIRMW / ThisSpc%MW_g)  )
   ENDIF
