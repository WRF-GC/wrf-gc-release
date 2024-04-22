.. _hco-ext:

################
HEMCO extensions
################

========
Overview
========

Emission inventories sometimes include dynamic source types and
nonlinear scale factors that have functional dependencies on local
environmental variables such as wind speed or temperature, which are
best calculated online during execution of the model. HEMCO includes a
suite of additional modules (extensions) that perform online emission
calculations for a variety of sources (see list below). Extensions are
enabled in section :ref:`Extension Switches <hco-cfg-ext-switches>`
of :ref:`the HEMCO configuration file <hco-cfg>`.

.. _hco-ext-list:

==================
List of extensions
==================

The full list of available extensions is given below. Extensions can be
selected individually in the :ref:`Extension Switches
<hco-cfg-ext-switches>` section of the :ref:`hco-cfg`, as can the species to
be considered.

.. option:: DustAlk

   - **Species**: DSTAL1, DSTAL2, DSTAL3, DSTAL4
   - **Reference**: Fairlie et al (check)

.. option:: DustDead

   Emissions of mineral dust from the DEAD dust mobilization model.

   - **Species**: DST1, DST2, DST3, DST4
   - **Reference**: :cite:t:`Zender_et_al._2003`

.. option:: DustGinoux

   Emissions of mineral dust from the P. Ginoux dust mobilization model.

   - **Species**: DST1, DST2, DST3, DST4
   - **Reference**: :cite:t:`Ginoux_et_al._2001`

.. option:: FINN

   Biomass burning emissions from the FINN model.

   - **Species**: NO, CO, ALK4, ACET, MEK, ALD2, PRPE, C2H2, C2H4,
     C3H8, CH2O, C2H6, SO2, NH3, BCPI, BCPO, OCPI, OCPO, GLYC, HAC,
     SOAP
   - **Reference**: :cite:t:`Wiedinmyer_et_al._2014`

.. option:: GC_Rn-Pb-Be

   Emissions of radionuclide species as used in the `GEOS-Chem
   <https://geos-chem.readthedocs.io>`_ model.

   - **Species**: Rn222, Be7, Be7Strat, Be10, Be10Strat

   .. option:: ZHANG_Rn222

      If :option:`ZHANG_Rn222` is :literal:`on`, then Rn222 emissions
      will be computed according to :cite:t:`Zhang_et_al._2021`.

      If :option:`ZHANG_Rn222` is :literal:`off`, then Rn222 emissions
      will be computed according to :cite:t:`Jacob_et_al._1997`.

.. option:: GFED

   Biomass burning emissions from the GFED model.

   - **Version**: GFED3 and GFED4 are available.
   - **Species**: NO, CO, ALK4, ACET, MEK, ALD2, PRPE, C2H2, C2H4, C3H8, CH2O
     C2H6, SO2, NH3, BCPO, BCPI, OCPO, OCPI, POG1, POG2, MTPA, BENZ, TOLU, XYLE
     NAP, EOH, MOH, SOAP,
   - **Reference**: :cite:t:`van_der_Werf_et_al._2010`

.. option:: Inorg_Iodine

   - **Species**:  HOI, I2
   - **Reference**:  TBD

.. option:: LightNOx

   Emissions of NOx from lightning.

   - **Species**: NO
   - **Species**: :cite:`Murray_et_al._2012`

.. option:: MEGAN

   Biogenic VOC emissions.

   - **Version**: 2.1
   - **Species:** ISOP, ACET, PRPE, C2H4, ALD2, CO, OCPI, MONX, MTPA, MTPO,
     LIMO, SESQ
   - **Reference:** :cite:t:`Guenther_et_al._2012`

.. option:: PARANOx

   Plume model for ship emissions.

   - **Species**: NO, NO2, O3, HNO3
   - **Reference**: :cite:t:`Vinken_et_al._2011`

.. option:: SeaFlux

   Air-sea exchange.

   - Species: DMS, ACET, ALD2, MENO3, ETNO3, MOH
   - References: :cite:t:`Johnson_2010`, :cite:t:`Nightingale_et_al._2000`

.. option:: SeaSalt

   Sea salt aerosol emission.

   - **Species**: SALA, SALC, SALACL, SALCCL, SALAAL, SALCAL, BrSALA,
     BrSALC, MOPO, MOPI
   - **References**: :cite:t:`Jaegle_et_al._2011`, :cite:t:`Gong_2003`

.. option:: SoilNOx

   Emissons of NOx from soils and fertilizers.

   - **Species**: NO
   - **Reference**: :cite:t:`Hudman_et_al._2012`


.. option:: Volcano

   Emissions of volcanic SO2 from AEROCOM.

   - **Species**: SO2
   - **Reference**:


.. option:: TOMAS_Jeagle

   Size-resolved sea salt emissions for `TOMAS aerosol microphysics
   <http://wiki.geos-chem.org/TOMAS_aerosol_microphysics>`_
   simulations.

   - **Species**: SS1, SS2, SS3, SS4, SS5, SS6, SS7, SS8, SS9, SS10,
     SS11, SS12, SS13, SS14, SS15, SS16, SS17, SS18, SS19, SS20, SS21,
     SS22, SS23, SS24, SS25, SS26, SS27, SS28, SS29, SS30, SS31, SS32,
     SS33, SS34, SS35, SS36, SS37, SS38, SS39, SS40
   - **Reference**: :cite:t:`Jaegle_et_al._2011`

.. option:: TOMAS_DustDead

   Size-resolved dust emissions for `TOMAS aerosol microphysics
   <http://wiki.geos-chem.org/TOMAS_aerosol_microphysics>`_
   simulations.

   - **Species**: DUST1, DUST2, DUST3, DUST4, DUST5, DUST6, DUST7,
     DUST8, DUST9, DUST10, DUST11, DUST12, DUST13, DUST14, DUST15,
     DUST16, DUST17, DUST18, DUST19, DUST20, DUST21, DUST22, DUST23,
     DUST24, DUST25, DUST26, DUST27, DUST28, DUST29, DUST30, DUST31,
     DUST32, DUST33, DUST34, DUST35, DUST36, DUST37, DUST38, DUST39,
     DUST40
   - **Reference**: :cite:t:`Zender_et_al._2003`


.. _hco-ext-gridded-data:

============
Gridded data
============

HEMCO can host all environmentally independent data sets (e.g. source
functions) used by the extensions. The environmental variables are
either provided by the atmospheric model or directly read from file
through the HEMCO configuration file. Entries in :ref:`the HEMCO
configuration file <hco-cfg>` file are given priority over fields
passed down from the atmospheric model, i.e. if the HEMCO
configuration file contains an entry for a given environmental
variable, this field will be used instead of the field provided by the
atmospheric model. The field name provided in the HEMCO configuration
file must exactly match the name of the HEMCO environmental parameter.

To use the NCEP reanalysis monthly surface wind fields
(http:, , www.esrl.noaa.gov, psd, data, gridded, data.ncep.reanalysis.derived.surface.html)
in all HEMCO extensions, add the following two lines to the
:ref:`Base Emissions <hco-cfg-base>` section of :ref:`the HEMCO
configuration file <hco-cfg>`:

.. code-block:: kconfig

   * U10M /path/to/uwnd.mon.mean.nc uwnd 1948-2014/1-12/1/0 C xy m/s * - 1 1
   * V10M /path/to/vwnd.mon.mean.nc vwnd 1948-2014/1-12/1/0 C xy m/s * - 1 1

This will use these wind fields for all emission calculations, even if
the atmospheric model uses a different set of wind fields.

It is legal to assign scale factors (and masks) to the environmental
variables read through :ref:`the HEMCO configuration file
<hco-cfg>`. This is particularly attractive for sensitivity
studies. For example, a scale factor of 1.1 can be assigned to the
NCEP surface wind fields to study the sensitivity of emissions on a
10% increase in wind speed:

In the :ref:`Base Emissions <hco-cfg-base>` section:

.. code-block:: kconfig

   * U10M /path/to/uwnd.mon.mean.nc uwnd 1948-2014/1-12/1/0 C xy m/s * 123 1 1
   * V10M /path/to/vwnd.mon.mean.nc vwnd 1948-2014/1-12/1/0 C xy m/s * 123 1 1

In the :ref:`Scale Factors <hco-cfg-scalefac>` section:

.. code-block:: kconfig

   123 SURFWIND_SCALE 1.1 - - - xy 1 1

As for any other entry in the HEMCO configuration file, spatially
uniform values can be set directly in the HEMCO configuration file. For
example, a spatially uniform, but monthly varying surface albedo can be
specified by adding the following entry to the :ref:`Base Emissions
<hco-cfg-base>` section of :ref:`the HEMCO configuration file <hco-cfg>`:

.. code-block:: kconfig

   * ALBD 0.7/0.65/0.6/0.5/0.5/0.4/0.45/0.5/0.55/0.6/0.6/0.7 - 2000/1-12/1/0 C xy 1 * - 1 1

.. _hco-ext-env-fields:

==================================
Environmental fields used by HEMCO
==================================

The following fields can be passed from the atmospheric model to HEMCO
for use by the various extensions:

.. option:: AIR

   Air mass.

   - **Dim**: xyz
   - **Units**: kg
   - **Used by**: :option:`GC_Rn-Pb-Be`, :option:`PARANOx`

.. option:: AIRVOL

   Air volume (i.e. volume of grid box).

   - **Dim**: xyz
   - **Units**: kg
   - **Used by**: :option:`PARANOx`

.. option:: ALBD

   Surface albedo.

   - **Dim**: xy
   - **Units**: unitless
   - **Used by**: :option:`SoilNOx`, :option:`SeaFlux`

.. option:: CLDFRC

   Cloud fraction

   - **Dim**: xy
   - **Units**: unitless
   - **Used by**: :option:`MEGAN`

.. option:: CNV_MFC

   Convective mass flux.

   - **Dim**: xyz
   - **Units**: kg/m2/s
   - **Used by**: :option:`LightNOx`

.. option:: FRAC_OF_PBL

   Fraction of grid box within the planetary boundary layer (PBL).

   - **Dim**: xyz
   - **Units**: unitless
   - **Used by**: :option:`PARANOx`, :option:`SeaFlux`

.. option:: FRCLND

   Land fraction

   - **Dim**: xy
   - **Units**: unitless
   - **Used by**: :option:`GC_Rn-Pb-Be`, :option:`SeaFlux`

.. option:: GWETROOT

   Root soil moisture.


   - **Dim**: xy
   - **Units**: unitless
   - **Used by**: :option:`MEGAN`

.. option:: GWETTOP

   Top soil moisture.

   - **Dim**: xy
   - **Units**: unitless
   - **Used by**: :option:`MEGAN`

.. option:: HNO3

   HNO3 mass.

   - **Dim**: xyz
   - **Units**: kg
   - **Used by**: :option:`PARANOx`

.. option:: JO1D

   Photolysis J-value for O1D.

   - **Dim**: xy
   - **Units**: 1/s
   - **Used by**: :option:`PARANOx`

.. option:: JNO2

   Photolysis J-value for NO2.

   - **Dim**: xy
   - **Units**: 1/s
   - **Used by**: :option:`PARANOx`

.. option:: LAI

   Leaf area index.

   - **Dim**: xy
   - **Units**: cm2 leaf/cm2 grid box
   - **Used by**: :option:`MEGAN`

.. option:: NO

   NO mass.

   - **Dim**: xyz
   - **Units**: kg
   - **Used by**: :option:`PARANOx`

.. option:: NO2

   NO2 mass.

   - **Dim**: xyz
   - **Units**: kg
   - **Used by**: :option:`PARANOx`

.. option:: O3

   O3 mass.

   - **Dim**: xyz
   - **Units**: kg
   - **Used by**: :option:`PARANOx`

.. option:: PARDF

   Diffuse photosynthetic active radiation

   - **Dim**: xy
   - **Units**: W/m2
   - **Used by**: :option:`MEGAN`

.. option:: PARDR

   Direct photosynthetic active radiation

   - **Dim**: xy
   - **Units**: W/m2
   - **Used by**: :option:`MEGAN`

.. option:: RADSWG

   Short-wave incident surface radiation

   - **Dim**: xy
   - **Units**: W/m2
   - **Used by**: :option:`SoilNOx`

.. option:: SNOWHGT

   Snow height (mm of H2O equivalent).

   - **Dim**: xy
   - **Units**: kg H2O/m2
   - **Used by**: :option:`DustDead`, :option:`TOMAS_DustDead`

.. option:: SPHU

   Specific humidity

   - **Dim**: xyz
   - **Units**: kg H2O/kg air
   - **Used by**: :option:`DustDead`, :option:`PARANOx`,
     :option:`TOMAS_DustDead`

.. option:: SZAFACT

   Cosine of the solar zenith angle.

   - **Dim**: xy
   - **Units**: unitless
   - **Used by**: :option:`MEGAN`

.. option:: TK

   Temperature.

   - **Dim**: xyz
   - **Units**: K
   - **Used by**: :option:`DustDead`, :option:`LightNOx`,
     :option:`TOMAS_DustDead`

.. option:: TROPP

   Tropopause pressure.

   - **Dim**: xy
   - **Units**: Pa
   - **Used by**: :option:`GC_Rn-Pb-Be`, :option:`LightNOx`

.. option:: TSKIN

   Surface skin temperature

   - **Dim**: xy
   - **Units**: K
   - **Used by**: :option:`SeaFlux`, :option:`SeaSalt`

.. option:: U10M

   E/W wind speed @ 10 meters above surface.

   - **Dim**: xy
   - **Units**: m/s
   - **Used by**:  :option:`DustAlk`,  :option:`DustDead`,
     :option:`DustGinoux`, :option:`PARANOx`, :option:`SeaFlux`,
     :option:`SeaSalt`, :option:`SoilNOx`, :option:`TOMAS_DustDead`,
     :option:`TOMAS_Jeagle`

.. option:: USTAR

   Friction velocity.

   - **Dim**: xy
   - **Units**: m/s
   - **Used by**: :option:`DustDead`, :option:`TOMAS_DustDead`

.. option:: V10M

   N/S wind speed @ 10 meters above surface.

   - **Dim**: xy
   - **Units**: m/s
   - **Used by**:  :option:`DustAlk`,  :option:`DustDead`,
     :option:`DustGinoux`, :option:`PARANOx`, :option:`SeaFlux`,
     :option:`SeaSalt`, :option:`SoilNOx`, :option:`TOMAS_DustDead`,
     :option:`TOMAS_Jeagle`

.. option:: WLI

   Water-land-ice flags (:literal:`0` = water, :literal:`1` = land,
   :literal:`2` =  ice).

   - **Dim**: xy
   - **Units**: unitless
   - **Used by**: Almost every extension

.. option:: Z0

   Roughness height.

   - **Dim**: xy
   - **Units**: m
   - **Used by**: :option:`DustDead`, :option:`TOMAS_DustDead`

.. _hco-ext-rst-vars:

=================
Restart variables
=================

Some extensions rely on restart variables, i.e. variables that are
highly dependent on historical information such as previous-day leaf
area index or soil NOx pulsing factor. During a simulation run, the
extensions continuously archive all necessary information and update
estart variables accordingly. The updated variables become
automatically written into the HEMCO restart file
(:file:`HEMCO_restart.YYYYMMDDhhmmss.nc`) at the end of a
simulation. The fields from this file can then be read through the
HEMCO configuration file to resume the simulation at this date ("warm"
restart). For example, the soil NOx restart variables can be made
available to the soil NOx extension by adding the following lines to
the :ref:`Base Emissions section <hco-cfg-base>` of :ref:`the HEMCO
configuration file <hco-cfg>`.

.. code-block:: kconfig

   104 PFACTOR         ./HEMCO_restart.$YYYY$MM$DD$HH00.nc  PFACTOR       $YYYY/$MM/$DD/$HH E xy  unitless NO - 1 1
   104 DRYPERIOD       ./HEMCO_restart.$YYYY$MM$DD$HH00.nc  DRYPERIOD     $YYYY/$MM/$DD/$HH E xy  unitless NO - 1 1
   104 GWET_PREV       ./HEMCO_restart.$YYYY$MM$DD$HH00.nc  GWET_PREV     $YYYY/$MM/$DD/$HH E xy  unitless NO - 1 1
   104 DEP_RESERVOIR   ./HEMCO_restart.$YYYY$MM$DD$HH00.nc  DEP_RESERVOIR $YYYY/$MM/$DD/$HH E xy  unitless NO - 1 1

Many restart variables are very time and date-dependent. It is therefore
recommended to set the time slice selection flag to E to ensure that
only data is read that exactly matches the simulation start date (also
see :ref:`hco-cfg-base`.  HEMCO will perform a "cold start" if no
restart field can be found for a given simulation start date,
e.g. default values will be used for those restart variables.

.. _built_in_tools_for_scalingmasking:

==================================
Built-in tools for scaling/masking
==================================

HEMCO has built-in tools to facilitate the application of both uniform
and spatiotemporal :ref:`scale factors <hco-cfg-scalefac>` to
emissions calculated by the extensions. At this point, not all
extensions take advantage of these tools yet. A list of extensions
that support the built-in scaling tools are given below.

For extensions that support the built-in scaling tools, you can specify
the uniform and/or spatiotemporal scale factors to be applied to the
extension species of interest in section :ref:`hco-cfg-ext-switches`
:ref:`the HEMCO configuration file <hco-cfg>`.

For example, to uniformly scale GFED CO by a factor of 1.05 and GFED NO
emissions by a factor of 1.2, add the following two lines to the HEMCO
configuration file (highlighted in GREEN):

.. code-block:: kconfig

   111    GFED              : on    CO/NO/ACET/ALK4
      --> GFED3             :       false
      --> GFED4             :       true
      --> GFED_daily        :       false
      --> GFED_3hourly      :       false
      --> Scaling_CO        :       1.05
      --> Scaling_NO        :       1.20

Similarly, a spatiotemporal field to be applied to the species of
interest can be defined via setting :literal:`ScaleField`, e.g.

.. code-block:: kconfig

   111     GFED              : on    CO/NO/ACET/ALK4
       --> GFED3             :       false
       --> GFED4             :       true
       --> GFED_daily        :       false
       --> GFED_3hourly      :       false
       --> Scaling_CO        :       1.05
       --> Scaling_NO        :       1.20
       --> ScaleField_NO     :       GFED_SCALEFIELD_NO

The corresponding scale field needs be defined in section
:ref:`hco-cfg-base` . A simple example would be a monthly
varying scale factor for GFED NO emissions:

.. code-block:: kconfig

   111 GFED_SCALEFIELD_NO   0.9/1.1/1.3/1.4/1.6/1.7/1.7/1.8/1.5/1.3/0.9/0.8 - 2000/1-12/1/0 C xy unitless * - 1 1

It is legal to apply scale factors and/or masks to the extension scale
fields (in the same way as the 'regular' base emission fields). A more
sophisticated example on how to scale soil NOx emissions is given in
HEMCO examples.

.. _hco-ext-scale-mask:

==============================================
Extensions supporting built-in scaling/masking
==============================================

The following extensions currently support the built-in scaling/masking
tools: :option:`SoilNOx`, :option:`GFED`, :option:`FINN`.

===========================
Adding new HEMCO extensions
===========================

All HEMCO extensions are called through the extension interface
routines in :file:`HEMCO/Extensions/hcox_driver_mod.F90`:
:code:`HCOX_INIT`, :code:`HCOX_RUN`, :code:`HCOX_FINAL`. For
every new extension, a corresponding subroutine call needs to be added
to those three routines.  You will quickly see that these calls only
take a few arguments, most importantly the HEMCO state object
:code:`HcoState` and the extensions state object :code:`ExtState`.

:code:`ExtState` is defined in
:file:`HEMCO/src/Extensions/hcox_state_mod.F90`. It contains logical
switches for each extension as well as pointers to any external data
(such as met fields). For a new extension, you'll have to add a new
logical switch to the Ext_State object. If you need external data that
is not yet included in ExtState, you will also have to add those
(including the pointer associations in subroutine
:code:`SET_EXTOPT_FIELDS` in
:file:`GeosCore/hco_interface_gc_mod.F90`.

The initialization call (:code:`HCOX_XXX_INIT`) should be used to
initialize all extension variables and to read all settings from the
HEMCO configuration file. There are a number of helper routines in
:file:`HEMCO/src/Extensions/hco_extlist_mod.F90` to do this:

- :code:`GetExtNr( ExtName )` returns the extension number for the
  given extension name. Will return –1 if extension is turned off!

- :code:`GetExtOpt( ExtNr, Attribute, Value, RC )` can be used to read
  any additional extension options (logical switches, path and names
  of csv-tables, etc.). Note that value can be of various types
  (:code:`logical`, :code:`character`, :code:`double`,...).

- :code:`GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )`
  matches the extension species names (as defined in the configuration
  file) to the species defined in HEMCO state (i.e. to all available
  HEMCO species). A value of –1 is returned if the given species is
  not defined in HEMCO.

All :code:`ExtState` variables used by this extension should be
updated. This includes the logical switch and all external data needed
by the extension. For example, if the extension needs temperature
data, this pointer should be activated by setting
:code:`ExtState%TK%DoUse = .TRUE.`

The run call (:code:`HCOX_XXX_RUN`) calculates the 2D fluxes and
passes them to HcoState via subroutine :code:`HCO_EmisAdd( HcoState,
Flux, HcoID, RC)`. External data is assessed through :code:`ExtState`
(e.g. :code:`ExtState%TX%Arr%Val(I,J,L)`), and any data automatically
read from netCDF files (through the HEMCO interface) can be obtained
through :code:`EmisList_GetDataArr( am_I_Root, FieldName, Pointer,
RC )` The body of the run routine is typically just the code of the
original module.

It's probably easiest to start from an existing extension (or the
:file:`Custom` extension template) and to add any modifications as is
needed.
