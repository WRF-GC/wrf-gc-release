.. _edit-hco-cfg:

##############
Basic examples
##############

.. note::

   The following sections contain simple HEMCO configuration file
   examples for demonstration purposes.  If you are using HEMCO with
   an external model, then your HEMCO configuration file may be more
   complex than the examples shown below.

All emission calculation settings are specified in :ref:`the HEMCO
configuration file <hco-cfg>`, which is named :file:`HEMCO_Config.rc`.

Modification of the HEMCO source code (and recompilation) is only
required if new extensions are added, or to use HEMCO in a new model
environent (see sections :ref:`hco-hood` and :ref:`hco-hood-int`).

In the sections that follow, we provide some basic examples that
demonstrate how to modify the configuration file to customize your
HEMCO simulation.

.. _edit-hco-cfg-ex1:

=============================================
Example 1: Add global anthropogenic emissions
=============================================

Suppose monthly global anthropogenic CO emissions from the **MACCity**
inventory :cite:`Lamarque_et_al._2010` are stored in file
:file:`MACCity.nc` as variable :literal:`CO`. The following HEMCO
configuration file then simulates CO emissions with gridded
hourly scale factors applied to it (the latter taken from variable
:literal:`factor` of file :file:`hourly.nc`).

The horizontal grid and simulation datetimes employed by HEMCO depends
on the HEMCO-to-model interface. If HEMCO is coupled to an external
model (such as `GEOS-Chem <https://geos-chem.readthedocs.io>`_) these
values are taken from the chemistry model. If run standalone, the grid
specification and desired datetimes need be specified as described in
Interfaces.

.. code-block:: kconfig

   ###############################################################################
   ### BEGIN SECTION SETTINGS
   ###############################################################################
   ROOT:                        /dir/to/data
   Logfile:                     HEMCO.log
   DiagnFile:                   HEMCO_Diagn.rc
   DiagnPrefix:                 HEMCO_diagnostics
   Wildcard:                    *
   Separator:                   /
   Unit tolerance:              1
   Negative values:             0
   Only unitless scale factors: false
   Verbose:                     0
   Warnings:                    1

   ### END SECTION SETTINGS ###

   ###############################################################################
   ### BEGIN SECTION EXTENSION SWITCHES
   ###############################################################################
   # ExtNr ExtName           on/off  Species
   0       Base              : on    *
       --> MACCITY           :       true

   ### END SECTION EXTENSION SWITCHES ###

   ###############################################################################
   ### BEGIN SECTION BASE EMISSIONS
   ###############################################################################
   # ExtNr Name sourceFile sourceVar sourceTime C/R/E SrcDim SrcUnit Species ScalIDs Cat Hier

   (((MACCITY
   0 MACCITY_CO $ROOT/MACCity.nc  CO 1980-2014/1-12/1/0 C xy kg/m2/s CO 500 1 1
   )))MACCITY

   ### END SECTION BASE EMISSIONS ###

   ###############################################################################
   ### BEGIN SECTION SCALE FACTORS
   ###############################################################################
   # ScalID Name srcFile srcVar srcTime  CRE Dim Unit Oper

   500 HOURLY_SCALFACT $ROOT/hourly.nc factor 2000/1/1/0-23 C xy 1 1

   ### END SECTION SCALE FACTORS ###

   ###############################################################################
   ### BEGIN SECTION MASKS
   ###############################################################################

   ### END SECTION MASKS ###

The various attributes are explained in more detail in the
:ref:`hco-cfg-base` and :ref:`hco-cfg-scalefac` sections.

.. note::

   We have used an index of 500 for :literal:`HOURLY_SCALFACT` in
   order to reduce confusion with the :literal:`Cat` and
   :literal:`Hier` values.

As described in :ref:`hco-cfg-data-coll`, all of the files
contained between the brackets :literal:`(((MACCITY` and
:literal:`)))MACCITY` will be read if you set the switch

.. code-block:: text

   --> MACCITY           :       true

These files will be ignored if you set

.. code-block::

   --> MACCITY           :       false

This is a quick way to shut off individual emissions inventories without
having to manually comment out many lines of code. You can add a set of
brackets, with a corresponding true/false switch, for each emissions
inventory that you add to the configuration file.

.. _edit-hco-cfg-ex2:

=====================================
Example 2: Overlay regional emissions
=====================================

To add regional monthly anthropogenic CO emissions from the EMEP
European inventory :cite:`Vestreng_et_al._2009` (in file
:file:`EMEP.nc`)  to the simulation, modify the configuration file as
follows:

.. code-block:: kconfig

    ###############################################################################
    #### BEGIN SECTION EXTENSION SWITCHES
    ###############################################################################
    # ExtNr ExtName           on/off  Species
    0       Base              : on    *
        --> MACCITY           :       true
        --> EMEP              :       true

    ### END SECTION EXTENSION SWITCHES ###

    ###############################################################################
    ### BEGIN SECTION BASE EMISSIONS
    ###############################################################################
    #ExtNr Name srcFile srcVar srcTime CRE Dim Unit Species ScalIDs Cat Hier

    (((MACCITY
    0 MACCITY_CO $ROOT/MACCity.nc CO 1980-2014/1-12/1/0 C xy kg/m2/s CO  500      1 1
    )))MACCITY

    (((EMEP
    0 EMEP_CO    $ROOT/EMEP.nc    CO 2000-2014/1-12/1/0 C xy kg/m2/s CO  500/1001 1 2
    )))EMEP

    ### END SECTION BASE EMISSIONS###

    ###############################################################################
    ### BEGIN SECTION SCALE FACTORS
    ###############################################################################
    #ScalID Name srcFile srcVar srcTime CRE Dim Unit Oper

    500 HOURLY_SCALFACT $ROOT/hourly.nc factor 2000/1/1/0-23 C xy 1 1

    ### END SECTION SCALE FACTORS ###

    ###############################################################################
    ### BEGIN SECTION MASKS
    ###############################################################################
    #ScalID Name srcFile srcVar srcTime CRE Dim Unit Oper Box

    1001 MASK_EUROPE $ROOT/mask_europe.nc MASK 2000/1/1/0 C xy 1 1 -30/30/45/70

    ### END SECTION MASKS ###

For now, we have omitted the **Settings section**  because nothing has
changed since :ref:`the previous example <edit-hco-cfg-ex1>`.

Note the increased hierarchy (:literal:`2`) of the regional EMEP
inventory compared to the global MACCity emissions (:literal:`1`) in
column :option:`Hier`. This will cause the EMEP emissions to replace
the MACCity emissions in the region where EMEP is defined, which is
specified by the MASK_EUROPE variable.

.. _edit-hco-cfg-ex3:

=============================================
Example 3: Adding the AEIC aircraft emissions
=============================================

To add aircraft emissions from the AEIC inventory
:cite:`Stettler_et_al._2011`, available in file :file:`AEIC.nc`,
modify the :ref:`configuration file <hco-cfg>` accordingly:

.. code-block :: kconfig

   ###############################################################################
   #### BEGIN SECTION EXTENSION SWITCHES
   ###############################################################################
   # ExtNr ExtName           on/off  Species
   0       Base              : on    *
       --> MACCITY           :       true
       --> EMEP              :       true
       --> AEIC              :       true
   ### END SECTION EXTENSION SWITCHES ###

   ###############################################################################
   #### BEGIN SECTION BASE EMISSIONS
   ###############################################################################
   #ExtNr Name srcFile srcVar srcTime CRE Dim Unit Species ScalIDs Cat Hier

   (((MACCITY
   0 MACCITY_CO $ROOT/MACCity.nc CO 1980-2014/1-12/1/0 C xy  kg/m2/s CO 500        1 1
   )))MACCITY

   (((EMEP
   0 EMEP_CO    $ROOT/EMEP.nc    CO 2000-2014/1-12/1/0 C xy  kg/m2/s CO 500 1/1001 1 2
   )))EMEP

   (((AEIC
   0 AEIC_CO    $ROOT/AEIC.nc    CO 2005/1-12/1/0      C xyz kg/m2/s CO -          2 1
   )))AEIC

   ### END SECTION BASE EMISSIONS ###

Note the change in the emission category (column :option:`Cat`) from
:literal:`1` to :literal:`2`.  In this example, category 1 represents
anthropogenic emissions and category 2 represents aircraft emissions.

.. _edit-hco-cfg-ex4:

========================================
Example 4: Add biomass burning emissions
========================================

GFED4 biomass burning emissions (Giglio et al, 2013), which are
implemented as a HEMCO Extension, can be added to the simulation by:

#. Adding the corresponding extension to section **Extension
   Switches**
#. Adding all the input data needed by GFED4 to section **Base
   Emissions**.

The extension number defined in the **Extension Switches** section
must match the corresponding :option:`ExtNr` entry in the Base
Emissions section (in this example, :literal:`111`).

.. code-block:: kconfig

   ###############################################################################
   #### BEGIN SECTION EXTENSION SWITCHES
   ###############################################################################
   # ExtNr ExtName           on/off  Species
   0       Base              : on    *
       --> MACCITY           :       true
       --> EMEP              :       true
       --> AEIC              :       true
   #------------------------------------------------------------------------------
   111     GFED              : on    CO
       --> GFED3             :       false
       --> GFED4             :       true
       --> GFED_daily        :       false
       --> GFED_3hourly      :       false
       --> Scaling_CO        :       1.05

   ### END SECTION EXTENSION SWITCHES ###

   ###############################################################################
   #### BEGIN SECTION BASE EMISSIONS
   ###############################################################################
   #ExtNr Name srcFile srcVar srcTime CRE Dim Unit Species ScalIDs Cat Hier

   (((MACCITY
   0 MACCITY_CO $ROOT/MACCity.nc  CO 1980-2014/1-12/1/0 C xy  kg/m2/s CO 500      1 1
   )))MACCITY

   (((EMEP
   0 EMEP_CO    $ROOT/EMEP.nc     CO 2000-2014/1-12/1/0 C xy  kg/m2/s CO 500/1001 1 2
   )))EMEP

   (((AEIC
   0 AEIC_CO    $ROOT/AEIC.nc     CO 2005/1-12/1/0      C xyz kg/m2/s CO -        2 1
   )))AEIC

   ###############################################################################
   ###  BEGIN SECTION EXTENSION DATA (subsection of BASE EMISSIONS SECTION
   ###
   ### These fields are needed by the extensions listed above. The assigned ExtNr
   ### must match the ExtNr entry in section 'Extension switches'. These fields
   ### are only read if the extension is enabled.  The fields are imported by the
   ### extensions by field name.  The name given here must match the name used
   ### in the extension's source code.
   ###############################################################################

   # --- GFED biomass burning emissions (Extension 111) ---
   111 GFED_HUMTROP    $ROOT/GFED3/v2014-10/GFED3_humtropmap.nc              humtrop           2000/1/1/0             C xy 1         * - 1 1

   (((GFED3
   111 GFED_WDL        $ROOT/GFED3/v2014-10/GFED3_gen.1x1.$YYYY.nc           GFED3_BB__WDL_DM  1997-2011/1-12/01/0    C xy kgDM/m2/s * - 1 1
   111 GFED_AGW        $ROOT/GFED3/v2014-10/GFED3_gen.1x1.$YYYY.nc           GFED3_BB__AGW_DM  1997-2011/1-12/01/0    C xy kgDM/m2/s * - 1 1
   111 GFED_DEF        $ROOT/GFED3/v2014-10/GFED3_gen.1x1.$YYYY.nc           GFED3_BB__DEF_DM  1997-2011/1-12/01/0    C xy kgDM/m2/s * - 1 1
   111 GFED_FOR        $ROOT/GFED3/v2014-10/GFED3_gen.1x1.$YYYY.nc           GFED3_BB__FOR_DM  1997-2011/1-12/01/0    C xy kgDM/m2/s * - 1 1
   111 GFED_PET        $ROOT/GFED3/v2014-10/GFED3_gen.1x1.$YYYY.nc           GFED3_BB__PET_DM  1997-2011/1-12/01/0    C xy kgDM/m2/s * - 1 1
   111 GFED_SAV        $ROOT/GFED3/v2014-10/GFED3_gen.1x1.$YYYY.nc           GFED3_BB__SAV_DM  1997-2011/1-12/01/0    C xy kgDM/m2/s * - 1 1
   )))GFED3

   (((GFED4
   111 GFED_WDL        $ROOT/GFED4/v2015-03/GFED4_gen.025x025.$YYYY.nc       WDL_DM            2000-2013/1-12/01/0    C xy kg/m2/s   * - 1 1
   111 GFED_AGW        $ROOT/GFED4/v2015-03/GFED4_gen.025x025.$YYYY.nc       AGW_DM            2000-2013/1-12/01/0    C xy kg/m2/s   * - 1 1
   111 GFED_DEF        $ROOT/GFED4/v2015-03/GFED4_gen.025x025.$YYYY.nc       DEF_DM            2000-2013/1-12/01/0    C xy kg/m2/s   * - 1 1
   111 GFED_FOR        $ROOT/GFED4/v2015-03/GFED4_gen.025x025.$YYYY.nc       FOR_DM            2000-2013/1-12/01/0    C xy kg/m2/s   * - 1 1
   111 GFED_PET        $ROOT/GFED4/v2015-03/GFED4_gen.025x025.$YYYY.nc       PET_DM            2000-2013/1-12/01/0    C xy kg/m2/s   * - 1 1
   111 GFED_SAV        $ROOT/GFED4/v2015-03/GFED4_gen.025x025.$YYYY.nc       SAV_DM            2000-2013/1-12/01/0    C xy kg/m2/s   * - 1 1
   )))GFED4

   (((GFED_daily
   111 GFED_FRAC_DAY   $ROOT/GFED3/v2014-10/GFED3_dailyfrac_gen.1x1.$YYYY.nc GFED3_BB__DAYFRAC 2002-2011/1-12/1-31/0  C xy 1         * - 1 1
   )))GFED_daily

   (((GFED_3hourly
   111 GFED_FRAC_3HOUR $ROOT/GFED3/v2014-10/GFED3_3hrfrac_gen.1x1.$YYYY.nc   GFED3_BB__HRFRAC  2002-2011/1-12/01/0-23 C xy 1         * - 1 1
   )))GFED_3hourly

   ### END SECTION BASE EMISSIONS ###

As in the previous examples, the tags beginning with :literal:`(((` and
:literal:`)))` denote options that can be toggled on or off in the
Extension Switches section. For example, if you wanted to use GFED3
biomass emissions instead of GFED4, you would set the switch for GFED3
to true and the switch for GFED4 to false.

Scale factors and other extension options (e.g. :literal:`Scaling_CO`)
can be specified in the Extension Switches section.

.. _edit-hco-cfg-ex5:

===============================================
Example 5: Tell HEMCO to use additional species
===============================================

The HEMCO configuration file can hold emission specifications of as
many species as desired. For example, to add anthropogenic NO
emissions from the MACCity inventory, modify the HEMCO configuration
file as shown:

.. code-block:: kconfig

   ###############################################################################
   #### BEGIN SECTION BASE EMISSIONS
   ###############################################################################
   #ExtNr Name srcFile srcVar srcTime CRE Dim Unit Species ScalIDs Cat Hier

   (((MACCITY
   0 MACCITY_CO $ROOT/MACCity.nc CO 1980-2014/1-12/1/0 C xy kg/m2/s CO 500 1 1
   0 MACCITY_NO $ROOT/MACCity.nc NO 1980-2014/1-12/1/0 C xy kg/m2/s NO 500 1 1
   )))MACCITY

To include NO in GFED, we can just add NO to the list of species that
GFED will process in the Extension Switches section.

.. code-block:: kconfig

   ###############################################################################
   #### BEGIN SECTION EXTENSION SWITCHES
   ###############################################################################
   # ExtNr ExtName           on/off  Species
   0       Base              : on    *
       --> MACCITY           :       true
       --> EMEP              :       true
       --> AEIC              :       true
   #------------------------------------------------------------------------------
   111     GFED              : on    CO/NO
       --> GFED3             :       false
       --> GFED4             :       true
       --> GFED_daily        :       false
       --> GFED_3hourly      :       false
       --> Scaling_CO        :       1.05

Finally, let's add sulfate emissions to the simulation. Emissions of
SO4 are approximated from the MACCity SO2 data, assuming that SO4
constitutes 3.1% of the SO2 emissions. The final configuration file
now looks like this:

.. code-block:: kconfig

   ###############################################################################
   #### BEGIN SECTION SETTINGS
   ###############################################################################
   ROOT:                        /dir/to/data
   Logfile:                     HEMCO.log
   DiagnFile:                   HEMCO_Diagn.rc
   DiagnPrefix:                 HEMCO_diagnostics
   Wildcard:                    *
   Separator:                   /
   Unit tolerance:              1
   Negative values:             0
   Only unitless scale factors: false
   Verbose:                     0
   Warnings:                    1

   ### END SECTION SETTINGS ###

   ###############################################################################
   ### BEGIN SECTION EXTENSION SWITCHES
   ###############################################################################
   # ExtNr ExtName           on/off  Species
   0       Base              : on    *
       --> MACCITY           :       true
       --> EMEP              :       true
       --> AEIC              :       true
   #------------------------------------------------------------------------------
   111     GFED              : on    CO/NO/SO2
       --> GFED3             :       false
       --> GFED4             :       true
       --> GFED_daily        :       false
       --> GFED_3hourly      :       false
       --> Scaling_CO        :       1.05

   ### END SECTION EXTENSION SWITCHES ###

   ###############################################################################
   #### BEGIN SECTION BASE EMISSIONS
   ###############################################################################
   #ExtNr Name srcFile srcVar srcTime CRE Dim Unit Species ScalIDs Cat Hier
   (((MACCITY
   0 MACCITY_CO  $ROOT/MACCity.nc CO  1980-2014/1-12/1/0 C xy  kg/m2/s CO  500     1 1
   0 MACCITY_NO  $ROOT/MACCity.nc NO  1980-2014/1-12/1/0 C xy  kg/m2/s NO  500     1 1
   0 MACCITY_SO2 $ROOT/MACCity.nc SO2 1980-2014/1-12/1/0 C xy  kg/m2/s SO2 -       1 1
   0 MACCITY_SO4 -                -   -                  - -   -       SO4 600     1 1
   )))MACCITY

   (((EMEP
   0 EMEP_CO     $ROOT/EMEP.nc     CO 2000-2014/1-12/1/0 C xy  kg/m2/s CO 500/1001 1 2
   )))EMEP

   (((AEIC
   0 AEIC_CO     $ROOT/AEIC.nc     CO 2005/1-12/1/0      C xyz kg/m2/s CO -        2 1
   )))AEIC

   ###############################################################################
   ###  BEGIN SECTION EXTENSION DATA (subsection of BASE EMISSIONS SECTION
   ###
   ### These fields are needed by the extensions listed above. The assigned ExtNr
   ### must match the ExtNr entry in section 'Extension switches'. These fields
   ### are only read if the extension is enabled.  The fields are imported by the
   ### extensions by field name.  The name given here must match the name used
   ### in the extension's source code.
   ##############################################################################

   # --- GFED biomass burning emissions (Extension 111) ---
   111 GFED_HUMTROP    $ROOT/GFED3/v2014-10/GFED3_humtropmap.nc              humtrop           2000/1/1/0             C xy 1         * - 1 1

   (((GFED3
   111 GFED_WDL        $ROOT/GFED3/v2014-10/GFED3_gen.1x1.$YYYY.nc           GFED3_BB__WDL_DM  1997-2011/1-12/01/0    C xy kgDM/m2/s * - 1 1
   111 GFED_AGW        $ROOT/GFED3/v2014-10/GFED3_gen.1x1.$YYYY.nc           GFED3_BB__AGW_DM  1997-2011/1-12/01/0    C xy kgDM/m2/s * - 1 1
   111 GFED_DEF        $ROOT/GFED3/v2014-10/GFED3_gen.1x1.$YYYY.nc           GFED3_BB__DEF_DM  1997-2011/1-12/01/0    C xy kgDM/m2/s * - 1 1
   111 GFED_FOR        $ROOT/GFED3/v2014-10/GFED3_gen.1x1.$YYYY.nc           GFED3_BB__FOR_DM  1997-2011/1-12/01/0    C xy kgDM/m2/s * - 1 1
   111 GFED_PET        $ROOT/GFED3/v2014-10/GFED3_gen.1x1.$YYYY.nc           GFED3_BB__PET_DM  1997-2011/1-12/01/0    C xy kgDM/m2/s * - 1 1
   111 GFED_SAV        $ROOT/GFED3/v2014-10/GFED3_gen.1x1.$YYYY.nc           GFED3_BB__SAV_DM  1997-2011/1-12/01/0    C xy kgDM/m2/s * - 1 1
   )))GFED3

   (((GFED4
   111 GFED_WDL        $ROOT/GFED4/v2015-03/GFED4_gen.025x025.$YYYY.nc       WDL_DM            2000-2013/1-12/01/0    C xy kg/m2/s   * - 1 1
   111 GFED_AGW        $ROOT/GFED4/v2015-03/GFED4_gen.025x025.$YYYY.nc       AGW_DM            2000-2013/1-12/01/0    C xy kg/m2/s   * - 1 1
   111 GFED_DEF        $ROOT/GFED4/v2015-03/GFED4_gen.025x025.$YYYY.nc       DEF_DM            2000-2013/1-12/01/0    C xy kg/m2/s   * - 1 1
   111 GFED_FOR        $ROOT/GFED4/v2015-03/GFED4_gen.025x025.$YYYY.nc       FOR_DM            2000-2013/1-12/01/0    C xy kg/m2/s   * - 1 1
   111 GFED_PET        $ROOT/GFED4/v2015-03/GFED4_gen.025x025.$YYYY.nc       PET_DM            2000-2013/1-12/01/0    C xy kg/m2/s   * - 1 1
   111 GFED_SAV        $ROOT/GFED4/v2015-03/GFED4_gen.025x025.$YYYY.nc       SAV_DM            2000-2013/1-12/01/0    C xy kg/m2/s   * - 1 1
   )))GFED4

   (((GFED_daily
   111 GFED_FRAC_DAY   $ROOT/GFED3/v2014-10/GFED3_dailyfrac_gen.1x1.$YYYY.nc GFED3_BB__DAYFRAC 2002-2011/1-12/1-31/0  C xy 1         * - 1 1
   )))GFED_daily

   (((GFED_3hourly
   111 GFED_FRAC_3HOUR $ROOT/GFED3/v2014-10/GFED3_3hrfrac_gen.1x1.$YYYY.nc   GFED3_BB__HRFRAC  2002-2011/1-12/01/0-23 C xy 1         * - 1 1
   )))GFED_3hourly

   ### END SECTION BASE EMISSIONS ###

   ###############################################################################
   #### BEGIN SECTION SCALE FACTORS
   ###############################################################################
   # ScalID Name srcFile srcVar srcTime CRE Dim Unit Oper

   500 HOURLY_SCALFACT $ROOT/hourly.nc factor  2000/1/1/0-23 C xy 1 1
   600 SO2toSO4        0.031           -       -             - -  1 1

   ### END SECTION SCALE FACTORS ###

   ###############################################################################
   #### BEGIN SECTION MASKS
   ###############################################################################
   #ScalID Name srcFile srcVar srcTime CRE Dim Unit Oper Box

   1001 MASK_EUROPE $ROOT/mask_europe.nc MASK 2000/1/1/0 C xy 1 1 -30/30/45/70

   ### END SECTION MASKS ###

.. _edit-hco-cfg-ex6:

======================================================================================
Example 6: Add inventories that do not separate out biofuels and/or trash emissions
======================================================================================

Several emissions inventories (e.g. CEDS and EDGAR) lump biofuels
and/or and trash emissions together with anthropogenic emissions. For
inventories such as these, HEMCO allows you to specify up to 3
multiple categories for each species listing in the HEMCO
configuration file. All of the emissions will go into the first listed
category, and the other listed categories will be set to zero.

In this example, all NO emissions from the EDGAR inventory power
sector will be placed into the the anthropogenic emissions category
(:literal:`Cat=1`), while the biofuel emissions category (Cat=2) will
be set to zero.

.. code-block:: kconfig

   0 EDGAR_NO_POW EDGAR_v43.NOx.POW.0.1x0.1.nc emi_nox 1970-2010/1/1/0 C xy kg/m2/s NO 1201/25/115  1/2  2

In this example, all NO emissions from CEDS inventory agriculture
sector will be placed into the the anthropogenic emissions category
(:literal:`Cat=1`), while the biofuel emissions category
(:literal:`Cat=2`) and trash emissions category (:literal:`Cat=12`)
will be set to zero.

.. code-block:: kconfig

   0 CEDS_NO_AGR NO-em-anthro_CMIP_CEDS_$YYYY.nc  NO_agr 1750-2014/1-12/1/0 C xy kg/m2/s NO  25 1/2/12 5
