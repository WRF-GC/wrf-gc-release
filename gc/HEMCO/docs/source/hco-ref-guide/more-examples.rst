.. |br| raw:: html

   <br />

.. _cfg-ex:

###########################
More configuration examples
###########################

.. _cfg-ex-scl:

=====================
Scale factor examples
=====================

.. _cfg-ex-scl-shapefile:

Scale (or zero) emissions with a shapefile country mask
-------------------------------------------------------

HEMCO has the ability to define country-specific scale factors. To
utilize this feature, you must first specify a mask file in the
**NON-EMISSIONS DATA** section of :ref:`the HEMCO configuration file
<hco-cfg>`, such as:

.. code-block:: kconfig

   #==============================================================================
   # --- Country mask file ---
   #==============================================================================
   * COUNTRY_MASK /path/to/file/countrymask_0.1x0.1.nc CountryID 2000/1/1/0 C xy count * - 1 1

The mask file specified above was created from a shapefile obtained
from the `GADM database <http://www.gadm.org>`_. The country mask
netCDF file (`countrymask_0.1x0.1.nc
<http://geoschemdata.wustl.edu/ExtData/HEMCO/MASKS/v2014-07/countrymask_0.1x0.1.nc>`_
) identifies countries by their ISO 3166-1 numeric code. Countries and
their ISO3166-1-numeric codes are listed in the `country_codes.csv
<http://geoschemdata.wustl.edu/ExtData/HEMCO/MASKS/v2014-07/country_codes.csv>`_
file.

The country-specific scale factors can be specified in a separate
ASCII file ending with with the suffix :literal:`.txt.` The container
name of the mask file (e.g. :literal:`COUNTRY_MASK`) must be given in
the first line of the file. The following lines define the
country-specific scale factors. ID 0 is reserved for the default
values that are applied to all countries with no specific values
listed. An example :file:`scalefactor.txt` file is provided below:

.. code-block:: kconfig

   # Country mask field name
   COUNTRY_MASK

   # Country data
   # Name   | ID  | Scale factor
   DEFAULT    0     1.0
   CHINA      156   0.95
   INDIA      356   1.10
   KOREA      410   0.0

The scale factor(s) listed are interpreted by HEMCO the same way as
other scale factors. Multiple values separated by :literal:`/` are
interpreted as temporally changing values:

  - 7 values = Sun, Mon, ..., Sat;
  - 12 values = Jan, Feb, ..., Dec;
  - 24 values = 12am, 1am, ..., 11pm (local time!).

The country-specific scale factors would then be defined in the
:ref:`Scale Factors <hco-cfg-scalefac>` section of :ref:`the HEMCO
configuration file <hco-cfg>` as:

.. code-block:: kconfig

   501 SCALE_COUNTRY /path/to/file/scalefactor.txt  - - - xy count 1

The scale factors can the be applied to the emission field(s) that you
wish to scale. For example:

.. code-block:: kconfig

   0 MIX_NO_IND MIX_Asia_NO.generic.025x025.nc NO_INDUSTRY 2008-2010/1-12/1/0 C xy kg/m2/s NO  1/27/25/1006/ 501 1/2 45

These steps can also be used to scale emissions for different regions
(e.g. provinces, states) by providing HEMCO with a mask file
containing the regions to be scaled.


.. _cfg-ex-scl-rec-mask:

Scale (or zero) emissions with a rectangular mask
-------------------------------------------------

.. important::

   If you are using HEMCO versions prior to 3.5.0, you may encounter a
   bug when trying to follow this example. See Github issue:
   https://github.com/geoschem/HEMCO/issues/153 for a workaround.

Another way to scale all emissions over a country (or set them to
zero) is to apply a rectangular mask.

For example, to set all emissions over Australia and surrounding
islands to zero, add this line to the :ref:`hco-cfg-masks` section of
:ref:`the HEMCO configuration file <hco-cfg>`:

.. code-block:: kconfig

    1010 AUS_MASK 105.0/-46.0/160.0/-10.0 - 2000/1/1/0 C xy 1 1 105/-46/160/–10

Here you directly provide the lower left and upper right corner of the
mask region mask instead of a netCDF file:
:literal:`lon1/lat1/lon2/lat2` You can then combine this mask with
a scale factor of zero to eliminate any emissions over that area.

In :ref:`Base emissions <hco-cfg-base>`

.. code-block:: kconfig

    0 HTAP_NO_IND /path/to/HTAP_NO_INDUSTRY.generic.01x01.nc emi_no 2008-2010/1-12/1/0 C xy kg/m2/s NO 1/27/25/501 1/2 4

In :ref:`Scale Factors <hco-cfg-scalefac>`:

.. code-block:: kconfig

   501 SCALE_AUS 0.0 - - - xy unitless 1 1010

In :ref:`hco-cfg-masks`:

.. code-block:: kconfig

   # Defines a rectangular region that should cover AUS + surrounding islands
   1010 AUS_MASK 105.0/-46.0/160.0/-10.0 – 2000/1/1/0 C xy 1 1 105.0/-46.0/160.0/-10.0

.. _cfg-ex-scl-spc:

Scale emissions by species
--------------------------

You may define uniform scale factors for single species that
apply across all emission inventories, sectors and extensions. These
scale factors can be set in the :ref:`Settings <hco-cfg-settings>`
section of :ref:`the HEMCO configuration file <hco-cfg>`, using the
:literal:`EmissScale_<species-name>`, where :literal:`<species-name>`
denotes the name of a HEMCO species such as :literal:`CO`,
:literal:`CH4`, :literal:`NO`, etc.

For instance, to scale all NO emissions by 50% add the line
:literal:`EmisScale_NO` to the :ref:`Settings <hco-cfg-settings>`
section of the :ref:`the HEMCO configuration file <hco-cfg>`:

.. code-block:: kconfig

   ###############################################################################
   ### BEGIN SECTION SETTINGS
   ###############################################################################

   ROOT:                        /path/to/HEMCO/data/directory
   Logfile:                     HEMCO.log
   ... etc ...
   EmisScale_NO                 1.5

   ### END SECTION SETTINGS ###

.. _cfg-ex-scl-zero-spc:

Zero emissions of selected species
----------------------------------

To zero emissions of a given species (e.g. NO) from any inventory
listed under :ref:`Base Emissions <hco-cfg-base>`, do the following:

#. Create your own scale factor and assign value 0.0 to it. This must
   go into the :ref:`Scale Factors <hco-cfg-scalefac>` section of
   :ref:`the HEMCO configuration file <hco-cfg>`:

   .. code-block:: kconfig

      400 ZERO 0.0 - - - xy 1 1

#. Apply this scale factor to all of the emissions entries in the
   HEMCO configuration file that you would like to zero out.  For
   example:

   .. code-block:: kconfig

      0 MIX_NO_IND  /path/to/MIX_Asia_NO.generic.025x025.nc NO_INDUSTRY  2008-2010/1-12/1/0 C xy kg/m2/s  NO  1/27/25/400/1006 1/2 45

This can be a useful way to set the emissions of some species to zero
for sensitivity study purposes.

.. note::

   All scale factors should be listed before masks.

.. _cfg-ex-ext-global:

Scale extension emissions globally by species
---------------------------------------------

You may pass a global scale factor to the :ref:`hco-ext`.  For
example, to double soil NO emissions everywhere, add the
:literal:`Scaling_NO` to the section for the :option:`SoilNOx`
extension.  This is located in the :ref:`Extension Switches
<hco-cfg-ext-switches>` section of :ref:`the HEMCO configuration file
<hco-cfg>`, as shown below:

.. code-block:: kconfig

   104     SoilNOx           : on    NO
       --> Use fertilizer NOx:       true
       --> Scaling_NO        :       2.0

.. _cfg-ex-summer-nox:

Scale summertime soil NOx emisions over the US
----------------------------------------------

It is possible to pass uniform and/or spatiotemporal scale factors to
some of the extensions, including :option:`SoilNOx`.

For instance, suppose you want to halve summertime soil NOx emissions
over the continental US. You can do this by defining a scale field
(here, :literal:`SOILNOX_SCALE`) to the :option:`SoilNOx` emission
field in the :ref:`Extension Switches <hco-cfg-ext-switches>` section
of :ref:`the HEMCO configuration file <hco-cfg>`:

.. code-block:: kconfig

   104 SOILNOX_ARID    /path/to/soilNOx.climate.generic.05x05.nc  ARID     2000/1/1/0 C xy unitless NO -   1 1
   104 SOILNOX_NONARID /path/to/soilNOx.climate.generic.05x05.nc  NON_ARID 2000/1/1/0 C xy unitless NO -   1 1
   104 SOILNOX_SCALE   1.0                                        -        2000/1/1/0 C xy unitless *  333 1 1

:literal:`SOILNOX_SCALE` is just a dummy scale factor with a global
uniform value of 1.0.   The actual temporal scaling over
the US is done via scale factor :literal:`333` assigned to this
field. This approach ensures that all :option:`SoilNOx` emissions
outside of the US remain intact.

The next step is to define scale factor :literal:`333` (named
:literal:`SOILNOX\_SCALE`) in the :ref:`Scale Factors
<hco-cfg-scalefac>` section of the :ref:`configuration file <hco-cfg>`:

.. code-block:: kconfig

   # Scale factor to scale US soil NOx emissions by a factor of 0.5 in month June-August
   333 SOILNOX_SCALE 1.0/1.0/1.0/1.0/1.0/0.5/0.5/0.5/1.0/1.0/1.0/1.0 - 2000/1-12/1/0 - xy 1 1 5000

Scale factor :literal:`SOILNOX_SCALE` defines a monthly varying scale
factor, with all scale factors being 1.0 except for months
June-August, where the scale factor becomes 0.5. The last column of
the :literal:`SOILNOX_SCALE` entry assigns mask number :literal:`5000`
to this scale factor. This ensures that the scale factor will only be
applied over the region spanned by mask :literal:`5000`. This musk
mast be defined in the :ref:`hco-cfg-masks` section of :ref:`the HEMCO
configuration file <hco-cfg>`:

.. code-block:: kconfig

   1005 USA_MASK      /path/to/usa.mask.nei2005.geos.1x1.nc  MASK 2000/1/1/0 C xy 1 1 -165/10/-40/90
   5000 SOILNOX_MASK   -106.3/37.0/-93.8/49.0                 -    -         - xy 1 1 -106.3/37.0/-93.8/49.0

In this example, mask :literal:`5000` is defined as the region between
106.3 - 93.8 degrees west and 37.0 - 49.0 degrees north. If you want
to apply the soil NOx scaling over the entire US, you can also just
refer to the existing USA mask, e.g.:

.. code-block:: kconfig

   # Scale factor to scale US soil NOx emissions by a factor of 0.5 in month June-August.
   333 SOILNOX_SCALE 1.0/1.0/1.0/1.0/1.0/0.5/0.5/0.5/1.0/1.0/1.0/1.0 - 2000/1-12/1/0 - xy 1 1 1005

.. _cfg-ex-mask:

==================
Mask file examples
==================

Exercise care in defining mask regions
--------------------------------------

In an effort to reduce I/O HEMCO ignores any emission entries that are
deemed "irrelevant" because there is another (global) emission entry
for the same species and emission category (:option:`Cat`), but higher
hierarchy (:option:`Hier`).

For instance, suppose you have the following two fields defined under
:ref:`Base Emissions <hco-cfg-base>`:

.. code-block:: kconfig

    0 TEST_1 file.nc var 2000/1/1/0 C xy 1 1 CO - 1 1
    0 TEST_2 file.nc var 2000/1/1/0 C xy 1 1 CO - 1 2

In this case, during initialization HEMCO determines that
:literal:`TEST_1` is obsolete because it will always be overwritten by
:literal:`TEST_2` because of its higher hierarchy. But if there is a
mask assigned to an emission inventory, HEMCO uses the
provided mask domain to determine whether this inventory has
to be treated as "global" or not.

Going back to the example above, let's add a mask to :literal:`TEST_2`:

.. code-block:: kconfig

   0 TEST_1 file.nc var 2000/1/1/0 C xy 1 1 CO -    1 1
   0 TEST_2 file.nc var 2000/1/1/0 C xy 1 1 CO 1000 1 2

and let´s define the following :ref:`mask <hco-cfg-masks>`:

.. code-block:: kconfig

   1000 TEST_MASK mask.nc var 2000/1/1/0 C xy 1 1 -180/180/-90/90

HEMCO uses the mask range (:literal:`180/180/-90/90`) to define the
extension of this mask. If that range covers the entire HEMCO grid
domain, it considers every emission inventory linked with this mask as
¨global¨. In our example, :literal:`TEST_2` would still be considered
global because the mask extends over the entire globe, and
:literal:`TEST_1` is thus ignored by HEMCO.

However, changing the mask domain to something smaller will tell HEMCO
that :literal:`TEST_2` is not global, and that it cannot drop
:literal:`TEST_1` because of that:

.. code-block:: kconfig

   1000 TEST_MASK mask.nc var 2000/1/1/0 C xy 1 1 -90/180/-45/45

Long story short: if you set the mask range to a domain that is
somewhat smaller than your simulation window, things work just
fine. But if you set the range to something bigger, HEMCO will start
ignoring emission files.

.. _cfg-ex-mask-frac:

Preserve fractional values when masking emissions
-------------------------------------------------

Question from a HEMCO user:

    I see that when the mask files are regridded they are remapped to
    0 or 1 via regular rounding. Unfortunately, this method will not
    work well for my application, because the region I am trying to
    zero out is a small region inside the 4x5 grid cell and thus the
    current mask will not change the emissions on a
    :math:`4^{\circ}{\times}5^{\circ}` scale.

    I was wondering whether it would be possible/straightforward to
    modify the mask regridding method such that
    :math:`4^{\circ}{\times}5^{\circ}` emissions scale will
    scale with the fraction of the gird cell that is masked (e.g., if
    a quarter of the grid cells in one of the
    :math:`4^{\circ}{\times}5^{\circ}` grid are masked, the emissions
    will scale down by 25%).

For this application, it may better to define your mask file in the
:ref:`Scale Factors <hco-cfg-scalefac>` section of :ref:`the HEMCO
configuration file <hco-cfg>`.

By defining a mask in the :ref:`hco-cfg-masks` section, HEMCO
identifies the data container type as MASK and treats the data as
binary.  Long story short:

.. code-block:: kconfig

   ###############################################################################
   ### BEGIN SECTION MASKS
   ###############################################################################

   If your mask file is currently defined here ...

   ### END SECTION MASKS ###

If you instead move that line to the SECTION SCALE FACTORS then HEMCO
will treat the mask as type SCAL. I believe that would preserve the
regridded value (in your example 0.25) and apply that to the emissions
in a 4x5 grid box.

.. code-block:: kconfig

   ###############################################################################
   ### BEGIN SECTION SCALE FACTORS
   ###############################################################################

   ... put your mask file here instead ...

   ### END SECTION SCALE FACTORS ###

.. _cfg-ex-mask-tagged:

Create emissions for geographically tagged species
--------------------------------------------------

.. important::

   Tagging emissions by geographic regions is currently supported only
   for :ref:`base emissions <hco-cfg-base>` but not for emissions
   computed by :ref:`hco-ext`. We hope to add this capability into a
   future HEMCO version.

If you are using HEMCO interfaced to an external model, and need to
create emissions for geographically tagged species, follow thse steps.

#. Define masks for your geographic regions in the :ref:`hco-cfg-masks`
   secton of :ref:`the HEMCO configuration file <hco-cfg>`:

   .. code-block:: kconfig

      #==============================================================================
      # Country/region masks
      #==============================================================================
      1001 MASK_1  -30/30/45/70    - 2000/1/1/0 C xy 1 1 -30/30/45/70
      1002 MASK_2  -118/17/-95/33  - 2000/1/1/0 C xy 1 1 -118/17/-95/33
      1003 MASK_3  my_mask_file.nc - 2000/1/1/0 C xy 1 1 105/-46/160/–10

      # ... etc ...

   If your mask regions are rectangular, you can specify the
   longitude and latitude at the box corners (such as was done for
   :literal:`MASK_1` and :literal:`MASK_2`).  You may also read a mask
   definition from a netCDF file (as was done for :literal:`MASK_3`).

#. In the :ref:`Base Emissions <hco-cfg-base>` section of :ref:`the
   HEMCO configuration file <hco-cfg>`, add extra entries for tagged
   species underneath the entry for the global species, such as:

   .. code-block:: kconfig

      #==============================================================================
      # --- EDGAR v4.2 emissions, various sectors ---
      #==============================================================================
      (((EDGAR

      ### Gas and oil ###
      0 CH4_GAS__1B2a    v42_CH4.0.1x0.1.nc  ch4_1B2a  2004-2008/1/1/0 C xy kg/m2/s CH4   -    1 1
      0 CH4_GAS__1b2a_a  -                   -         -               - -  -       CH4_a 1001 1 1
      0 CH4_GAS__1b2a_b  -                   -         -               - -  -       CH4_b 1002 1 1
      0 CH4_GAS__1b2a_c  -                   -         -               - -  -       CH4_c 1003 1 1
      # ... etc ...

      ### Coal mines ###
      0 CH4_COAL__1B1    v42_CH4.0.1x0.1.nc  ch4_1B1   2004-2008/1/1/0 C xy kg/m2/s CH4   -    2 1
      0 CH4_COAL__1B1_a  -                   -         -               - -  -       CH4_a 1001 2 1
      0 CH4_COAL__1B1_b  -                   -         -               - -  -       CH4_b 1002 2 1
      0 CH4_COAL__1B1_c  -                   -         -               - -  -       CH4_c 1003 2 1
      # ... etc ...``


This will put the total emissions into your CH4 tracer (tracer #1). It
will then also apply the regional masks to the total emissions and
then store them into tagged species (i.e. :literal:`CH4_a`,
:literal:`CH4_b`, and :literal:`CH4_c`).  These tagged species must
also be defined in your external model with the same names.

.. _cfg-ex-ext:

=========================
HEMCO extensions examples
=========================

.. _cfg-ex-ext-fix-megan:

Fix MEGAN extension emissions to a specified year
-------------------------------------------------

Question submitted by a HEMCO user:

   Is it possible to fix :option:`MEGAN` emissions to a given year? I know
   this works for many other :ref:`base emissions <hco-cfg-base>`
   inventories, but MEGAN emissions are dependent on environmental
   variables.

Your best option may be to run the HEMCO standalone and save out
MEGAN emissions for the desired year.  Then, in a subsequent run, you
can read in the :ref:`HEMCO diagnostic output <hco-diag>` files
containing the archived :option:`MEGAN` emissions.

#. Run the HEMCO standalone model. Make sure the following entries
   to your :file:`HEMCO_Diagn.rc` file:

   .. code-block:: kconfig

      EmisISOP_Biogenic  ISOP   108    -1  -1   2   kg/m2/s  ISOP_emissions_from_biogenic_sources
      EmisISOP_Biogenic  ISOP   108    -1  -1   2   kg/m2/s  ISOP_emissions_from_biogenic_sources
      EmisALD2_Biogenic  ALD2   108    -1  -1   2   kg/m2/s  ALD2_emissions_from_biogenic_sources
      # ... etc for other MEGAN species ...

   In the above entries, :literal:`108` tells HEMCO to get the
   emissions from the :option:`MEGAN` extension, which is listed in
   the :ref:`Extension Switches <hco-cfg-ext-switches>` section of the
   :ref:`configuration file <hco-cfg>` with :option:`ExtNr` 108.

#. Add the following lines in the :ref:`Settings <hco-cfg-settings>`
   section of :ref:`the HEMCO configuration file <hco-cfg>`:

   .. code-block:: kconfig

      DiagnFile:                   HEMCO_Diagn.rc
      DiagnPrefix:                 HEMCO_diagnostics
      DiagnFreq:                   Monthly

   For more information, see the sections on :option:`DiagnFile`,
   :option:`DiagnPrefix`, :option:`DiagnFreq`.

#. Turn off the MEGAN extension in the :ref:`Extension Switches
   <hco-cfg-ext-switches>` section of the configuration file.

   .. code-block:: kconfig

      108     MEGAN                  : off   ISOP/ACET/PRPE/...etc additional species...

#. Add entries for reading the fixed MEGAN emission that were archived
   in Step 1 under :ref:`Base Emissions <hco-cfg-base>`.  For example:

   .. code-block:: kconfig

      0 MEGAN_ISOP /path/to/HEMCO_diagnostic.2016$MM010000.nc EmisISOP_Biogenic 2016/1-12/1/1/0 C xy kg/m2/s ISOP - 4 1

   .. note::

      HEMCO category :literal:`Cat = 4` is reserved for biogenic emissions.

#. Run HEMCO in either standalone mode, or coupled to an external
   model, dependingon your application.

.. _cfg-ex-ext-emit-2d-levels:

Add 2D emissions into specific levels
-------------------------------------

HEMCO can emit emissions into a layer other than the surface layer.
For example:

.. code-block:: kconfig

   0 EMEP_CO EMEP.nc CO 2000-2014/1-12/1/0 C xyL5 kg/m2/s CO 1/1001 1 2

will release the :literal:`EMEP_CO` into level 5 instead of
level 1. Theoretically, you could create a separate HEMCO entry for
every emission level (under :ref:`Base Emissions <hco-cfg-base>`:

.. code-block:: kconfig

   0 EMEP_CO_L1 EMEP.nc CO 2000-2014/1-12/1/0 C xyL1 kg/m2/s CO 1 150/1001 1 2
   0 EMEP_CO_L2 EMEP.nc CO 2000-2014/1-12/1/0 C xyL2 kg/m2/s CO 1 151/1001 1 2
   0 EMEP_CO_L3 EMEP.nc CO 2000-2014/1-12/1/0 C xyL3 kg/m2/s CO 1 152/1001 1 2

and assign :ref:`Scale Factors <hco-cfg-scalefac>` (e.g. 150, 151,
152) to specify the fraction of EMEP emissions to be added into each level:

.. code-block:: kconfig

   151 EMEP_LEV1_FRAC 0.5 - - - xy 1 1
   152 EMEP_LEV2_FRAC 0.1 - - - xy 1 1
   153 EMEP_LEV3_FRAC 0.1 - - - xy 1 1``

But this approach is somewhat cumbersome. Also, this won’t give you
the possibility to specifically emit a fraction above the PBL given
that the PBL height is variable over time.

Use this notation (under :ref:`Base Emissions <hco-cfg-base>`) to tell
HEMCO that you would like EMEP emissins to be added into levels 1 through 3:

.. code-block:: kconfig

   0 EMEP_CO_L1 EMEP.nc CO 2000-2014/1-12/1/0 C xyL=1:3 kg/m2/s CO 1 1001 1 2

The emissions are then spread across the lowest 3 model levels based
upon the model level thicknesses.

Instead of specifying the model levels, you may also specify the
altitude in meters or use :literal:`PBL` for the planetary boundary
layer:

.. code-block:: kconfig

   # Emit from surface up to 2500 meters
   0 EMEP_CO_L1 EMEP.nc CO 2000-2014/1-12/1/0 C xyL=1:2500m kg/m2/s C 1001 1 2

   # Emit between 1000 and 5000 meters altitude
   0 EMEP_CO_L1 EMEP.nc CO 2000-2014/1-12/1/0 C xyL=1000m:5000m kg/m2/s CO 1 1001 1 2

   # Emit between 5000 meters altitude and model level 17
   0 EMEP_CO_L1 EMEP.nc CO 2000-2014/1-12/1/0 C xyL=500m:17 kg/m2/s CO 1 1001 1 2

   # Emit from the surface to the PBL top
   0 EMEP_CO_L1 EMEP.nc CO 2000-2014/1-12/1/0 C xyL=1:PBL kg/m2/s CO 1 1001 1 2

HEMCO can also read the emission levvel from an external source
(e.g. netCDF file) that is listed as a scale factor.  This field can
then be referred to using its scale factor ID.  As an example, let's
assume daily varying emission heights for 2009-2010 are archived in
:file:`emis_heights.nc` as variable :literal:`emish` in units of
:literal:`m`. available for years 2009 to 2010). You can then define a
:ref:`Scale Factor <hco-cfg-scalefac>` such as:

.. code-block:: kconfig

   300 EMIT_HEIGHT emis_heights.nc emish 2009-2010/1-12/1-31/0 C xy m 1

and refer to this scale factor as the upper bound of the injection
height under :ref:`Base Emissions <hco-cfg-base>`:

.. code-block:: kconfig

   0 GFAS_CO GFAS_201606.nc cofire 2009-2010/1-12/1-31/0 C xyL=1:scal300 kg/m2/s CO - 5 3

It should be noted that HEMCO always regrids the fields to the model
grid before doing any data operations. If the emission height file is
very spotty and contains a lot of zeros the averaged injection heights
may be too low. In this case it may be required to set all zeros to
missing values (which are ignored by HEMCO) to achieve the desired result.

.. _cfg-ex-ext-fix-vert-dist-2d:

Vertically distributing emissions
---------------------------------

In HEMCO 3.0.0 and later versions, the capability to vertically
allocate emissions has been added. To achieve this, HEMCO first copies
emissions to all levels when dimensions :literal:`xyL*` are specified.
Scale factors can then be applied to determine distribute the
emissions vertically.

For example, let's assume that we have a file :file:`vert_alloc.nc`
containing the ratio of emissions to apply to each level for CEDS
energy, industry, and ship emissions.  We must add the following
entries to under the :ref:`Scale Factors <hco-cfg-scalefac>` section
of the :ref:`the HEMCO configuration file <hco-cfg>`:

.. code-block:: kconfig

   #==============================================================================
   # --- CEDS vertical partitioning ---
   #==============================================================================
   (((CEDS
   315 ENERGY_LEVS   vert_alloc.nc g_energy   2017/1/1/0 C xyz 1 1
   316 INDUSTRY_LEVS vert_alloc.nc g_industry 2017/1/1/0 C xyz 1 1
   317 SHIP_LEVS     vert_alloc.nc cmv_c3     2017/1/1/0 C xyz 1 1
   )))CEDS

These scale factors are then applied to the :literal:`CEDS_*_ENE`,
:literal:`CEDS_*_IND`,  and :literal:`CEDS_*_SHIP` fields that are
listed under :ref:`Base Emissions <hco-cfg-base>`.  These fields are
2D in the CEDS data files, but we now can specify dimensions
:literal:`xyL*` instead of :literal:`xy` to tell HEMCO to copy the
field into each emissions level:

.. code-block:: kconfig

   0 CEDS_CO_ENE CO-em-total-anthro_CEDS_$YYYY.nc  CO_ene  1970-2017/1-12/1/0 C xyL* kg/m2/s CO 26/37/35/315 1  5
   0 CEDS_CO_IND CO-em-total-anthro_CEDS_$YYYY.nc  CO_ind  1970-2017/1-12/1/0 C xyL* kg/m2/s CO 26/316       1  5
   0 CEDS_CO_SHP CO-em-total-anthro_CEDS_$YYYY.nc  CO_shp  1970-2017/1-12/1/0 C xyL*`kg/m2/s CO 26/317       10 5

.. _cfg-ex-other-math:

=================================
Mathematical expressions examples
=================================

You may define mathematical expressions in :ref:`the HEMCO
configuration file <hco-cfg>`.  Similar to uniform values, these must
be placed in in the :option:`sourceFile` column.  All expressions are
evaluated during run-time. They can be used e.g. to model an
oscillating emission source. All mathematical expressions must contain
at least one time-dependent variable that is evaluated
on-the-fly. Mathematical expressions are specified by using the prefix
:literal:`MATH:`, followed by the mathematical expression. The
expression is a combination of variables, mathematical operations, and
constants (e.g. :literal:`MATH:5.0+2.5\*sin(HH)`.

.. _cfg-ex-other-math-vars:

Supported variables and operators
---------------------------------

The following variable names and mathematical operations are currently
supported:

**Variable names**

- :literal:`YYYY` (current year)
- :literal:`MM` (current month)
- :literal:`DD` (current day)
- :literal:`HH` (current hour)
- :literal:`NN` (current minute)
- :literal:`SS` (current second)
- :literal:`SS` (current second)
- :literal:`DOY` (day of year)
- :literal:`DOM` (days in current month)
- :literal:`WD` (Weekday: 0=Sun, 1=Mon .. 7=Sat)
- :literal:`LH` (hour in local time)
- :literal:`PI` (the constant PI)

**Basic mathematical operators:** + - / * ^ ( )

**Advanced mathematical functions**: *sin*, *cos*, *tan*,
*asin*, *acos*, *atan*, *sinh*, *cosh*, *tanh*, *sind*,
*cosd*, *tand*,  *log*, *log10*, *nint*, *anint*,
*aint*, *exp*, *sqrt*, *abs*, *floor*. The names refer to
the equivalent Fortran functions.

.. important::

   When using mathematical expressions, we recommend setting the
   :option:`sourceTime` attribute to :literal:`*`, especially if you
   are using the short-term variables (:literal:`HH`, :literal:`NN`,
   :literal:`SS`, :literal:`LH`).  This will ensure that your
   expression will get evaluated on every emission time step.


.. _cfg-ex-other-math-sine:

Example: Define a sinusoidal source
-----------------------------------

To define a sine-wave emission source of NO with an oscillation
frequency of 24 hours, add the following line to section :ref:`Base
Emissions <hco-cfg-base>` in :ref:`the HEMCO configuration file
<hco-cfg>`.  Place the mathematical expression under the
:option:`sourceFile` column (i.e. the 3rd column):

.. code-block:: kconfig

   0 SINE_NO  MATH:sin(HH/12*PI) - * C xy kg/m2/s NO - 1 500

This defines an emission category (:option:`Cat`) of :literal:`1` and
hierarchy (:option:`Hier`) of :literal:`500`.  No scale factors are
applied.

.. important::

   Mathematical expressions can produce negative emissions, which by
   default cause HEMCO to stop with an error. Negative emissions can
   be enabled by setting :literal:`Negative values: 2` in the
   :ref:`Settings <hco-cfg-settings>` section of :ref:`the HEMCO
   configuration file <hco-cfg>`.

In order to avoid negative values, you may specify an offset, as is
shown below:

.. code-block:: kconfig

   0 SINE_NO  MATH:2.0+sin(HH/12*PI) - * C xy kg/m2/s NO - 1 500

.. _cfg-ex-other:

==============
Other examples
==============

.. _cfg-ex-other-passive:

Assign emissions to passive species in an external model
--------------------------------------------------------

The HEMCO passive species module allows you to run a suite of passive
species alongside any simulation, i.e. it works with all simulation
types. To use the passive species within GEOS-Chem, follow these steps:

Let's assume you are using HEMCO in an external model, and that you
have two passive species named :literal:`PASV1` and :literal:`PASV2`
that have constant emissions fluxes.  Add the following entries to the
:ref:`Base Emissions <hco-cfg-base>` section of :ref:`the HEMCO
configuration file <hco-cfg>`:

.. code-block:: kconfig

   # Assign PASV1 a flux of 0.001 kg/m2/s
   0 PASV1_Flux 1.0e-3  - - - xy kg/m2/s PASV1 - 1 1

   # Assign PASV2 a flux of 1e-9 kg/m2/s
   0 PASV2_Flux 1.0e-9  - - - xy kg/m2/s PASV2 - 1 1

   # ... etc for additional species ...

To define emissions for passive species that are geographically
tagged, simply assign corresponding mask values in the third-to-last
column:

.. code-block:: kconfig

   0 PASV1_Flux 1.0e-3  - - - xy kg/m2/s PASV1 1000 1 1
   0 PASV2_Flux 1.0e-9  - - - xy kg/m2/s PASV2 1001 1 1

   # ... etc for additional species...

Here, 1000 and 1001 refer to :ref:`mask definitions <hco-cfg-masks>`
in :ref:`the HEMCO configuration file <hco-cfg>`.


Next, request HEMCO diagnostic output.  Define the following entries
in the :ref:`diagnostics configuration file <hco-diag-configfile>` (aka
:file:`HEMCO_Diagn.rc`):

.. code-block:: kconfig

   # Name       Spec  ExtNr Cat Hier Dim Unit     Longname
   PASV1_TOTAL  PASV1 -1    -1  -1   2   kg/m2/s  PASV1_emission_flux
   PASV2_TOTAL  PASV2 -1    -1  -1   2   kg/m2/s  PASV2_emission_flux

   # ... etc for additional species ...

To activate these diagnostics, you must specify values for
:option:`DiagnFile` and :option:`DiagnFreq` in the :ref:`Settings
<hco-cfg-settings>` section of :ref:`the HEMCO configuration file
<hco-cfg>`:

.. code-block:: kconfig

   DiagnFile:                 HEMCO_Diagn.rc
   DiagnFreq:                 00000000 003000

The :option:`DiagnFile` option tells HEMCO to read the diagnostic
definitions in the file that you specify (the default is
:file:`HEMCO_Diagn.rc`).  Use :option:`DiagnFreq` to specify the
diagnostic frequency (i.e. the interval at which diagnostics
output will be created).
