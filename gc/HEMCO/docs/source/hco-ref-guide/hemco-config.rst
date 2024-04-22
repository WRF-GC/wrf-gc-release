.. |br| raw:: html

   <br />

.. _hco-cfg:

############################
The HEMCO configuration file
############################

The HEMCO Configuration file is composed of several sections:
:ref:`Settings <hco-cfg-settings>`,
:ref:`Base Emissions <hco-cfg-base>`,
:ref:`Scale Factors, <hco-cfg-scalefac>`, and
:ref:`Masks <hco-cfg-masks>`.

An overview of the structure and key formats of the HEMCO configuration file
can be found in Figure 2 of :cite:t:`Lin_et_al._2021`:

.. image:: ../_static/lin-et-al-2021-fig2.png
   :width: 600

.. _hco-cfg-settings:

========
Settings
========

Parameters and variables used by HEMCO are defined in between these
comment lines:

.. code-block:: kconfig

   ###############################################################################
   ### BEGIN SECTION SETTINGS
   ###############################################################################

   settings go here

   ### END SECTION SETTINGS ###

The order within the settings section is irrelevant. Many of these
settings are optional, and default values will be used if not set.

.. _hco-cfg-settings-general:

General simulation settings
---------------------------

These settings control HEMCO simulation options.


.. option:: ROOT

   Root folder containing emissions inventories and other data to be
   read by HEMCO.

.. option:: METDIR

   Root folder of meteorology data files that are needed for HEMCO
   extensions.  Usually this is a subdirectory of :option:`ROOT`.

.. option:: MODEL

   If present, the :literal:`$MODEL` token will be set to the
   value specified.

   If omitted, this value is determined based on compiler switches.

.. option:: RES

   If present, the :literal:`$RES` token will be set to the value
   specified.

   If omitted, this value is determined based on compiler switches.

.. option:: LogFile

   Path and name of the  output log file (which is typically named
   :file:`HEMCO.log`).  If set to the :option:`Wildcard` character,
   all HEMCO output is written to **stdout** (i.e. the screen).

.. option:: Unit tolerance

   Integer value denoting the tolerance against differences between
   the units set in the :ref:`HEMCO configuration file <hco-cfg>`
   and data units found in the source file.  Allowable values are"

   .. option:: 0

      No tolerance.  A units mismatch will halt a HEMCO simulation.
      mismatch).

   .. option:: 1

      Medium tolerance. A units mismatch will print a warning message
      but not halt a HEMCO simulation.  **(Default setting)**

   .. option:: 2

      High tolerance.  A units mismatch will be ignored.

.. option:: Negative values

   Integer value that defines how negative values are handled.

   .. option:: 0

      No negative values are allowed.  **(Default setting)**

   .. option:: 1

     All negative values are set to zero and a warning is given.

   .. option:: 2

      Negative values are kept as they are.

.. option:: Verbose

   Integer value that controls the amount of additional information
   printed to the HEMCO log file.  Allowable values are :literal:`0`
   (no additional output) to :literal:`3` (lots of additional output).
   Setting  :literal:`3` is useful for debugging.

   **Default setting**: :literal:`0`.

.. option:: Warnings

   Integer value that controls the amount of warnings printed
   to the HEMCO log file.  Allowable values are :literal:`0` (no
   warnings) to :literal:`3` (all warnings).

   **Default setting**: :literal:`1` (only severe warnings).

.. option:: Wildcard

   Wildcard character.  On Unix/Linux, this should be set to :literal:`*`.

.. option:: Separator

   Separator symbol. On Unix/4Linux systems, this should be set to
   :literal:`/`.

.. option:: Mask fractions

   If :literal:`true`, the fractional mask values are taken
   into account. This means that mask values can take any value
   between 0.0 and 1.0.

   If :literal:`false`, masks are binary, and grid boxes are
   100% inside or outside of a mask region.

   **Default setting:** :literal:`false`

.. option:: PBL dry deposition

   If :literal:`true`, it is assumed that dry deposition occurs over
   the entire boundary layer. In this case,  extensions that  include
   loss terms (e.g. air-sea exchange) will calculate a loss term for
   every grid box that is partly within the planetary boundary layer.

   If :literal:`false`, a loss term is calculated for the surface
   layer only.

   **Default setting:** :literal:`false`

.. _hco-cfg-settings-emissions:

Emissions settings
------------------

The following options can be used to hold emissions constant over a
year, month, day, or hour, and to scale emissions to a given value:

.. option:: Emission year

   If present, this emission year will be used regardless of the model
   simulation year.

   If omitted, the emission year will be set to the model simulation
   year.

.. option:: Emission month

   If present, this emission month will be used regardless of the model
   simulation month.

   If omitted, the emission month will be set to the model simulation
   month.

.. option:: Emission day

   If present, this emission day will be used regardless of the model
   simulation day.

   If omitted, the emission day will be set to the model simulation
   day.

.. option:: Emission hour

   If present, this emission month will be used regardless of the model
   simulation hour.

   If omitted, the emisison month will be set to the model simulation
   hour.

.. option:: EmissScale_<species-name>

   Optional argument to define a uniform scale factor that will be
   applied across all inventories, categories, hierarchies, and
   extensions.  Can be set for every species individually, e.g.

   .. code-block:: kconfig

      EmisScale_NO: 1.5
      EmisScale_CO: 2.0

   Scales all NO emissions by 50% and doubles CO emissions.

.. _hco-cfg-settings-diagnostics:

Diagnostics settings
--------------------

The following options control archival of diagnostic quantities.  For
more information about HEMCO diagnostics, please see the
:ref:`hco-diag` section.

.. option:: DiagnFile

   Specifies the configuration file for the HEMCO default diagnostics
   collection. This is usually named :file:`HEMCO_Diagn.rc`.  This
   file contains a list of fields to be added to the default
   collection.

   Each line of the diagnostics definition file
   represents a diagnostics container. It expects the following 7 entries
   (all on the same line):

   #. Container name (character)
   #. HEMCO species (character)
   #. Extension number (integer)
   #. Emission category (integer)
   #. Emission hierarchy (integer)
   #. Space dimension (2 or 3)
   #. Output unit (character)
   #. Long name of diagnostic, for the netCDF :literal:`long_name`
      variable attribute (character)

   .. note::

      If you are not sure what the container name, extension number,
      category, and hierarchy are for a given diagnostic, you can set
      :literal:`Verbose` to 3 in the HEMCO configuration file, and run a
      very short simulation (a couple of model hours). Then you can look
      at the output in the :file:`HEMCO.log` file to determine what these
      values should be.

   Please see the :ref:`Default diagnostics collection
   <hco-diag-default>` section for more information about the
   configuration file (e.g. :file:`HEMCO_Diagn.rc`).

.. option:: DiagnFreq

   This setting (located in the HEMCO configuration file) specifies
   the output frequency of the :ref:`Default  <hco-diag-default>`
   collection.  Allowable values are:

   .. option:: Always

      Archives diagnostics on each time step.

   .. option:: Hourly

      Sets the diagnostic time period to 1 hour.

   .. option:: Daily

      Sets the diagnostic time period to 1 day.

   .. option:: Monthly

      Sets the diagnostic time period to 1 hour.

   .. option:: Annually

      Sets the diagnostic time period to 1 year.

   .. option:: End

      Sets the diagnostic time period so that output will only happen
      at the end of the simulation.

   .. option:: YYYYMMDD hhmnss

      Sets the diagnostic time period to an interval specified by a
      15-digit string with year-month-day, hour-minute-second.  For
      example:

      - :code:`00010000 000000` will generate diagnostic output once
	per year.
      - :code:`00000001 000000` will generate diagnostic output once
	per day.
      - :code:`00000000 020000` will generate diagnostic output every
	2 hours.
      - etc.

.. option:: DiagnPrefix

   Specifies the name of the diagnostic files to be created.  For
   example:

   .. code-block::

      DiagnPrefix: ./OutputDir/HEMCO_diagnostics

   will create HEMCO diagnostics files in the :file:`OutputDir/`
   subdirectory of the run directory, and all files will begin with
   the text :file:`HEMCO_diagnostics`.

.. option:: DiagnRefTime

   This option must be explicity added to the HEMCO configuration
   file.

   By default, the value of the :literal:`time:units` attribute in the
   :file:`HEMCO_diagnostics.*.nc` files will be  :literal:`hours since
   YYYY-MM-DD hh:mn:ss`, where :literal:`YYYY-MM-DD hh:mn:ss` is the
   diagnostics datetime.  This default value can be overridden and set
   to a fixed datetime by setting :option:`DiagnRefTime` in the HEMCO
   configuration file.  For example:

   .. code-block:: console

      DiagnRefTime: hours since 1985-01-01 00:00:00

   will set the :literal:`time:units` attribute to :literal:`hours since
   1985-01-01 00:00:00`.

.. option:: DiagNoLevDim

   This option must be explicity added to the HEMCO configuration
   file. If omitted, the default behavior will be :code:`false`.

   If :literal:`true`, the created :file:`HEMCO_diagnostics*.nc` files
   will contain dimensions :literal:`(time,lat,lon)`.     But if at least
   one of the diagnostic quantities has a :literal:`lev` dimension,
   then the created files will have :literal:`(time,lev,lat,lon)`
   dimensions.

   If :literal:`false`, the :file:`HEMCO_diagnostics.*.nc` files will
   always contain dimensions :literal:`(time,lev,lat,lon)`.

.. option:: DiagnTimeStamp

   This option must be explicity added to the HEMCO configuration
   file.  If omitted, the default behavior will be :option:`End`.

   Allowable values are:

   .. option:: End

      Uses the date and time at the end of the diagnostics time window
      to timestamp diagnostic files. With this option, a 1-hour
      simulation from :code:`20220101 000000` to :code:`20220101
      010000` will create a diagnostic file named
      :file:`HEMCO_Diagnostics.202201010100.nc`.

   .. option:: Start

      Uses the date and time at the start of the diagnostics time
      window to timestamp diagnostic files.  With this option, a
      1-hour simulation from :code:`20220101 000000` to
      :code:`20220101 010000` will create a diagnostic file named
      :file:`HEMCO_Diagnostics.202201010000.nc`.

   .. option:: Mid

      Uses the date and time at the midpoint of the diagnostics time
      window to timestamp diagnostic files. With this option, a 1-hour
      simulation from :code:`20220101  000000` to :code:`20220101
      010000` will create a diagnostic file named
      :file:`HEMCO_diagnostics.202201010030.nc`.

.. _hco-cfg-settings-standalone:

HEMCO standalone simulation settings
------------------------------------

In standalone mode, the three simulation description files also need be
specified:

.. option:: GridFile

   Path and name of the grid description file, which is usually named
   :file:`HEMCO_sa_Grid.rc`.

.. option:: SpecFile

   Path and name of the species description file, which is usually named
   :file:`HEMCO_sa_Spec.rc`.

.. option:: GridFile

   Path and name of the time description file, which is usually named
   :file:`HEMCO_sa_Time.rc`.

.. _hco-cfg-settings-usrdef:

User-defined tokens
-------------------

Users can specify any additional token in the **Settings** section
section. The token name/value pair must be separated by the colon (:)
sign. For example, adding the following line to the settings section
would register token :literal:`$ENS` (and assign value 3 to it):

.. code-block:: kconfig

   ENS: 3

User-defined tokens can be used the same way as the built-in tokens
(:literal:`$ROOT`, :literal:`$RES`, :literal:`YYYY`, etc.). See
:literal:`sourceFile` in the Base emissions for more details about
tokens.

.. important::

    User-defined token names must not contain numbers or
    special characters such as :literal:`.`, :literal:`_`,
    :literal:`-`, or :literal:`x`.

.. _hco-cfg-ext-switches:

==================
Extension switches
==================

HEMCO performs automatic emission calculations using all fields that
belong to the :ref:`base emisisons extension <hco-cfg-base>`. Additional
emissions that depend on environmental parameter such as wind speed or
air temperature--and/or that use non-linear parameterizations--are
calculated through :ref:`hco-ext`.  A list of currently implemented
extensions in HEMCO is given in `Keller et al. (2014)
<#References>`__. To add new extensions to HEMCO, modifications of the
source code are required, as described further in :ref:`hco-hood`.

The first section of the configuration file lists all available
extensions and whether they shall be used or not. For each extension,
the following attributes need to be specified:

.. option:: ExtNr

   Extension number associated with this field. All
   :ref:`base emissions <hco-cfg-base>` should have extension number
   zero.  The extension number` of the data listed in section
   :ref:`hco-ext` data must match with the corresponding extension
   number.

   The extension number can be set to the wildcard character. In that
   case, the field is read by HEMCO (if the assigned species name
   matches any of the HEMCO species, see :option:`Species` below) but
   not used for emission calculation. This is particularly useful if
   HEMCO is only used for data I/O but not for emission calculation.

.. option:: ExtName

   Extension name.

.. option:: Toggle

   If :literal:`on`, the extension will be used.

   If :literal:`off`, the extension will not be used.

.. option:: Species

   List of species to be used by this extension. Multiple species are
   separated by the :option:`Separator` symbol
   (e.g. :literal:`/`). All listed species must be supported by the
   given extension.

   - For example, the soil NO emissions extension only supports one
     species (NO). An error will be raised if additional species are
     listed.

Additional extension-specific settings can also be specified in the
'Extensions Settings' section (see also an example in
:ref:`edit-hco-cfg` and the definition of
:ref:`hco-cfg-data-coll`. These settings must immediately follow the
extension definition.

HEMCO expects an extension with extension number zero, denoted the
:ref:`base emisisons extension <hco-cfg-base>` extension. All emission
fields linked to the base extension will be used for automatic
emission calculation. Fields assigned to any other extension number
will not be inlcuded in the base emissions calculation, but they are
still read/regridded by HEMCO (and can be made available readily
anywhere in the model code). These data are only read if the
corresponding extension is enabled.

All species to be used by HEMCO must be listed in column
:option:`Species` of the base extension switch. In particular, all
species used by any of the other extensions must also be listed as
base species, otherwise they will not be recognized. It is possible
(and recommended) to use the :option:`Wildcard` character, in which
case HEMCO automatically determines what species to use by matching
the atmospheric model species names with the species names assigned to
the base emission fields and/or any emission extension.

The environmental fields (wind speed, temperature, etc.) required by the
extensions are either passed from the atmospheric model or read through
the HEMCO configuration file, as described in :ref:`hco-ext`.

.. _hco-cfg-base:

==============
Base emissions
==============

The BASE EMISSIONS section lists all base emission fields and how they
are linked to :ref:`scale factors <hco-cfg-scalefac>`. Base emissions
settings must be included between these comment lines:

.. code-block:: kconfig

   ###############################################################################
   ### BEGIN SECTION BASE EMISSIONS
   ###############################################################################
   settings go here

   ### END SECTION BASE EMISSIONS ###

The :option:`ExtNr` field is defined in :ref:`hco-cfg-ext-switches`.

Other attributes that need to be defined for each base emissions entry
are:

.. option:: Name

   Descriptive field identification name. Two consecutive underscore
   characters (:literal:`__`) can be used to attach a 'tag' to a
   name. This is only of relevance if multiple base emission fields
   share the same species, category, hierarchy, and scale factors. In
   this case, emission calculation can be optimized by assigning field
   names that onlydiffer by its tag to those fields
   (e.g. :literal:`DATA__SECTOR1`, :literal:`DATA__SECTOR2`, etc.).

   For fields assigned to extensions other than the base extension
   (:literal:`ExtNr = 0`), the field names are prescribed and must not
   be modified because the data is identified by these extensions by
   name.

.. option:: sourceFile

   Path and name of the input file.

   Name tokens can be provided that become evaluated during
   runtime. For example, to use the root directory specified in the
   :ref:`Section settings section  <hco-cfg-settings>`, the
   :literal:`$ROOT` token can be used.  Similarly the token
   :literal:`$CFDIR` refers to the location of the configuration
   file. This allows users to reference data relative to the
   location of the configuration file. For instance, if the
   data  is located in subfolder :literal:`data` of the same directory
   as the configuration file, the file name can be set to
   :literal:`$CFDIR/data/filename.nc`.

   Similarly, the **date tokens** :literal:`$YYYY`, :literal:`$MM`,
   :literal:`$DD`, :literal:`$HH`, and :literal:`$MN` can be used to
   refer to the the current valid year, month, day, hour, and
   minute, respectively. These values are determined
   from the current simulation datetime and the :option:`sourceTime`
   specification for this entry.

   The tokens :literal:`$MODEL` and :literal:`$RES` refer to the
   meteorological model (:option:`MODEL`) and resolution
   (:option:`RES`). These tokens can be set explicitly in the settings
   section. In `GEOS-Chem <https://geos-chem.readthedocs.io>`_ they
   are set to compiler-flag specific values if not set in the settings
   section.  Any token defined in the settings section can be used to
   construct a part of the file name (see :ref:`hco-cfg-settings-usrdef`).

   As an alternative to an input file, **geospatial uniform values**
   can directly be specified in the configuration file (see e.g. scale
   factor :literal:`SO2toSO4` in :ref:`edit-hco-cfg`). If multiple
   values are provided (separated by the separator    character), they
   are interpreted as different time slices. In this case, the
   :option:`sourceTime` attribute can be used to specify the times
   associated with the individual slices. If no time attribute is set,
   HEMCO attempts to determine the time slices from the number of data
   values: 7 values are interpreted as weekday (Sun, Mon, ..., Sat); 12
   values as month (Jan, ..., Dec); 24 values as hour-of-day (12am,
   1am, ..., 11pm).

   Uniform values can be combined with **mathematical expressions**,
   e.g. to model a sine-wave emission source. Mathematical
   expressions must be labeled :literal:`MATH:`, followed by
   the expression, e.g. :literal:`MATH:2.0+sin(HH/12*PI)`.

   **Country-specific data** can be provided through an ASCII file
   (:literal:`.txt`). More details on this option are given in the
   Input File Format section.

   If this entry is **left empty** (:literal:`-`), the filename from
   the preceding entry is taken, and the next 5 attributes will be
   ignored (see entry :literal:`MACCITY_SO4` in :ref:`edit-hco-cfg`.

.. option:: sourceVar

   Source file variable of interest. Leave empty (:literal:`-`) if
   values are directly set through the :option:`sourceFile` attribute
   or if :option:`sourceFile` is empty.

.. option:: sourceTime

   This attribute defines the time slices to be used and the data
   refresh frequency. The format is
   :literal:`year/month/day/hour`. Accepted are discrete dates for
   time-independent data (e.g. :literal:`2000/1/1/0`) and time ranges
   for temporally changing fields
   (e.g. :literal:`1980-2007/1-12/1-31/0-23`). Data will automatically
   become updated as soon as the simulation date enters a new time
   interval.

   The provided time attribute determines the data refresh
   frequency. It does not need to correspond to the datetimes of the
   input file.

   - For example, if the input file contains daily data of
     year 2005 and the time attribute is set to :literal:`2005/1/1/0`,
     the file  will be read just once (at the beginning  of the
     simulation) and the data of Jan 1, 2005 is used throughout the
     simulation. |br|
     |br|
   - If the time attribute is set to :literal:`2005/1-12/1/0`, the
     data is updated on every month, using the first day data of the
     given month. For instance, if the simulation starts on July 15,
     the data of July 1,2005 are used until August 1, at which point
     the  data will be refreshed to values from August 1, 2005. |br|
     |br|
   - A time attribute of :literal:`2005/1-12/1-31/0` will make
     sure that the input data are refreshed daily to the current day's
     data. |br|
     |br|
   - Finally, if the time attribute is set to
     :literal:`2005/1-12/1-31/0-23`, the data file is read every
     simulation hour, but the same daily data is used throughout the
     day (since there are no hourly data in the file). Providing too
     high update frequencies is not recommended unless the data
     interpolation option is enabled (see below).

   If the provided time attributes do not match a datetime of the
   input file, the **most likely** time slice is selected. The most
   likely time slice is determined based on the specified source time
   attribute, the datetimes available in the input file, and the
   current simulation date. In most cases, this is just the closest
   available time slice that lies in the past.

   - For example, if a file contains annual data from 2005 to 2010 and
     the source time attribute is set to :literal:`2005-2010/1-12/1/0`,
     the data of 2005 is used for all simulation months in 2005. |br|
     |br|
   - More complex datetime selections occur for files with
     discontinuous time slices, e.g. a file with monthly data for
     year 2005, 2010, 2020, and 2050. In this case, if the time
     attribute is set to :literal:`2005-2020/1-12/1/0`, the monthly
     values of 2005 are (re-)used for all years between 2005 and 2010,
     the monthly values of 2010 are used for simulation years 2010 -
     2020, etc.

   It is possible to use tokens :literal:`$YYYY`, :literal:`$MM`,
   :literal:`$DD`, and :literal:`$HH`, which will automatically be
   replaced by the current simulation date. Weekly data (e.g. data
   changing by the day of the week) can be indicated by setting the
   day attribute to :literal:`WD` (the wildcard character will work,
   too, but is not recommended). Weekly data needs to consist of at
   least seven time slices - in increments of one day - representing
   data for every weekday starting on Sunday. It is possible to store
   multiple weekly data, e.g. for every month of a year:
   :literal:`2000/1-12/WD/0`. These data must contain  time slices for
   the first seven days of every month, with the first day per month
   representing Sunday data, then followed by Monday,
   etc. (irrespective of the real weekdays of the given month). If the
   wildcard character is used for the days, the data will be
   interpreted if (and only if) there are exactly seven time
   slices. See the Input File Format section for more details. Default
   behavior is to interpret weekly data as 'local time', i.e. token
   :literal:`WD` assumes that the provided values are in local
   time. It is possible to use weekly data referenced to UTC time
   using token :literal:`UTCWD`.

   Similar to the weekday option, there is an option to indicate
   hourly data that represents local time: :literal:`LH`. If using
   this flag, all hourly data of a given time interval (day, month,
   year) are read into memory and the local hour is picked at every
   location. A downside of this is that all hourly time slices in
   memory are updated based on UTC time. For instance, if a file holds
   local hourly data for every day of the year, the source time
   attribute can be set to :literal:`2011/1-12/1-31/LH`. On every new
   day (according to UTC time), this will read all 24 hourly time
   slices of that UTC day and use those hourly data for the next 24
   hours. For the US, for instance, this results in the wrong daily
   data being used for the last 6-9 hours of the day (when UTC time is
   one day ahead of local US time).

   There is a difference between source time attributes
   :literal:`2005-2008/$MM/1/0` and :literal:`2005-2008/1-12/1/0`. In
   the first case, the file will be updated annually, while the update
   frequency is monthly in the second case. The token :literal:`$MM`
   simply indicates that the current simulation month shall be used
   whenever the file is updated, but it doesn’t imply a refresh
   interval. Thus, if the source time attribute is set to
   :literal:`$YYYY/$MM/$DD/$HH`, the file will be read only once and
   the data of the simulation start date is taken (and used throughout
   the simulation). For uniform values directly set in the
   configuration file, all time attributes but one must be fixed,
   e.g. valid entries are :literal:`1990-2007/1/1/0` or
   :literal:`2000/1-12/1/1`, but not :literal:`1990-2007/1-12/1/1`.

   .. note::

      All data read from netCDF file are assumed to be in UTC time,
      except for weekday data that are always assumed to be in local
      time. Data read from the configuration file and/or from ASCII are
      always assumed to be in local time.

   It is legal to keep different time slices in different files,
   e.g. monthly data of multiple years can be stored in files
   :file:`file_200501.nc`, :file:`file_200502.nc`, ...,
   :file:`file_200712.nc`.  By setting the source file attribute to
   :file:`file_$YYYY$MM.nc` and the source time attribute to
   :file:`2005-2007/1-12/1/0`, data of :file:`file_200501.nc` is used
   for simulation dates of January 2005 (or any January of a previous
   year), etc. The individual files can also contain only a subset of
   the provided data range, e.g. all monthly files of a year can be
   stored in one file: :file:`file_2005.nc`, :file:`file_2006.nc`,
   :file:`file_2007.nc`. In this case, the source file name should be
   set to :file:`file_$YYYY`, but the source time attribute should
   still be :literal:`2005-2007/1-12/1/0` to indicate that the field
   shall be updated monthly.

   This attribute can be set to the wildcard character (:literal:`*`), which
   will force the file to be updated on every HEMCO time step.

   File reference time can be shifted by a fixed amount by adding an
   optional fifth element to the time stamp attribute. For instance,
   consider the case where 3-hourly averages are provided in
   individual files with centered time stamps, e.g.:
   :file:`file.yyyymmdd_0130z.nc`, :file:`file.yyyymmdd_0430z.nc`,
   ..., :file:`file.yyymmdd_2230z.nc`. To read these files **at the
   beginning of their time intervals**, the time stamp can be shifted by
   90 minutes: :literal:`2000-2016/1-12/1-31/0-23/+90minutes`.  At
   time 00z, HEMCO will then read file 0130z and keep using this file
   until 03z, when it switches to file 0430z. Similarly, it is
   possible to shift the file reference time by any number of years,
   months, days, or hours. Time shifts can be forward or backward in
   time (use :literal:`-` sign to shift backwards).

.. option:: CRE

   Controls the time slice selection if the simulation date is outside
   the range provided in attribute source time (see above). The
   following options are available:

   .. option:: C

      **Cycling:**  Data are interpreted asclimatology and recycled
      once the end of the last time slice is reached. For instance, if
      the input data contains monthly data of year 2000, and the
      source time attribute is set to :literal:`2000/1-12/1/0 C`, the
      same monthly data will be re-used every year.

      If the input data spans multiple years (e.g. monthly data from
      2000-2003), the closest available year will be used outside of
      the available range (e.g. the monthly data of 2003 is used for
      all simulation years after 2003).

   .. option:: CS

      **Cycling, Skip:** Data are interpreted as climatology and recycled
      once the end of the last time slice is reached. Data that aren't
      found are skipped. This is useful when certain fields aren't found
      in a restart file and, in that case, those fields will be
      initialized to default values.

   .. option:: CY

      **Cycling, Use Simulation Year:**, Same as :option:`C`, except
      don't allow :envvar:`Emission year` setting to override year value.

   .. option:: CYS

      **Cycling, Use Simulation Year, Skip:**  Same as :option:`CS`,
      except don't allow :envvar:`Emission year` setting to override year
      value.

   .. option:: R

      **Range:** Data are only considered as long as the simulation
      time is within the time range specified in attribute :option:`sourceTime`.
      The provided range does not necessarily need to match the time
      stamps of the input file. If it is outside of the range of the
      netCDF time stamps, the closest available date will be used.

      For instance, if a file contains data for years 2003 to 2010 and
      the  provided range is set to :literal:`2006-2010/1/1/0 R`, the file
      will only be considered between simulation years 2006-2010. For
      simulation years 2006 through 2009, the corresponding field on
      the file is used. For all years beyond 2009, data of year 2010
      is used. If the simulation date is outside the provided time
      range, the data is ignored but HEMCO does not return an error -
      the field is simply treated as empty (a corresponding warning is
      issued in the HEMCO log file).

      - Example: if the source time attribute is set to
        :literal:`2000-2002/1-12/1/0 R`, the data will be used for
        simulation years 2000 to 2002 and ignored  for all other years.

   .. option:: RA

      **Range, Averaging Otherwise:** Combination of flags :option:`R`
      and :option:`A`. As long as the simulation year is within the
      specified year range, HEMCO will use just the data from that
      particular year. As soon as the simulation year is outside the
      specified year range, HEMCO will use the data averaged over the
      specified years.

      - Consider the case where the emission file contains
        monthly data for years 2005-2010. Setting the time attribute to
        :literal:`2005-2010/1-12/1/0 R` will ensure that this data is
        only used within simulation years 2005 to 2010 and ignored
	outside of it. |br|
	|br|

      - When setting the time attribute to
	:literal:`2005-2010/1-12/1/0 A`, HEMCO will always use the
        2005-2010 averaged monthly values, even for simulation years 2005
        to 2010.

      - A time attribute of :literal:`2005-2010/1-12/1/0 RA` will make
	sure that HEMCO uses the monthly data of the current year if
	the simulation  year is between 2005 and  2010, and the
	2005-2010 average for simulation  years before and after 2005
	and 2010,  respectively.

   .. option:: RF

      **Range, Forced:**  Same as ``R``, but HEMCO stops with an error
      if the simulation date is outside the provided range.

   .. option:: RY

      **Range, Use Simulation Year:** Same as :option:`R`, except
      don't allow :envvar:`Emission year` to override year value.

   .. option:: E

      **Exact:**  Fields are only used if the time stamp on the field
      exactly matches the current simulation datetime. In all other
      cases, data is ignored but HEMCO does not return an error.

      - For example, if :option:`sourceTime` is set to
        :literal:`2000-2013/1-12/1-31/0 E`, every time the simulation
        enters a new day HEMCO will attempt to find a data field for
	the current simulation date.  If no such field can be found on
	the file, the data is ignored (and a warning is
	prompted). This setting is particularly useful for data that
	is highly sensitive to date and time, e.g. restart variables.

   .. option:: EF

      **Exact, Forced:** Same as :option:`E`, but HEMCO stops with an
      error if no data field can be found for the current simulation
      date and time.

   .. option:: EC

      **Exact, Read/Query Contiuously.**.

   .. option:: ECF

      **Exact, Read/Query Continuously, Forced.**

   .. option:: EFYO

      **Exact, Forced, Simulation Year, Once**: Same as :option:`EF`,
      with the following additions:

      - :envvar:`Y`: HEMCO will stop thie simulation if the simulation
	year does not match the year in the file timestamp.
      - :envvar:`O`: HEMCO will only read the file once.

      This setting is typically only used for model restart files
      (such as `GEOS-Chem Classic restart files
      <https://geos-chem.readthedocs.io/en/stable/gcc-guide/04-data/restart-files-gc.html>`_).
      This ensures that the simulation will stop unless the restart
      file timestamp matches the simulation start date and time.

      .. attention::

	 Consider changing the time cycle flag from :option:`EFYO` to
	 :option:`CYS` if you would like your simulation to read a
	 data file (such as a simulation restart file) whose file
	 timestamp differs from the simulaton start date and time.

   .. option:: EY

      **Exact, Use Smulation Year:** Same as :option:`E`, except don't
      allow :envvar:`Emission year` setting to override year value.

   .. option:: A

      **Averaging:** Tells HEMCO to average the data over the
      specified range of years.

      - For instance, setting :option:`sourceTime` to
	:literal:`1990-2010/1-12/1/0 A` will cause HEMCO to calculate
        monthly means between 1990 to 2010 and use these regardless of
	the current simulation date.

      The data from the different years can be spread out over multiple
      files. For example, it is legal to use the averaging flag in
      combination with files that use year tokens such as
      :literal:`file_$YYYY.nc`.

   .. option:: I

      **Interpolation:** Data fields are interpolated in time. As an
      example, let's assume a file contains annual data for years
      2005, 2010, 2020, and 2050. If :option:`sourceTime` is set to
      :literal:`2005-2050/1/1/0 I`, data becomes interpolated between
      the two closest years every time we enter a new simulation
      year. If the simulation starts on January 2004, he value of 2005
      is used for years 2004 and 2005. At the beginning of 2006, the
      used data is calculated as a weighted mean for the 2005 and 2010
      data, with 0.8 weight given to 2005 and 0.2 weight given to 2010
      values. Once the simulation year changes to 2007, the weights
      hange to 0.6 for 2005 and 0.4 for 2010, etc. The interpolation
      frequency is determined by :option:`sourceTime` the source time
      attribute.

      For example, setting the source time attribute to
      :literal:`2005-2050/1-12/1/0 I` would result in a recalculation
      of the weights on every new simulation month. Interpolation
      works in a very similar manner for discontinuous monthly,daily,
      and hourly data. For instance if a file contains monthly data of
      2005, 2010, 2020, and 2050 and the source time attribute is set
      to :literal:`2005-2050/1-12/1/0 I`, the field is recalculated
      every month using the two bracketing fields of the given month:
      July 2007 values are calculated from July 2005 and July 2010
      data (with weights of 0.6 and 0.4, respectively), etc.

      Data interpolation also works between multiple files. For
      instance, if monthly data are stored in files
      :literal`file_200501.nc`, :file:`file_200502.nc`, etc., a
      combination of source file name :file:`file_$YYYY$MM.nc` and
      :option:`sourceTime` attribute :literal:`2005-2007/1-12/1-31/0
      :literal:I` will result in daily data interpolation between the two
      bracketing files, e.g. if the simulation day is July 15, 2005,
      the fields current values are calculated from files
      :literal:`file_200507.nc` and :literal:`file_200508.nc`,
      respectively.

      Data interpolation across multiple files also works if there are
      file 'gaps', for example if there is a file only every three
      hours: :file:`file_20120101_0000.nc`,
      :file:`file_20120101_0300.nc`,  etc. Hourly data interpolation
      between those files can be achieved by setting source file to
      :file:file_$YYYY$MM$DD_$HH00.nc`, and :option:`sourceTime` to
      :literal:`2000-2015/1-12/1-31/0-23 I` (or whatever the covered
      year range is).

.. option:: SrcDim

    Spatial dimension of input data (:literal:`xy` for horizontal
    data; :literal:`xyz` for 3-dimensional data).

    The :option:`SrcDim` attribute accepts an integer number as
    vertical coordinate to indicate the number of vertical levels to
    be read, as well as the direction of the vertical axis. For
    example, to use the lowest 5 levels of the input data only, set
    :option:`SrcDim` to :literal:`xy5`.   This will place the lowest 5
    levels of the input  data into HEMCO levels 1 to 5. To use the
    topmost  5 levels of the input data, set :option:`SrcDim` to
    :literal:`xy-5`. The minus sign will force the vertical  axis to
    be flipped, i.e. the 5 topmost levels will be placed into HEMCO
    levels 1 to 5 (in reversed order, so that the topmost level of the
    input data will be placed in HEMCO lev el 1, etc.).

    The :option:`SrcDim` attribute can also be used to indicate the
    level into which 2D data shall be released by setting  the
    vertical coordinate to `:literal:`LX```, with :literal:`X` being
    the release level. For instance, to emit a 2D field into level 5,
    set :option:`SrcDim` to :literal:`xyL5`.

    HEMCO can has two options to specify the emission injection
    height:

    #. The vertical height can be given as model level (default) or in
       meters, e.g. to emit a source at 2000m:
       :literal:`xyL=2000m`. |br|
       |br|

    #. For 2D fields it is legal to define a range of levels, in which
       case the  emissions are uniformly distributed across these
       levels (maintaining the original total emissions).
       Examples for this are:

       - :literal:`xyL=1:5`: Emit into levels 1-5;
       - :literal:`xyL=2:5000m` Emit between model level 2 and 5000m;
       - :literal:`xyL=1:PBL`: Emit from the surface up to the PBL top.

    HEMCO can also get the injection height information from an
    external source (i.e. netCDF file). For now, these heights are
    expected to be in meters. The injection height data must be
    listed as a scale factor and can then be referenced in the
    :option:`SrcDim` setting.

    HEMCO can read read netCDF files with an arbitrary additional
    dimension. For these files, the name of the additional dimension
    and  the desired dimension index must be specified as part of the
    :option:`SrcDim` attribute.

    - For example, to read a file that contains 3D ensemble data
      (with the individual ensemble runs as additional dimension
      :literal:`ensemble`), set :option:`SrcDim` to
      :literal:`xyz+"ensemble=3`  to indicate that you wish to read
      the third ensemble member. You may also use a
      :ref:`user-defined token <hco-cfg-settings-usrdef>` for the
      dimension index to be used, e.g. :literal:`xyz+"ensemble=$ENS"`.

    .. note::

       Arbitrary additional dimensions are currently not supported in
       a high-performance environment that uses the ESMF/MAPL
       input/output libraries.

.. option:: SrcUnit

   Units of the data.

.. option:: Species

   HEMCO emission species name. Emissions will be added to this
   species. All HEMCO emission species are defined at the beginning of
   the simulation (see the Interfaces section) If the species name
   does not match any of the HEMCO species, the field is ignored
   altogether.

   The species name can be set to the wildcard character, in which
   case the field is always read by HEMCO but no species is assigned
   to it. This can be useful for extensions that import some
   (species-independent) fields by name.

The three entries below only take effect for fields that are assigned
to the base extension (:literal:`ExtNr = 0`), e.g. that are used for
automatic  emission calculation. They are used by HEMCO to determine
how the final emission fields are assembled from all provided data fields.

.. option:: ScalIDs

   Identification numbers of all scale factors and masks that shall be
   applied to this base emission field. Multiple entries must be
   separated by the separator character. The :option:`ScalIDs` must
   correspond to the numbers provided in the :ref:`hco-cfg-scalefac`
   and :ref:`hco-cfg-masks` sections.

.. option:: Cat

   Emission category. Used to distinguish different, independent
   emission sources. Emissions of different categories are always
   added to each other.

   Up to three emission categories can be assigned to each entry
   (separated by the separator character).  Emissions are always
   entirely written into the first listed category, while emissions of
   zero are used for any other assigned category.

   In practice, the only time when more than one emissions category
   needs to be specified is when an :ref:`inventory does not separate
   between anthropogenic, biofuels, and/or trash emissions
   <edit-hco-cfg-ex6>`

   For example, the CEDS inventory uses categories :literal:`1/2/12`
   because CEDS lumps both biofuel emissions and trash emissions with
   anthropogenic Because. The :literal:`1/2/12` category designation
   means "Put everything into the first listed category
   (1=anthropogenic), and set the other listed categories (2=biofuels,
   12=trash) to zero.

.. option:: Hier

   Emission hierarchy. Used to prioritize emission fields within the
   same emission category.  Emissions of higher hierarchy overwrite
   lower hierarchy data. Fields are only considered within their
   defined domain, i.e. regional inventories are only considered
   within their mask boundaries.


.. _hco-cfg-scalefac:

=============
Scale factors
=============

The SCALE FACTORS section of the configuration file lists all scale
factors applied to the base emission field. Scale factors that are not
used by any of the base emission fields are ignored. Scale factors can
represent:

#. Temporal emission variations including diurnal, seasonal, or
   interannual variability;
#. Regional masks that restrict the applicability of the base inventory
   to a given region; or
#. Species-specific scale factors, e.g., to split lumped organic
   compound emissions into individual species.

This sample snippet of the HEMCO configuration file shows how scale
factors can either be read from a netCDF file or listed as a set of
values.

.. code-block:: kconfig

   ###############################################################################
   ### BEGIN SECTION SCALE FACTORS
   ###############################################################################
   # ScalID Name srcFile srcVar srcTime CRE Dim Unit Oper

   # %%% Hourly factors, read from disk %%%
   1 HOURLY_SCALFACT hourly.nc                                      factor 2000/1/1/0-23 C xy 1 1

   # %%% Scaling SO2 to SO4 (molar ratio) %%%
   2 SO2toSO4        0.031                                          -      -             - -  1 1

   # %%% Daily scale factors, list 7 entries %%%
   20 GEIA_DOW_NOX   0.784/1.0706/1.0706/1.0706/1.0706/1.0706/0.863 -      -             - xy 1 1

   ### END SECTION SCALE FACTORS ###

Options :option:`sourceFile`, :option:`sourceVar`,
:option:`sourceTime`, :option:`CRE`, :option:`SrcDim`, and
:option:`SrcUnit`, are described in :ref:`hco-cfg-base`.

Other scale factor options not previously described are:

Scale factor options not previously described are:

.. option:: ScalID

   Scale factor identification number. Used to link the scale factors
   to the base emissions through the corresponding ScalIDs attribute
   in the :ref`hco-cfg-base`.

.. option:: Oper

   Scale factor operator. Determines the operation performed on the
   scale factor.  Possible values are:

   -  1 for multiplication (Emission = Base \* Scale);
   -  -1 for division (Emission = Base / Scale);
   -  2 for squared (Emission = Base \* Scale**2).

.. option:: MaskID

   Optional.  ScalID of a mask field. This optional value can be used
   if a scale factor shall only be used over a given region. The
   provided MaskID must have a corresponding entry in the
   :ref:`Masks section <hco-cfg-masks>` of the configuration file.

.. note::

   Scale factors are assumed to be :literal:`unitless` (aka
   :literal:`1`) and no automatic unit conversion is performed.

.. _hco-cfg-masks:

=====
Masks
=====

This section lists all masks used by HEMCO. Masks are binary scale
factors (1 inside the mask region, 0 outside). If masks are regridded,
the remapped mask values (1 and 0) are determined through regular
rounding, i.e. a remapped mask value of 0.49 will be set to 0 while 0.5
will be set to 1.

The MASKS section in the HEMCO configuration file will look similar to
this (it will vary depending on the type of GEOS-Chem simulation you are
using):

.. code-block:: kconfig

   ###############################################################################
   ### BEGIN SECTION MASKS
   ###############################################################################
   # ScalID Name sourceFile sourceVar sourceTime CRE SrcDim SrcUnit Oper Lon1/Lat1/Lon2/Lat2

   #==============================================================================
   # Country/region masks
   #==============================================================================
   1000 EMEP_MASK   EMEP_mask.geos.1x1.20151222.nc           MASK     2000/1/1/0 C xy unitless 1 -30/30/45/70
   1002 CANADA_MASK Canada_mask.geos.1x1.nc                  MASK     2000/1/1/0 C xy unitless 1 -141/40/-52/85
   1003 SEASIA_MASK SE_Asia_mask.generic.1x1.nc              MASK     2000/1/1/0 C xy unitless 1  60/-12/153/55
   1004 NA_MASK     NA_mask.geos.1x1.nc                      MASK     2000/1/1/0 C xy unitless 1 -165/10/-40/90
   1005 USA_MASK    usa.mask.nei2005.geos.1x1.nc             MASK     2000/1/1/0 C xy unitless 1 -165/10/-40/90
   1006 ASIA_MASK   MIX_Asia_mask.generic.025x025.nc         MASK     2000/1/1/0 C xy unitless 1 46/-12/180/82
   1007 NEI11_MASK  USA_LANDMASK_NEI2011_0.1x0.1.20160921.nc LANDMASK 2000/1/1/0 C xy 1        1 -140/20/-50/60
   1008 USA_BOX     -129/25/-63/49                           -        2000/1/1/0 C xy 1        1 -129/25/-63/49

  ### END SECTION MASKS ###

The required attributes for mask fields are described below:

Options :option:`ScalID` and :option:`Oper` are described in
:ref:`hco-cfg-scalefac`.

Options :option:`Name`, :option:`sourceFile`, :option:`sourceVar`,
:option:`sourceTime`, :option:`CRE`, :option:`SrcDim`, and
:option:`SrcUnit`, are described in :ref:`hco-cfg-base`.

The :envvar:`Box` option is deprecated.

Instead of specifying the :option:`sourceFile` and :option:`sourceVar`
fields, you can directly provide the lower left and upper right box
coordinates: :literal:`Lon1/Lat1/Lon2/Lat2` . Longitudes must be in
degrees east, latitudes in degrees north. Only grid boxes whose mid
points are within the specified mask boundaries.  You may also specify
a single grid point (:literal:`Lon1/Lat1/Lon1/Lat1/`).

.. _hco-cfg-masks-caveat:

Caveat for simulations using cropped horizontal grids
-----------------------------------------------------

Consider the following combination of global and regional
emissions inventories:

In the :ref:`Base Emissions <hco-cfg-base>` section:

.. code-block:: kconfig

   0 GLOBAL_INV_SPC1   ...  SPC1 -     1 5
   0 INVENTORY_1_SPC1  ...  SPC1 1001  1 56
   0 INVENTORY_2_SPC1  ...  SPC1 1002  1 55

In the :ref:`Masks <hco-cfg-masks>` section:

.. code-block:: kconfig

   1001 REGION_1_MASK ... 1 1 70/10/140/60
   1002 REGION_2_MASK ... 1 1 46/-12/180/82

For clarity, we have omitted the various elements in these entries of
:file:`HEMCO_Config.rc` that are irrelevant to this issue.

With this setup, we should expect the following behavior:

   #. Species :literal:`SPC1` should be emitted globally from inventory
      :literal:`GLOBAL_INV` (hierarchy = 5). |br|
      |br|

   #. Regional emissions of :literal:`SPC1` from
      :literal:`INVENTORY_1` (hierarchy = 56) should overwrite global
      emissions in the region specified by :literal:`REGION_1_MASK`. |br|
      |br|

   #. Likewise, regional emissions of :literal:`SPC1` from
      :literal:`INVENTORY_2` (hierarchy = 55) should overwrite global
      emissions in the region specified by :literal:`REGION_2_MASK`. |br|
      |br|

   #. In the locations where :literal:`REGION_2_MASK` intersects
      :literal:`REGION_1_MASK`, emissions from
      :literal:`INVENTORY_1` will be applied. This is because
      :literal:`INVENTORY_1` has a higher hierarchy (56) than
      :literal:`INVENTORY_2` (55).

When running simulations that use cropped grids, one or both of the
boundaries specified for the masks (:literal:`70/10/140/60` and
:literal:`46/-12/180/82`) in :file:`HEMCO_Config.rc` can potentially
extend beyond the bounds of the simulation domain.  If this should
happen, HEMCO would treat the regional inventories as if they were
global, the emissions for the highest hierarchy (i.e.,
:literal:`INVENTORY_1`) would be applied globally. Inventories with
lower hierarchies would be ignored.

.. tip::

   Check the HEMCO log output for messages to make sure that none of
   your desired emissions have been skipped.

The solution is to make the boundaries of each defined mask region at least
a little bit smaller than the boundaries of the nested domain. This
involves inspecting the mask itself to make sure that no relevant
gridboxes will be excluded.

For example, assuming the simulation domain extends from 70E to 140E
in longitude, using this mask definition:

.. code-block:: kconfig

   1001 REGION_1_MASK ... 1 1 70/10/136/60

would prevent :literal:`INVENTORY_1` from being mistakely treated as a
global inventory. We hope to add improved error checking for this
condition into a future HEMCO version.

.. _hco-cfg-data-coll:

================
Data collections
================

The fields listed in :ref:`the HEMCO configuration file <hco-cfg>` data
collections. Collections can be enabled/disabled in section extension
switches. Only fields that are part of an enabled collection will be
used by HEMCO.

The beginning and end of a collection is indicated by an opening and
closing bracket, respectively: :literal:(((CollectionName` and
:literal:`)))CollectionName`. These brackets must be on individual lines
immediately preceeding / following the first/last entry of a collection.
The same collection bracket can be used as many times as needed.

The collections are enabled/disabled in the Extension Switches section
(see :ref:`Extension Switches <hco-cfg-ext-switches>`). Each
collection name must be provided as an extension setting and can then
be readily enabled/disabled:

.. code-block:: kconfig

   ###############################################################################
   #### BEGIN SECTION EXTENSION SWITCHES
   ###############################################################################
   # ExtNr ExtName           on/off  Species
   0       Base              : on    *
       --> MACCITY           :       true
       --> EMEP              :       true
       --> AEIC              :       true

   ### END SECTION EXTENSION SWITCHES

   ###############################################################################
   ### BEGIN SECTION BASE EMISSIONS
   ###############################################################################
   ExtNr Name srcFile srcVar srcTime CRE Dim Unit Species ScalIDs Cat Hier

   (((MACCITY
   0 MACCITY_CO MACCity.nc  CO 1980-2014/1-12/1/0 C xy  kg/m2/s CO 500      1 1
   )))MACCITY

   (((EMEP
   0 EMEP_CO    EMEP.nc     CO 2000-2014/1-12/1/0 C xy  kg/m2/s CO 500/1001 1 2
   )))EMEP

   (((AEIC
   0 AEIC_CO    AEIC.nc     CO 2005/1-12/1/0      C xyz kg/m2/s CO -        2 1
   )))AEIC

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

.. _hco-cfg-data-coll-ext-names:

Extension names
---------------

The collection brackets also work with :ref:`extension names
<hco-ext-list>`, e.g. data can be included/excluded based on
extensions. This is particularly useful to include an emission
inventory for standard emission calculation if (and only if) an
extension is not being used (see example below).

.. _hco-cfg-data-coll-undefined:

Undefined collections
---------------------

If, for a given collection, no corresponding entry is found in the
extensions section, it will be ignored. Collections are also ignored if
the collection is defined in an extension that is disabled. It is
recommended to list all collections under the base extension.

.. _hco-cfg-data-coll-exclude:

Exclude collections
-------------------

To use the opposite of a collection switch, .not. can be added in
front of an existing collection name. For instance, to read file
NOT_EMEP.nc only if EMEP is not being used:

.. code-block:: kconfig

   (((.not.EMEP
   0 NOT_EMEP_CO    $ROOT/NOT_EMEP.nc     CO 2000/1-12/1/0 C xy kg/m2/s CO 500/1001 1 2
   ))).not.EMEP

.. _hco-cfg-data-coll-combine:

Combine collections
-------------------

Multiple collections can be combined so that they are evaluated
together. This is achieved by linking collection names with .or..
For example, to use BOND biomass burning emissions only if both GFED and
FINN are not being used:

.. code-block:: kconfig

   (((.not.GFED.or.FINN
   0 BOND_BM_BCPI   $ROOT/BCOC_BOND/v2014-07/Bond_biomass.nc BC   2000/1-12/1/0 C xy kg/m2/s BCPI 70 2 1
   0 BOND_BM_BCPO   -                                        -    -             - -  -       BCPO 71 2 1
   0 BOND_BM_OCPI   $ROOT/BCOC_BOND/v2014-07/Bond_biomass.nc OC   2000/1-12/1/0 C xy kg/m2/s OCPI 72 2 1
   0 BOND_BM_OCPO   -                                        -    -             - -  -       OCPO 73 2 1
   0 BOND_BM_POA1   -                                        -    -             - -  -       POA1 74 2 1
   ))).not.GFED.or.FINN
