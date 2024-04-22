.. _hco-sa-sim-config:

######################
Configure a simulation
######################

.. note::

   Another useful resource for instructions on configuring HEMCO run
   directories is our `YouTube tutorial
   <https://www.youtube.com/watch?v=6Bup9V0ts6U&t=69s>`_.

Navigate to your new directory, and examine the contents:

.. code-block:: console

   $ cd /path/to/hemco/run/dir
   $ ls
   build/           HEMCO_Diagn.rc        HEMCO_sa_Spec.rc  README
   CodeDir@         HEMCO_sa_Config.rc    HEMCO_sa_Time.rc  rundir.version
   HEMCO_Config.rc  HEMCO_sa_Grid.4x5.rc  OutputDir/        runHEMCO.sh*

The following files can be modified to set up your HEMCO standalone simulation.

.. option:: HEMCO_sa_Config.rc

   Main configuration file for the HEMCO standalone simulation. This
   file points to the other configuration files used to set up your
   simulation (e.g. :option:`HEMCO_sa_Grid.4x5.rc`,
   :option:`HEMCO_sa_Time.rc`).

   This file typically references a :option:`HEMCO_Config.rc` file
   using

   .. code-block:: none

      >>>include HEMCO_Config.rc

   which contains the emissions settings. Settings in
   :option:`HEMCO_sa_Config.rc` will always override any settings in
   the included :option:`HEMCO_Config.rc` file.

.. option:: HEMCO_Config.rc

   Contains emissions settings. :option:`HEMCO_Config.rc` can be taken
   from a another model (such as GEOS-Chem), or can be built from a
   sample file.

   For more information on editing :option:`HEMCO_Config.rc`, please
   see the following chapters: :ref:`hco-cfg`, :ref:`edit-hco-cfg`,
   and :ref:`cfg-ex`.

   .. important::

      Make sure that the path to your data directory in the
      :option:`HEMCO_Config.rc` file is correct.  Otherwise, HEMCO
      standalone will not be able read data from disk.

.. option:: HEMCO_Diagn.rc

   Specifies which fields to save out to the HEMCO diagnostics file
   saved in :file:`OutputDir` by default. The frequency to save out
   diagnostics is controlled by the :option:`DiagnFreq` setting in
   :option:`HEMCO_sa_Config.rc`

   For more information, please see the chapter entitled
   :ref:`hco-diag-configfile`.

.. option:: HEMCO_sa_Grid.4x5.rc

   Defines the grid specification. Sample files are provided for 4.0 x
   5.0, 2.0 x 2.5, 0.5 x 0.625, and 0.25 x 0.3125 global grids in
   :file:`HEMCO/run/` and are automatically copied to the run
   directory based on options chosen when running
   :file:`createRunDir.sh`.  you choose to run with a custom grid or
   over a regional domain, you will need to modify this file
   manually.

.. option:: HEMCO_sa_Spec.rc

   Defines the species to include in the HEMCO standalone
   simulation. By default, the species in a GEOS-Chem full-chemistry
   simulation are defined. To include other species, you can modify
   this file by providing the species name, molecular weight, and
   other properties.

.. option:: HEMCO_sa_Time.rc

   Defines the start and end times of the HEMCO standalone simulation
   as well as the emissions timestep (s).

.. option:: runHEMCO.sh

   Sample run script for submitting a HEMCO standalone simulation via
   SLURM.
