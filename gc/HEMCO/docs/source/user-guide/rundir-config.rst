

Configuring a run directory
===========================

.. note::
   Another useful resource for instructions on configuring HEMCO run directories is our `YouTube tutorial <https://www.youtube.com/watch?v=6Bup9V0ts6U&t=69s>`_.

Navigate to your new run directory, and examine the contents:

.. code-block:: console

   $ cd /scratch/rundirs/hemco_4x5_merra2
   $ ls
   build/           HEMCO_Diagn.rc        HEMCO_sa_Spec.rc  README
   CodeDir@         HEMCO_sa_Config.rc    HEMCO_sa_Time.rc  rundir.version
   HEMCO_Config.rc  HEMCO_sa_Grid.4x5.rc  OutputDir/        runHEMCO.sh*
   
The following files can be modified to set up your HEMCO standalone simulation.
   
HEMCO_sa_Config.rc
   Main configuration file for the HEMCO standalone simulation. This file points to the other 
   configuration files used to set up your simulation (e.g. :file:`HEMCO_sa_Grid.4x5`, 
   :file:`HEMCO_sa_Time.rc`). This file typically references a HEMCO_Config.rc file using
   :literal:`>>>include HEMCO_Config.rc` which contains the emissions settings. Settings in
   :file:`HEMCO_sa_Config.rc` will always override any settings in the included
   :file:`HEMCO_Config.rc`.

HEMCO_Config.rc
   Contains emissions settings. This file is typically obtained
   from another model (e.g. GEOS-Chem).

HEMCO_Diagn.rc
   Specifies which fields to save out to the HEMCO diagnostics file saved in
   :file:`OutputDir` by default. The frequency to save out diagnostics is controlled 
   by the :literal:`DiagnFreq` setting in :file:`HEMCO_Config_sa.rc`
   
HEMCO_sa_Grid.4x5.rc
   Defines the grid specification. Sample files are provided for 4.0 x 5.0, 2.0 x 2.5,
   0.5 x 0.625, and 0.25 x 0.3125 global grids in :file:`HEMCO/run/` and are automatically copied
   to the run directory based on options chosen when running :file:`createRunDir.sh`. If
   you choose to run with a custom grid or over a regional domain, you will need to modify
   this file manually.

HEMCO_sa_Spec.rc
   Defines the species to include in the HEMCO standalone simulation. By default, the
   species in a GEOS-Chem full-chemistry simulation are defined. To include other species,
   you can modify this file by providing the species name, molecular weight, and other
   properties.

HEMCO_sa_Time.rc
   Defines the start and end times of the HEMCO standalone simulation as well as the emissions
   timestep (s).

runHEMCO.sh
   Sample run script for submitting a HEMCO standalone simulation via SLURM.