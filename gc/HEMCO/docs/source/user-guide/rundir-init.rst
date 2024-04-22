
Creating a Run Directory
========================

.. note::
   Another useful resource for HEMCO run directory creation instructions is our `YouTube tutorial <https://www.youtube.com/watch?v=6Bup9V0ts6U&t=69s>`_.

HEMCO run directories are created from within the source code.
A new run directory should be created for each different version of HEMCO you use. 
Git version information is logged to file :file:`rundir.version` within the run directory upon creation.

To create a run directory, navigate to the :file:`run/` subdirectory of the source code and execute shell script :file:`createRunDir.sh`.

.. code-block:: console

   $ cd HEMCO/run
   $ ./createRunDir.sh

During the course of script execution you will be asked a series of questions:

Enter ExtData path
------------------

The first time you create a HEMCO run directory on your system you will be prompted for a path to GEOS-Chem shared data directories,
which are also used by HEMCO.
The path should include the name of your :file:`ExtData/` directory and should not contain symbolic links. 
The path you enter will be stored in file :file:`~/.geoschem/config` in your home directory as environment variable :envvar:`GC_DATA_ROOT`. 
If that file does not already exist it will be created for you. 
When creating additional run directories you will only be prompted again if the file is missing or if the path within it is not valid.

.. code-block:: none

   -----------------------------------------------------------
   Enter path for ExtData:
   -----------------------------------------------------------

Choose meteorology source
-------------------------

Enter the integer number that is next to the input meteorology source you would like to use.

.. code-block:: none

   -----------------------------------------------------------
   Choose meteorology source:
   -----------------------------------------------------------
     1. MERRA-2 (Recommended)
     2. GEOS-FP

Choose horizontal resolution
----------------------------

Enter the integer number that is next to the horizontal resolution you would like to use.

.. code-block:: none

   -----------------------------------------------------------
   Choose horizontal resolution:
   -----------------------------------------------------------
     1. 4.0 x 5.0
     2. 2.0 x 2.5
     3. 0.5 x 0.625
     4. 0.25 x 0.3125
     5. Custom


Enter HEMCO_Config.rc path
--------------------------

Provide the path to a HEMCO_Config.rc file with your emissions settings. This is typically
obtained from another model (e.g. :file:`~/GEOS-Chem/run/HEMCO_Config.rc.templates/HEMCO_Config.rc.fullchem`)

.. code-block:: none

   -----------------------------------------------------------
   Enter path to the HEMCO_Config.rc file with your emissions settings.
   
   NOTE: This may be a HEMCO_Config.rc file from a GEOS-Chem run directory
   or a HEMCO_Config.template file from the GEOS-Chem source code repository.
   -----------------------------------------------------------

Enter run directory path
------------------------

Enter the target path where the run directory will be stored. You will be prompted to enter a new path if the one you enter does not exist.

.. code-block:: none

   -----------------------------------------------------------
   Enter path where the run directory will be created:
   -----------------------------------------------------------

Enter run directory name
------------------------

Enter the run directory name, or accept the default. You will be prompted for a new name if a run directory of the same name already exists at the target path.

.. code-block:: none

   -----------------------------------------------------------
   Enter run directory name, or press return to use default:
   
   NOTE: This will be a subfolder of the path you entered above.
   -----------------------------------------------------------

Enable version control (optional)
---------------------------------

Enter whether you would like your run directory tracked with git version control. 
With version control you can keep track of exactly what you changed relative to the original settings. 
This is useful for trouble-shooting as well as tracking run directory feature changes you wish to migrate back to the standard model.

.. code-block:: none

   -----------------------------------------------------------
   Do you want to track run directory changes with git? (y/n)
   -----------------------------------------------------------

If a run directory has successfully been created, you should see something like:

.. code-block:: none

   Created /scratch/rundirs/hemco_4x5_merra2