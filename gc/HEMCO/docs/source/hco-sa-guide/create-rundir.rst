.. _hco-sa-rundir:

######################
Create a run directory
######################

.. note::
   Another useful resource for :program:`HEMCO standalone` run
   directory creation instructions is our `YouTube tutorial
   <https://www.youtube.com/watch?v=6Bup9V0ts6U&t=69s>`_.

HEMCO standalone run directories are created from within the source code.
A new run directory should be created for each different version of
HEMCO you use. Git version information is logged to file
:file:`rundir.version` within the run directory upon creation.

To create a run directory, navigate to the :file:`run/` subdirectory
of the source code and execute shell script :file:`createRunDir.sh`.

.. code-block:: console

   $ cd HEMCO/run
   $ ./createRunDir.sh

During the course of script execution you will be asked a series of
questions:

==================
Enter ExtData path
==================

The first time you create a HEMCO standalone run directory on your
system you will be prompted for a path to the :file:`ExtData` folder,
which is the root data directory for HEMCO (as well as for `GEOS-Chem
<https://geos-chem.readthedocs.io>`_).

The path that you specify  should include the name of your
:file:`ExtData/` directory and should not contain symbolic links.  The
path you enter will be stored in file :file:`~/.geoschem/config` in
your home directory as environment variable :envvar:`GC_DATA_ROOT`. If
that file does not already exist it will be created for you. When
creating additional run directories you will only be prompted again if
the file is missing or if the path within it is not valid.

.. code-block:: console

   -----------------------------------------------------------
   Enter path for ExtData:
   -----------------------------------------------------------

=========================
Choose meteorology source
=========================

Enter the integer number that is next to the input meteorology source
you would like to use.

.. code-block:: console

   ===========================================================
   HEMCO STANDALONE RUN DIRECTORY CREATION
   ===========================================================

   -----------------------------------------------------------
   Choose meteorology source:
   -----------------------------------------------------------
     1. MERRA-2 (Recommended)
     2. GEOS-FP
     3. GISS ModelE2.1 (GCAP 2.0)

============================
Choose horizontal resolution
============================

Enter the integer number that is next to the horizontal resolution you
would like to use.

.. code-block:: console

   -----------------------------------------------------------
   Choose horizontal resolution:
   -----------------------------------------------------------
     1. 4.0 x 5.0
     2. 2.0 x 2.5
     3. 0.5 x 0.625
     4. 0.25 x 0.3125
     5. Custom

==========================
Enter HEMCO_Config.rc path
==========================

Provide the path to a :file:`HEMCO_Config.rc` file with your emissions
settings.

.. code-block:: console

   -----------------------------------------------------------
   Enter the file path to a HEMCO_Config.rc with your
   emissions settings.

    - This should be a HEMCO_Config.rc file from a
      pre-generated GEOS-Chem run directory and not a
      template config file from the GEOS-Chem repository.

    - If you do not have a pre-generated HEMCO_Config.rc file,
      type ./HEMCO_Config.rc.sample at the prompt below.
      This will copy a sample configuration file into your
      run directory.  You can then edit this configuration
      file with your preferred emission settings.
   -----------------------------------------------------------

If you have a pre-configured :file:`HEMCO_Config.rc` file available
(e.g. from a `GEOS_Chem <https://geos-chem.readthedocs.io>`_ run
directory), then then type the absolute path:

.. code-block:: console

   /path/to/my/HEMCO_Config.rc

If you do not have a :file:`HEMCO_Config.rc` template file handy, then
type:

.. code-block:: console

   ./HEMCO_Config.rc.sample

This will copy sample :file:`HEMCO_Config.rc` and
:file:`HEMCO_Diagn.rc` files to the run directory.  You can edit these
configuration files to include your preferred emission settings.

Refer to the :ref:`HEMCO Reference Guide <hco-ref-guide>`
for more information about how to edit :ref:`the HEMCO configuration
file <hco-cfg>`.

========================
Enter run directory path
========================

Enter the target path where the run directory will be stored. You will
be prompted to enter a new path if the one you enter does not exist.

.. code-block:: console

   -----------------------------------------------------------
   Enter path where the run directory will be created:
   -----------------------------------------------------------

========================
Enter run directory name
========================

Enter the run directory name, or accept the default. You will be
prompted for a new name if a run directory of the same name already
exists at the target path.

.. code-block:: console

   -----------------------------------------------------------
   Enter run directory name, or press return to use default:

   NOTE: This will be a subfolder of the path you entered above.
   -----------------------------------------------------------

If you press return, a default name such as :file:`hemco_4x5_merra2`,
:file:`hemco_2x25_geosfp`, etc. will be used.

=================================
Enable version control (optional)
=================================

Enter whether you would like your run directory tracked with
:ref:`hco-sa-soft-git` version control.  With version control you can
keep track of exactly what you changed relative to the original
settings. This is useful for trouble-shooting as well as tracking run
directory feature changes you wish to migrate back to a previous
version.

.. code-block:: console

   -----------------------------------------------------------
   Do you want to track run directory changes with git? (y/n)
   -----------------------------------------------------------

If a run directory has successfully been created, the name of the run
directory will be printed.  If you used the default run directory name
then you will see output similar to:

.. code-block:: console

   Created /path/to/hemco_4x5_merra2

etc.

======================
Run directory contents
======================

Navigate to the run directory that was just created and get a
directory listing:

.. code-block:: console

   $ cd hemco_4x5_merra2
   $ ls
   build/    HEMCO_Config.rc  HEMCO_sa_Config.rc    HEMCO_sa_Spec.rc  OutputDir/  rundir.version
   CodeDir@  HEMCO_Diagn.rc   HEMCO_sa_Grid.4x5.rc  HEMCO_sa_Time.rc  README      runHEMCO.sh*

:file:`build` is the folder is where you will :ref:`compile HEMCO
standalone <hco-sa-compile>`.

:file:`CodeDir` is a symbolic link back to the HEMCO source code.

:file:`OutputDir` is the folder where diagnostic outputs will be generated.

Files ending in :file:`.rc` are user-edtiable configuration files
that control HEMCO standalone simulation options.  We will discuss
these in more detail more in the :ref:`hco-sa-sim-config` chapter.

The :file:`rundir.version` file contains information about the Git
commit in the HEMCO source code corresponding to this run directory.
You will see output similar to this:

.. code-block:: console

   This run directory was created with /path/to/hemco/HEMCO/run/createRunDir.sh.

   HEMCO repository version information:

     Remote URL: git@github.com:geoschem/hemco.git
     Branch: dev
     Commit: Add fixes for generating HEMCO standalone run directory
     Date: Wed Jul 13 10:56:36 2022 -0400
     User: Melissa Sulprizio
     Hash: b29dac4

   Changes to the following run directory files are tracked by git:

    [master (root-commit) b8e694d] Initial run directory
    7 files changed, 477 insertions(+)
    create mode 100644 HEMCO_Config.rc
    create mode 100644 HEMCO_Diagn.rc
    create mode 100644 HEMCO_sa_Config.rc
    create mode 100644 HEMCO_sa_Grid.4x5.rc
    create mode 100644 HEMCO_sa_Spec.rc
