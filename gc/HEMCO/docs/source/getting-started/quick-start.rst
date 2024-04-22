

Quick Start
===========

This quickstart guide assumes your environment satisfies :ref:`HEMCO's requirements <software_requirements>`. 
This means you should load a compute environment such that programs like :program:`cmake` and :program:`mpirun`
are available, before continuing. You can find more detailed instructions in the user guide.

1. Clone HEMCO
--------------

Download the source code:

.. code-block:: console

   $ git clone https://github.com/geoschem/HEMCO.git ~/HEMCO
   $ cd ~/HEMCO

Checkout the HEMCO version that you want to use:

.. code-block:: console

   $ git checkout 3.0.0-rc.0

2. Create a run directory
-------------------------

Navigate to the :file:`run/` subdirectory. 
Create a run directory by running :file:`./createRunDir.sh` and answering the prompts:

.. code-block:: console

   $ cd run/
   $ ./createRunDir.sh

3. Configure your build
-----------------------

Create a build directory and :command:`cd` into it. 
A good name for this directory is :file:`build/`, and a good place for it is in the 
top-level of the source code:

.. code-block:: console

   $ mkdir ~/HEMCO/build
   $ cd ~/HEMCO/build

Initialize your build directory by running :program:`cmake` and passing it the path to your source code:

.. code-block:: console

   $ cmake ~/HEMCO

Now you can configure :ref:`build options <HEMCO_build_options>`. 
These are persistent settings that are saved to your build directory.
A common build option is :literal:`-DRUNDIR`. 
This option lets you specify one or more run directories that HEMCO is "installed" to when you do :command:`make install`. 
Configure your build so it installs HEMCO to the run directory you created in Step 2:

.. code-block:: console

   $ cmake . -DRUNDIR="/path/to/rundir"

.. note::
   The :literal:`.` in the :program:`cmake` command above is important. It tells CMake that your 
   current working directory (i.e., :literal:`.`) is your build directory.

4. Compile and install
----------------------

Compile HEMCO:

.. code-block:: console

   $ make -j

Next, install the compiled executable to your run directory (or directories):

.. code-block:: console

   $ make install

This copies :file:`build/bin/hemco_standalone` and supplemental files to your run directory. 

.. note::
   You can update build settings at any time:
   
   1. Navigate to your build directory.
   2. Update your build settings with :program:`cmake`. See 
   3. Recompile with :command:`make -j`. Note that the build system automatically figures out what (if any) files
      need to be recompiled.
   4. Install the rebuilt executable with :command:`make install`.


5. Configure your run directory
-------------------------------

Now, navigate to your run directory:

.. code-block:: console

   $ cd path/to/rundir

Simulation settings are configured in the :file:`.rc` files. The main configuration file
is :file:`HEMCO_sa_Config.rc`. The start end end time for your simulation can be modified in
:file:`HEMCO_sa_Time.rc`. The horizontal grid for your simulation can be modified in
:file:`HEMCO_sa_Grid.rc`. Emissions settings can be changed in the `HEMCO_Config.rc` file
that has been copied from another model (e.g. GEOS-Chem).

6. Run HEMCO
------------

HEMCO can be run interactively from within your run directory by typing:

.. code-block:: console

   $ ./hemco_standalone

You may also submit your HEMCO simulation as a batch job to a scheduler.  A sample run script
:file:`runHEMCO.sh` is included in your run directory. To submit a HEMCO simulation using
SLURM:

.. code-block:: console

   $ sbatch runHEMCO.sh

Those are the basics of using HEMCO! See the user guide, step-by-step guides, and reference pages
for more detailed instructions.