

Running HEMCO
=============

.. note::
   Another useful resource for instructions on running HEMCO is our `YouTube tutorial <https://www.youtube.com/watch?v=6Bup9V0ts6U&t=69s>`_.

Run interactively
-----------------

HEMCO may be run interactively at the command line by typing the following within your run directory

.. code-block:: console

   $ ./hemco_standalone

You may also specify the path to the HEMCO standalone configuration file using:

.. code-block:: console

   $ ./hemco_stantalone -c HEMCO_sa_Config.rc
   
If not specified, :file:`HEMCO_sa_Config.rc` will be used by default.

Run as batch job
----------------

Batch job run scripts will vary based on what job scheduler you have available. 
The example run script included in HEMCO run directories (:file:`runHEMCO.sh`) is for use
with SLURM. You may modify this file for your system and preferences as needed.

At the top of all batch job scripts are configurable run settings. 
Most critically are requested # cores, # nodes, time, and memory. 
Figuring out the optimal values for your run can take some trial and error. 

To submit a batch job using SLURM:

.. code-block:: console

   $ sbatch runHEMCO.sh
   
Standard output will be sent to a log file :file:`HEMCO_SA.log` once the job is started. 
Standard error will be sent to a file specific to your scheduler, e.g. :file:`slurm-jobid.out`
if using SLURM, unless you configure your run script to do otherwise.

If your computational cluster uses a different job scheduler, e.g. Grid Engine, LSF, or PBS, check with your IT staff or search the internet for how to configure and submit batch jobs. 
For each job scheduler, batch job configurable settings and acceptable formats are available on the internet and are often accessible from the command line. 
For example, type :command:`man sbatch` to scroll through options for SLURM, including various ways of specifying number of cores, time and memory requested.

Verify a successful run
-----------------------

There are several ways to verify that your run was successful.

1. NetCDF files are present in the :file:`OutputDir/` subdirectory
2. HEMCO log file :file:`HEMCO.log` ends with :literal:`HEMCO X.Y.Z FINISHED.`
3. Standard output file :file:`HEMCO_SA.log` ends with :literal:`HEMCO_STANDALONE FINISHED!`
4. The job scheduler log does not contain any error messages

If it looks like something went wrong, scan through the log files to determine where there may have been an error. Here are a few debugging tips:

* Review all of your configuration files to ensure you have proper setup
* Check to make sure you have downloaded all input files needed for your HEMCO standalone simulation

If you cannot figure out where the problem is please do not hesitate to create a GitHub issue.