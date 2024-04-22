.. _hco-sa-run:

################
Run a simulation
################

.. note::

   Another useful resource for instructions on running :program:`HEMCO` is our
   `YouTube tutorial
   <https://www.youtube.com/watch?v=6Bup9V0ts6U&t=69s>`_.

.. _hco-sa-run-int:

=================
Run interactively
=================

First, navigate to your run directory (if you aren't already there):

.. code-block:: console

   $ cd /path/to/hemco/run/dir

You can run HEMCO standalone interactively at the command line by typing:

.. code-block:: console

   $ ./hemco_standalone -c HEMCO_sa_Config.rc

where :literal:`-c` specifies the path to the
:option:`HEMCO_sa_Config.rc` configuraiton file.

.. _hco-sa-run-batch:

================
Run as batch job
================

Batch job run scripts will vary based on what job scheduler you have
available. The example run script included in HEMCO standalone run
directories (:file:`runHEMCO.sh`) is for use with SLURM. You may
modify this file for your system and preferences as needed.

At the top of all batch job scripts are configurable run
settings. Most critically are requested # cores, # nodes, time, and
memory. Figuring out the optimal values for your run can take some
trial and error.

To submit a batch job using SLURM:

.. code-block:: console

   $ sbatch runHEMCO.sh

Standard output will be sent to a log file :file:`HEMCO_SA.log` once
the job is started. Standard error will be sent to a file specific to
your scheduler, e.g. :file:`slurm-jobid.out` if using SLURM, unless
you configure your run script to do otherwise.

If your computational cluster uses a different job scheduler,
e.g. Grid Engine, LSF, or PBS, check with your IT staff or search the
internet for how to configure and submit batch jobs. For each job
scheduler, batch job configurable settings and acceptable formats are
available on the internet and are often accessible from the command
line. For example, type :command:`man sbatch` to scroll through
options for SLURM, including various ways of specifying number of
cores, time and memory requested.

.. _hco-sa-run-verify:

=======================
Verify a successful run
=======================

There are several ways to verify that your run was successful.

- :ref:`NetCDF <ncguide>` files are present in the :file:`OutputDir/`
  subdirectory;
- The HEMCO log file :file:`HEMCO.log` ends with :literal:`HEMCO X.Y.Z
  FINISHED.`;
- Standard output file :file:`HEMCO_SA.log` ends with
  :literal:`HEMCO_STANDALONE FINISHED!`;
- The job scheduler log does not contain any error messages

If it looks like something went wrong, scan through the log files to
determine where there may have been an error. Here are a few debugging
tips:

- Review all of your configuration files to ensure you have proper setup
- Check to make sure you have downloaded all input files needed for
  your HEMCO standalone simulation.

If you cannot figure out where the problem is please do not hesitate
to create a `GitHub issue
<https://github.com/geoschem/HEMCO/issues/new/choose/>`_.
