.. _hco-sa-login:

################################
Configure your login environment
################################

In this chapter, you will learn how to load the software packages that
you have created into your computational environment.  This will need
to be done each time you log in to your computer system.

.. tip::

   You may skip this section if you plan on using HEMCO standalone on
   an Amazon EC2 cloud instance.  When you initialize the EC2 instance
   with one of the pre-configured Amazon Machine Images (AMIs) all of
   the required software libraries will be automatically loaded.

An environment file does the following:

  1. Loads software libraries into your login environment.  This is
     often done with a module manager such as :command:`lmod`,
     :command:`spack`, or  :command:`environment-modules`.

  2. Stores settings for HEMCO and its dependent libraries in
     shell variables called `environment variables
     <https://www.networkworld.com/article/3215965/all-you-need-to-know-about-unix-environment-variables.html>`_.

Environment files allow you to easily switch between different sets of
libraries.  For example, you can keep one environment file to load the
Intel Compilers for HEMCO standalone and another to load
the GNU Compilers.

For general information about how libraries are loaded, see our
:ref:`Library Guide <libguide>` in the Supplemental Guides section.

We recommend that you place module load commands into a separate
**environment file**  rather than directly into your :file:`~/.bashrc`
or :file:`~/.bash_aliases` startup scripts.

.. _hco-sa-login-gnu:

================================================
Sample environment file for GNU 10.2.0 compilers
================================================

Below is a sample environment file from the Harvard Cannon computer
cluster.  This file will load software libraries built with the GNU
10.2.0 compilers.

Save the code below (with any appropriate modifications for your own
computer system) to a file named :file:`~/gnu102.env`.

.. code-block:: bash

   # Echo message if we are in a interactive (terminal) session
   if [[ $- = *i* ]] ; then
     echo "Loading modules for GEOS-Chem, please wait ..."
   fi

   #==============================================================================
   # Modules (specific to Cannon @ Harvard)
   #==============================================================================

   # Remove previously-loaded modules
   module purge

   # Load modules for GNU Compilers v10.2.0
   module load git/2.17.0-fasrc01
   module load gcc/10.2.0-fasrc01
   module load openmpi/4.1.0-fasrc01
   module load netcdf-fortran/4.5.3-fasrc03
   module load flex/2.6.4-fasrc01
   module load cmake/3.17.3-fasrc01

   #==============================================================================
   # Environment variables
   #==============================================================================

   # Parallelization settings
   export OMP_NUM_THREADS=8
   export OMP_STACKSIZE=500m

   # Make all files world-readable by default
   umask 022

   # Specify compilers
   export CC=gcc
   export CXX=g++
   export FC=gfortran

   # Netcdf variables for CMake
   # NETCDF_HOME and NETCDF_FORTRAN_HOME are automatically
   # defined by the "module load" commands on Cannon.
   export NETCDF_C_ROOT=${NETCDF_HOME}
   export NETCDF_FORTRAN_ROOT=${NETCDF_FORTRAN_HOME}

   # Set memory limits to max allowable
   ulimit -c unlimited              # coredumpsize
   ulimit -l unlimited              # memorylocked
   ulimit -u 50000                  # maxproc
   ulimit -v unlimited              # vmemoryuse
   ulimit -s unlimited              # stacksize

   # List modules loaded
   module list

.. tip::

   Ask your sysadmin how to load software libraries.  If you are using
   your institution's computer cluster, then chances are there will
   be a software module system installed, with commands similar to
   those listed above.

Then you can activate these seetings from the command line by typing:

.. code-block:: console

   $ source ~/gnu102.env

.. _hco-sa-login-intel:

==============================================
Sample environment file for Intel 19 compilers
==============================================

To load software libraries based on the Intel 19 compilers, we can
start from our :ref:`GNU 10.2.0 environment file <hco-sa-login-gnu>`
and add the proper :command:`module load` commands for Intel 19.

Add the code below (with the appropriate modifications for your
system) into a file named :file:`~/intel19.env`.

.. code-block:: bash

   # Echo message if we are in a interactive (terminal) session
   if [[ $- = *i* ]] ; then
     echo "Loading modules for GEOS-Chem, please wait ..."
   fi

   #==============================================================================
   # Modules (specific to Cannon @ Harvard)
   #==============================================================================

   # Remove previously-loaded modules
   module purge

   # Load modules for Intel compilers v19.0.4
   module load git/2.17.0-fasrc01
   module load intel/19.0.5-fasrc01
   module load openmpi/4.0.1-fasrc01
   module load netcdf-fortran/4.5.2-fasrc03
   module load flex/2.6.4-fasrc01
   module load cmake/3.17.3-fasrc01

   #==============================================================================
   # Environment variables
   #==============================================================================

   # Parallelization settings
   export OMP_NUM_THREADS=8
   export OMP_STACKSIZE=500m

   # Make all files world-readable by default
   umask 022

   # Specify compilers
   export CC=icc
   export CXX=icpc
   export FC=ifort

   # Netcdf variables for CMake
   # NETCDF_HOME and NETCDF_FORTRAN_HOME are automatically
   # defined by the "module load" commands on Cannon.
   export NETCDF_C_ROOT=${NETCDF_HOME}
   export NETCDF_FORTRAN_ROOT=${NETCDF_FORTRAN_HOME}

   # Set memory limits to max allowable
   ulimit -c unlimited              # coredumpsize
   ulimit -l unlimited              # memorylocked
   ulimit -u 50000                  # maxproc
   ulimit -v unlimited              # vmemoryuse
   ulimit -s unlimited              # stacksize

   # List modules loaded
   module list

.. tip::

   Ask your sysadmin how to load software libraries.  If you
   are using your institution's computer cluster, then chances
   are there will be a software module system installed, with
   commands similar to those listed above.

Then you can activate these seetings from the command line by typing:

.. code-block:: console

   $ source intel19.env

.. tip::

   Keep a separate environment file for each combination of
   modules that you will load.

.. _hco-sa-envvar-compilers:

=======================================
Set environment variables for compilers
=======================================

Add the following environment variables to your environment file to
specify the compilers that you wish to use:

.. table:: Environment variables that specify the choice of compiler
   :align: center

   +---------------+------------------+--------------------+-----------------+
   | Variable      | Specifies the:   | GNU name           | Intel name      |
   +===============+==================+====================+=================+
   | :envvar:`CC`  | C compiler       | :envvar:`gcc`      | :envvar:`icc`   |
   +---------------+------------------+--------------------+-----------------+
   | :envvar:`CXX` | C++ compiler     | :envvar:`g++`      | :envvar:`icpc`  |
   +---------------+------------------+--------------------+-----------------+
   | :envvar:`FC`  | Fortran compiler | :envvar:`gfortran` | :envvar:`ifort` |
   +---------------+------------------+--------------------+-----------------+

These environment variables should be defined in your
:ref:`environment file <hco-sa-login>`.

.. note::

   Only the Fortran compiler is needed to compile the HEMCO
   standalone.  But if you need to :ref:`manually install libraries
   <build-libraries-with-spack>`, you will also need the C and C++
   compilers.

.. _hco-sa-envvar-parallel:

=============================================
Set environment variables for parallelization
=============================================

The HEMCO standalone` uses `OpenMP parallelization
<Parallelizing_GEOS-Chem>`_, which is an implementation of
shared-memory (aka serial) parallelization.

.. important::

   OpenMP-parallelized programs cannot execute on more than 1
   computational node.  Most modern computational nodes typically
   contain  between 16 and 64 cores. Therefore, HEMCO standalone
   simulations will not be able to take advantage of more cores than
   these.

Add the following environment variables to your environment file to
control the OpenMP parallelization settings:

.. option:: OMP_NUM_THREADS

   The :envvar:`OMP_NUM_THREADS` environment variable sets the number of
   computational cores (aka threads) to use.

   For example, the command below will tell HEMCO standalone to use 8
   cores within parallel sections of code:

   .. code:: console

      $ export OMP_NUM_THREADS=8

.. option:: OMP_STACKSIZE

   In order to use HEMCO standalone with `OpenMP 
   parallelization <Parallelizing_GEOS-Chem>`_, you must request the
   maximum amount of stack memory in your login environment. (The
   stack memory is where local automatic variables and temporary
   :envvar:`!$OMP PRIVATE` variables will be created.) Add the
   following lines to your system startup file and to your GEOS-Chem
   run scripts:

   .. code-block:: bash

      ulimit -s unlimited
      export OMP_STACKSIZE=500m

   The :command:`ulimit -s unlimited` will tell the bash shell to use the
   maximum amount of stack memory that is available.

   The environment variable :envvar:`OMP_STACKSIZE` must also be set to a very
   large number. In this example, we are nominally requesting 500 MB of
   memory. But in practice, this will tell the GNU Fortran compiler to use
   the maximum amount of stack memory available on your system. The value
   **500m** is a good round number that is larger than the amount of stack
   memory on most computer clusters, but you can increase this if you wish.

.. _errors_caused_by_incorrect_settings:

=======================================
Fix errors caused by incorrect settings
=======================================

Be on the lookout for these errors:

  #. If :option:`OMP_NUM_THREADS` is set to 1, then your
     HEMCO standalone simulation will execute using only
     one computational core.  This will make your simulation take much
     longer than is necessary.

  #. If :option:`OMP_STACKSIZE` environment variable is not included
     in your environment file (or if it is set to a very low value),
     you might encounter a **segmentation fault**.  In this case,
     the HEMCO standalone "thinks" that it does not have
     enough memory to perform the simulation, even though sufficient
     memory may be present.
