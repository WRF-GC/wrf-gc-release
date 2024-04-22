.. _hco-sa-compile:

####################
Build the executable
####################

.. note::

   Another useful resource for HEMCO build instructions is our
   `YouTube tutorial
   <https://www.youtube.com/watch?v=6Bup9V0ts6U&t=69s>`_.

Once you have created a :ref:`run directory <hco-sa-rundir>`, you may
proceed to compile the HEMCO standalone source code into an executable
file.  You will compile HEMCO standalone from your :ref:`run directory
<hco-sa-rundir>`.

There are two steps to build HEMCO. The first step is to **configure your
build settings** with :ref:`hco-sa-soft-cmake`.  Build settings cover
options like enabling or disabling components or whether HEMCO should
be compiled in :option:`Debug` mode.

The second step is to **compile the source code into an executable**.
For this, you will use :ref:`make <hco-sa-soft-make>`, which builds the
executable according to your build settings.

.. _hco-sa-compile-builddir:

================================
Navigate to your build directory
================================

A subdirectory named :file:`build` is included in each :ref:`HEMCO
standalone run directory <hco-sa-rundir>` that you create. You can use
this directory (known as a **build directory**) to create the HEMCO
executable file.

You are not limited to using the build directory that is created
inside the :ref:`run directory <hco-sa-rundir>`.  In fact, you can
create as many build directories you wish in whatever location you
wish.  For example, if you want to compare HEMCO standalone
performance on both :ref:`GNU <hco-sa-soft-gnu>` and :ref:`Intel
<hco-sa-soft-intel>` compilers, you could create two different build
directories, one named :file:`build_gnu` and the other
:file:`build_intel`. Build directories do not necessarily need
to be kept in the :ref:`run directory <hco-sa-rundir>`, but it is
convenient to do so.

Each build directory is self-contained, so you can delete it at any
point to erase the HEMCO standalone build and its configuration.  Most
users will typically only need to build HEMCO standalone once, so we
recommend using the :file:`build` subdirectory of the :ref:`run
directory <hco-sa-rundir>` as the location to create the HEMCO
standalone exectuable.

.. important::

   There is one rule for build directories: **a build directory should
   be a new directory**.

In the example below, we will use the :file:`build` directory within
the :ref:`run directory <hco-sa-rundir>` to build the HEMCO standalone
executable.

Navigate to the run directory:

.. code-block:: console

   $ cd /path/to/hemco/run/dir

Then navigate to the :file:`build` folder within:

.. code-block:: console

   $ cd build

.. _ hco-sa-compile-init:

==============================
Initialize the build directory
==============================

Run :ref:`hco-sa-soft-cmake` to initialze the build directory.

.. code-block:: console

   $ cmake ../CodeDir -DRUNDIR=..

:file:`CodeDir` is a symbolic link to the HEMCO source code
directory.

The option :literal:`-DRUNDIR=..` specifies that the
directory where we will run HEMCO standalone is one level above
us. This makes sense as our :file:`build` folder is a subdirectory of
the run directory.  (More about :ref:`build options
<hco-sa-compile-options>` below:

You will see output similar to this:

.. code-block:: console

   -- The Fortran compiler identification is GNU 11.2.0
   -- Detecting Fortran compiler ABI info
   -- Detecting Fortran compiler ABI info - done
   -- Check for working Fortran compiler: /path/to/gfortran - skipped
   -- Checking whether /path/to/gfortran supports Fortran 90
   -- Checking whether /path/to/gfortran supports Fortran 90 - yes
   =================================================================
   HEMCO X.Y.Z
   Current status: X.Y.Z
   =================================================================
   -- Found OpenMP_Fortran: -fopenmp (found version "4.5")
   -- Found OpenMP: TRUE (found version "4.5")
   -- Found NetCDF: /path/to/netcdf/lib/libnetcdff.so
   -- Bootstrapping  /path/to/hemco/run/dir
   -- Settings:
     * OMP:          ON  OFF
     * USE_REAL8:    ON  OFF
   -- Configuring done
   -- Generating done
   -- Build files have been written to: /path/to/hemco/run/dir/build

In the above example output, the version number :literal:`X.Y.Z` will
refer to the actual HEMCO version number (e.g. :literal:`3.4.0`,
:literal:`3.5.0`, etc.).  Also the paths :file:`/path/to/...` in your
output instead be the actual paths to the compiler and libraries.

.. _hco-sa-compile-config:

======================
Configuring your build
======================

Build settings are controlled by :ref:`hco-sa-soft-cmake` commands
with the following form:

.. code-block:: console

   $ cmake . -D<NAME>=<VALUE>

where :literal:`<NAME>` is the name of the setting, and
:literal:`<VALUE>` is the value that you are assigning it. These
settings are persistent and saved in your build
directory. You can set multiple variables in a single command, and you
can run :ref:`hco-sa-soft-cmake` as many times as you need to
configure your desired settings.

.. note::

   The :literal:`.` argument is important. It is the path to your
   build directory which is :literal:`.` here.

HEMCO has no required build settings. You can find the complete list
of :ref:`HEMCO's build settings here <hco-sa-compile-options>`. The most
frequently used build setting is :literal:`RUNDIR` which lets you
specify one or more run directories where CMake will install
HEMCO. Here, "install" refers to copying the compiled executable, and
some supplemental files with build settings, to your run directories.

Since there are no required build settings, for this tutorial we will
stick with the default settings.

You should notice that when you run :ref:`hco-sa-soft-cmake` it ends with:

.. code-block:: console

   ...
   -- Configuring done
   -- Generating done
   -- Build files have been written to: /path/to/hemco/run/dir/build

This tells you the configuration was successful, and that you are
ready to compile.

.. _hco-sa-compile-compile:

========================
Compile HEMCO standalone
========================

Compile HEMCO standalone with this command

.. code-block:: console

   $ make -j

The :command:`-j` option will tell :ref:`hco-sa-soft-make` to compile
several source code files in parallel.  This reduces overall
compilation time.

Optionally, you can use the :literal:`VERBOSE=1` argument to see the
compiler commands.

This step creates :file:`./bin/hemco_standalone` which is the compiled
executable. You can copy this executable to your run directory
manually, or you can do

.. code-block:: console

   $ make install

which copies :file:`./bin/hemco_standalone` (and some supplemental
files) to the run directories specified in :option:`RUNDIR`.

Now you have compiled HEMCO!  You can now navigate back from the
:file:`build` folder to the run directory (which we remember is one
level higher):

.. code-block:: console

   $ cd ..

.. _hco-sa-compile-recompile:

=========================================
Recompile when you change the source code
=========================================

You need to recompile :program:`HEMCO` if you update a build setting
or make a modification to the source code. However, with
:ref:`hco-sa-soft-cmake`, you don't need to clean before
recompiling. The build system automatically figure out which
files need to be recompiled based on your modification. This is known
as incremental compiling.

To recompile HEMCO standalone, simply do

.. code-block:: console

   $ make -j
   $ make install

which will recompile HEMCO standalone and copy the new executable file
to the run directory.

.. _hco-sa-compile-options:

==============================
HEMCO standalone build options
==============================

.. option:: RUNDIR

   Paths to run directories where :command:`make install` installs
   HEMCO standalone. Multiple run directories can be specified by a
   semicolon separated list. A warning is issues if one of these
   directories does not look like a run directory.

   These paths can be relative paths or absolute paths. Relative paths
   are interpreted as relative to your build directory.

.. option:: CMAKE_BUILD_TYPE

   Specifies the build type.  Allowable values are:

   .. option:: Release

      The default option.  Compiles HEMCO standalone for speed.

   .. option:: Debug

      Compiles HEMCO standalone with several debugging flags turned
      on.  This may help you find common errors such as
      array-out-of-bounds, division-by-zero, or not-a-number.

      .. important::

         The additional error checks that are applied with
         :literal:`Debug` will cause HEMCO standalone to run much more
         slowly!  Do not use :literal:`Debug` for long production
         simulations.

.. option:: HEMCO_Fortran_FLAGS_<COMPILER_ID>

    Additional compiler options for HEMCO standalone for build type
    :literal:`<BUILD_TYPE>`.

   .. option:: <COMPILER_ID>

      Valid values are :literal:`GNU` and :literal:`Intel`.

.. option:: HEMCO_Fortran_FLAGS_<CMAKE_BUILD_TYPE>_<COMPILER_ID>

   Compiler options for HEMCO standalone for the given
   :option:`CMAKE_BUILD_TYPE`.

   .. option:: <COMPILER_ID>

      Valid values are :literal:`GNU` and :literal:`Intel`.
