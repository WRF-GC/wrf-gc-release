
Compiling HEMCO
===============

.. note::
    This user guide assumes you have loaded a computing environment that satisfies
    :ref:`HEMCO's software requirements <software_requirements>`.

.. note::
   Another useful resource for HEMCO build instructions is our `YouTube tutorial <https://www.youtube.com/watch?v=6Bup9V0ts6U&t=69s>`_.

There are two steps to build HEMCO. The first step is configuring your build. 
To configure your build you use :program:`cmake` to configure build settings. 
Build settings cover options like enabling or disabling components, specifying run directories
to install HEMCO to, or whether HEMCO should be compiled in Debug mode. 

The second step is compiling. To compile HEMCO you use :program:`make`. This
compiles HEMCO according to your build configuration.

Create your build directory
---------------------------

Create a build directory. This directory is going to be the working directory
for your build. The configuration and compile steps generate a 
bunch of build files, and this directory is going to store those. You can
think of a build directory as representing a HEMCO build. It stores configuration
settings, information about your system, and intermediate files from the compiler.

A build directory is self contained, so you can delete it at any point to erase 
the build and its configuration. You can have as many build directories as you 
would like. Most users only need one build directory, since they only build HEMCO
once; but, for example, if you were building HEMCO with Intel and GNU compilers to
compare performance, you would have two build directories: one for the Intel build,
and one for the GNU build. You can name your build directories whatever you want, but a good choice is :file:`build/`.
There is one rule for build directories: **a build directory should be a new directory**.

Create a build directory and initialize it. You initialize a build directory by
running :program:`cmake` with the path to the HEMCO source code. Here is an example
of creating a build directory in the top-level of the HEMCO source code:

.. code-block:: console

   $ cd HEMCO
   $ mkdir build
   $ cd build
   $ cmake ..
   -- The Fortran compiler identification is Intel 18.0.5.20180823
   -- Check for working Fortran compiler: /bin/intel64/ifort
   -- Check for working Fortran compiler: /bin/intel64/ifort  -- works
   ...
   -- Configuring done
   -- Generating done
   -- Build files have been written to: HEMCO/build


Configuring your build
----------------------

Build settings are controlled by :program:`cmake` commands with the following
form:

.. code-block:: none

    $ cmake . -D<NAME>="<VALUE>"

where :literal:`<NAME>` is the name of the setting, and :literal:`<VALUE>` is the
value that you are assigning it. These settings are persistent and saved in your build directory.
You can set multiple variables in a single command, and you can run :program:`cmake` as many times
as you need to configure your desired settings.

.. note:: 
   The :literal:`.` argument is important. It is the path to your build directory which
   is :literal:`.` here.

HEMCO has no required build settings. You can find the complete list of :ref:`HEMCO's build settings here <hemco_build_options>`.
The most frequently used build setting is :literal:`RUNDIR` which lets you specify one or more run directories
where CMake will install HEMCO. Here, "install" refers to copying the compiled executable, and some supplemental files
with build settings, to your run directories.

.. note::
    You can even update build settings after you compile HEMCO. Simply rerun :program:`make` and
    (optionally) :program:`make install`, and the build system will automatically figure out
    what needs to be recompiled.

Since there are no required build settings, for this tutorial we will stick with the
default settings. 

You should notice that when you run :program:`cmake` it ends with:

.. code-block:: console
   
   ...
   -- Configuring done
   -- Generating done
   -- Build files have been written to: HEMCO/build

This tells you the configuration was successful, and that you are ready to compile. 

Compile HEMCO
-------------

You compile HEMCO with:

.. code-block:: console

   $ make -j    # -j enables compiling in parallel

Optionally, you can use the :literal:`VERBOSE=1` argument to see the compiler commands.

This step creates :file:`./bin/hemco_standalone` which is the compiled executable. You can copy
this executable to your run directory manually, or you can do

.. code-block:: console
   
   $ make install

which copies :file:`./bin/hemco_standalone` (and some supplemental files) to 
the run directories specified in :ref:`RUNDIR <build_setting_rundir>`.

Now you have compiled HEMCO, and you are ready to move on to creating a run directory!

------------

Recompiling
-----------

You need to recompile HEMCO if you update a build setting or make a modification to the source code.
However, with CMake, you don't need to clean before recompiling. The build system automatically 
figure out which files need to be recompiled based on your modification. This is known as incremental compiling.

To recompile HEMCO, simply do 

.. code-block:: console
   
   $ make -j   # -j enables compiling in parallel

and optionally, do :command:`make install`.

------------

.. _hemco_build_options:

RUNDIR
   Paths to run directories where :command:`make install` installs HEMCO. Multiple
   run directories can be specified by a semicolon separated list. A warning is 
   issues if one of these directories does not look like a run directory.

   These paths can be relative paths or absolute paths. Relative paths are interpreted as relative to your build directory.

CMAKE_BUILD_TYPE
    The build type. Valid values are :literal:`Release` or :literal:`Debug`.
    Set this to :literal:`Debug` if you want to build in debug mode.

HEMCO_Fortran_FLAGS_<COMPILER_ID>
    Additional compiler options for HEMCO for build type :literal:`<BUILD_TYPE>`.
    
HEMCO_Fortran_FLAGS_<BUILD_TYPE>_<COMPILER_ID>
    Compiler options for HEMCO for all build types. Valid values for :literal:`<COMPILER_ID>` are :literal:`GNU` and
    :literal:`Intel`.
