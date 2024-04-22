Requirements
============

.. include:: <isonum.txt>

.. role:: raw-html(raw)
    :format: html

Hardware requirements
---------------------

Computer system requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before you can run HEMCO, you will need to have one the following items.

+----------------------------------+----------------------------------+
| Item                             | Description                      |
+==================================+==================================+
| EITHER                           | You will need a Unix operating   |
|                                  | system environment in order to   |
| A Unix-based computer system     | run HEMCO. Any flavor of         |
|                                  | Unix (e.g. CentOS, Ubuntu,       |
|                                  | Fedora, etc.) should work just   |
|                                  | fine.                            |
+----------------------------------+----------------------------------+
| OR                               | If your institution has          |
|                                  | computational resources (e.g. a  |
| An account on the `Amazon Web    | shared computer cluster with     |
| Services                         | many cores, sufficient disk      |
| cloud <http://cloud-             | storage and memory), then you    |
| gc.readthedocs.io/en/stable/>`__ | can run GEOS-Chem there. Contact |
|                                  | your IT staff for assistance.    |
|                                  | :raw-html:`<br/>`                |
|                                  | :raw-html:`<br/>`                |
|                                  | If your institution lacks        |
|                                  | computational resources (or if   |
|                                  | you need additional              |
|                                  | computational resources beyond   |
|                                  | what is available), then you     |
|                                  | should consider signing up for   |
|                                  | access to the Amazon Web         |
|                                  | Services cloud. Using the cloud  |
|                                  | has the following advantages:    |
|                                  |                                  |
|                                  | -  You can run GEOS-Chem without |
|                                  |    having to invest in local     |
|                                  |    hardware and maintenance      |
|                                  |    personnel.                    |
|                                  | -  You won't have to download    |
|                                  |    any meteorological fields or  |
|                                  |    emissions data. All of the    |
|                                  |    necessary data input for      |
|                                  |    GEOS-Chem will be available   |
|                                  |    on the cloud.                 |
|                                  | -  You can initialize your       |
|                                  |    computational environment     |
|                                  |    with all of the required      |
|                                  |    software (e.g. compilers,     |
|                                  |    libraries, utilities) that    |
|                                  |    you need for GEOS-Chem.       |
|                                  | -  Your GEOS-Chem runs will be   |
|                                  |    100% reproducible, because    |
|                                  |    you will initialize your      |
|                                  |    computational environment the |
|                                  |    same way every time.          |
|                                  | -  You will avoid GEOS-Chem      |
|                                  |    compilation errors due to     |
|                                  |    library incompatibilities.    |
|                                  | -  You will be charged for the   |
|                                  |    computational time that you   |
|                                  |    use, and if you download data |
|                                  |    off the cloud.                |
|                                  |                                  |
|                                  | GEOS-Chem 12.0.0                 |
|                                  | and later versions can be used   |
|                                  | on the Amazon Web Services cloud |
|                                  | computing platform. You can      |
|                                  | learn more about how to use      |
|                                  | GEOS-Chem on the cloud by        |
|                                  | `visiting this tutorial          |
|                                  | (cloud.geos-chem.org)            |
|                                  | <http://cloud.geos-chem.org>`__. |
+----------------------------------+----------------------------------+


Memory requirements
~~~~~~~~~~~~~~~~~~~

If you plan to run GEOS-Chem on a local computer system, please make
sure that your system has sufficient memory and disk space:


Software requirements
---------------------


Supported compilers for HEMCO
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The table below lists the supported compilers for HEMCO.

HEMCO is written in the Fortran programming language. However, you
will also need C and C++ compilers to install certain libraries (like
netCDF) on your system.

+----------------+----------------+----------------+----------------+
| Item           | Description    | Versions       | Best way to    |
|                |                |                | install        |
+================+================+================+================+
| `Intel         | **The Intel    |                | `Install from  |
| Compiler Suite | Compiler Suite |                | Intel <http    |
| (icc, icpc,    | is our         | The GCST       | s://software.i |
| ifort)         | recommended    | has tested     | ntel.com/conte |
|                | proprietary    | with these     | nt/www/us/en/d |
|                | compiler       | versions (but  | evelop/tools/o |
|                | collection.**  | others may     | neapi/componen |
|                |                | work as well): | ts/fortran-com |
|                | Intel          |                | piler.html>`__ |
|                | compilers      | -  19.0.5.281  | (requires      |
|                | produce        | -  18.0.5      | purchase of a  |
|                | well-optimized | -  17.0.4      | site license   |
|                | code that runs | -  15.0.0      | or a student   |
|                | extremely      | -  13.0.079    | license)       |
|                | efficiency on  | -  11.1.069    |                |
|                | machines with  |                |                |
|                | Intel CPUs.    |                |                |
|                | Many           |                |                |
|                | universities   |                |                |
|                | and            |                |                |
|                | institutions   |                |                |
|                | will have an   |                |                |
|                | Intel site     |                |                |
|                | license that   |                |                |
|                | allows you to  |                |                |
|                | use these      |                |                |
|                | compilers.     |                |                |
+----------------+----------------+----------------+----------------+
| `GNU Compiler  | **The GNU      |                | `Install via   |
| Collection     | Compiler       |                | Spack <http:   |
| (gcc, g++,     | Collection is  | The GCST       | //github.com/s |
| gfortra        | our            | has tested     | pack/spack>`__ |
| n) <GNU_Fortra | recommended    | with these     |                |
| n_compiler>`__ | open-source    | versions (but  |                |
|                | compiler       | others may     |                |
|                | collection.**  | work as well): |                |
|                |                |                |                |
|                | Because the    | -  10.2.0      |                |
|                | GNU Compiler   | -  9.3.0       |                |
|                | Collection is  | -  9.2.0       |                |
|                | free and open  | -  8.2.0       |                |
|                | source, this   | -  7.4.0       |                |
|                | is a good      | -  7.3.0       |                |
|                | choice if your | -  7.1.0       |                |
|                | institution    | -  6.2.0       |                |
|                | lacks an Intel |                |                |
|                | site license,  |                |                |
|                | or if you are  |                |                |
|                | running        |                |                |
|                | GEOS-Chem on   |                |                |
|                | the Amazon EC2 |                |                |
|                | cloud          |                |                |
|                | environment.   |                |                |
+----------------+----------------+----------------+----------------+


Required software packages for HEMCO
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+----------------------+----------------------+----------------------+
| Item                 | Description          | Best way to install  |
+======================+======================+======================+
| `Git (a source code  | GEOS-Chem source     | Direct install from  |
| management           | code can be          | `Git-SCM.com         |
| system) <#T          | downloaded using the | <https://git-s       |
| he_Git_source_code_m | Git source code      | cm.com/downloads>`__ |
| anagement_system>`__ | management system.   |                      |
|                      | GEOS-Chem software   |                      |
|                      | repositories are     |                      |
|                      | stored at the        |                      |
|                      | https://             |                      |
|                      | github.com/geoschem  |                      |
|                      | organization page.   |                      |
|                      | Please see our       |                      |
|                      | `Guide to using Git  |                      |
|                      | with                 |                      |
|                      | GEOS-Chem            |                      |
|                      | <Guide_to_using_Gi   |                      |
|                      | t_with_GEOS-Chem>`__ |                      |
|                      | for more information |                      |
|                      | about how to use Git |                      |
|                      | with GEOS-Chem.      |                      |
+----------------------+----------------------+----------------------+
| `CMake <ht           | CMake is software    | `Install with        |
| tps://cmake.org/>`__ | that directs how the | Spack <http://github |
|                      | GEOS-Chem source     | .com/spack/spack>`__ |
|                      | code is compiled     |                      |
|                      | into an executable.  |                      |
|                      |                      |                      |
|                      | -  CMake is optional |                      |
|                      |    for GEOS-Chem     |                      |
|                      |    versions 12.6.0   |                      |
|                      |    through 12.9.3.   |                      |
|                      | -  CMake is          |                      |
|                      |    **REQUIRED** for  |                      |
|                      |    GEOS-Chem         |                      |
|                      |    versions 13.0.0   |                      |
|                      |    and later.        |                      |
+----------------------+----------------------+----------------------+
| `GNU                 | GNU Make is software | `Install with        |
| Make <ht             | that can build       | Spack <http://github |
| tps://cmake.org/>`__ | executables from     | .com/spack/spack>`__ |
|                      | source code.         |                      |
|                      |                      |                      |
|                      | -  NOTE: While GNU   |                      |
|                      |    Make is not       |                      |
|                      |    required for      |                      |
|                      |    GEOS-Chem 13.0.0  |                      |
|                      |    and later, some   |                      |
|                      |    external          |                      |
|                      |    libraries that    |                      |
|                      |    you might need to |                      |
|                      |    build will        |                      |
|                      |    require GNU Make. |                      |
|                      |    Therefore it is   |                      |
|                      |    best to download  |                      |
|                      |    GNU Make along    |                      |
|                      |    with CMake.       |                      |
+----------------------+----------------------+----------------------+
| `netCDF and          | GEOS-Chem input and  | `Install with        |
| netCDF-Fortran <#Th  | output data files    | Spack <http://github |
| e_netCDF_library>`__ | use the netCDF file  | .com/spack/spack>`__ |
|                      | format. This is a    |                      |
| -  plus dependencies | self-describing file |                      |
|    (e.g. HDF5, zlib, | format that allows   |                      |
|    etc)              | metadata             |                      |
|                      | (descriptive text)   |                      |
|                      | to be stored         |                      |
|                      | alongside data       |                      |
|                      | values. Please see   |                      |
|                      | our `Guide to netCDF |                      |
|                      | in                   |                      |
|                      | GEO                  |                      |
|                      | S-Chem <Guide_to_net |                      |
|                      | CDF_in_GEOS-Chem>`__ |                      |
|                      | for more information |                      |
|                      | about the netCDF     |                      |
|                      | file format and      |                      |
|                      | software library.    |                      |
+----------------------+----------------------+----------------------+


Optional but recommended software packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+----------------------+----------------------+----------------------+
| Item                 | Description          | Best way to install  |
+======================+======================+======================+
| `GCPy <https://gcp   | GCPy is our          | Install from         |
| y.readthedocs.io>`__ | recommended python   | conda-forge          |
|                      | companion software   |                      |
|                      | to GEOS-Chem. While  |                      |
|                      | this is not a        |                      |
|                      | general-purpose      |                      |
|                      | plotting package, it |                      |
|                      | does contain many    |                      |
|                      | useful functions for |                      |
|                      | creating zonal mean  |                      |
|                      | and horizontal plots |                      |
|                      | from GEOS-Chem       |                      |
|                      | output. It also      |                      |
|                      | contains scripts to  |                      |
|                      | generate plots and   |                      |
|                      | tables from          |                      |
|                      | GEOS-Chem benchmark  |                      |
|                      | simulations.         |                      |
+----------------------+----------------------+----------------------+
| `gdb                 | The GNU debugger     | `Install with        |
| <https://www.gnu.    | (gdb) and its        | Spack <http://github |
| org/software/gdb>`__ | graphical interface  | .com/spack/spack>`__ |
| and                  | (cgdb) are very      |                      |
| `cgdb <https:/       | useful tools for     |                      |
| /cgdb.github.io/>`__ | tracking down the    |                      |
|                      | source of GEOS-Chem  |                      |
|                      | errors, such as      |                      |
|                      | segmentation faults, |                      |
|                      | out-of-bounds        |                      |
|                      | errors, etc.         |                      |
+----------------------+----------------------+----------------------+
| `ncview              | ncview is a netCDF   | `Install with        |
| <http://meteora.uc   | file viewer. While   | Spack <http://github |
| sd.edu/~pierce/ncvie | it does not produce  | .com/spack/spack>`__ |
| w_home_page.html>`__ | publication-quality  |                      |
|                      | output, ncview can   |                      |
|                      | let you easily       |                      |
|                      | examine the contents |                      |
|                      | of a netCDF data     |                      |
|                      | file (such as those  |                      |
|                      | which are input and  |                      |
|                      | output by            |                      |
|                      | GEOS-Chem). Ncview   |                      |
|                      | is very useful for   |                      |
|                      | debugging and        |                      |
|                      | development.         |                      |
+----------------------+----------------------+----------------------+
| `nco                 | NCO are the netCDF   | `Install with        |
| <http://meteora.uc   | operators. These are | Spack <http://github |
| sd.edu/~pierce/ncvie | very powerful        | .com/spack/spack>`__ |
| w_home_page.html>`__ | command-line tools   |                      |
|                      | for editing and      |                      |
|                      | manipulating data in |                      |
|                      | netCDF format.       |                      |
+----------------------+----------------------+----------------------+
| `cdo <https          | CDO are the Climate  | `Install with        |
| ://code.mpimet.mpg.d | Data Operators.      | Spack <http://github |
| e/projects/cdo/l>`__ | These are very       | .com/spack/spack>`__ |
|                      | powerful             |                      |
|                      | command-line tools   |                      |
|                      | for editing and      |                      |
|                      | manipulating data in |                      |
|                      | netCDF format.       |                      |
+----------------------+----------------------+----------------------+
| `The Kinetic         | KPP translates a     | `Install with        |
| PreProcessor (KPP)   | chemical mechanism   | Git <FlexChem        |
| chemical             | specification from   | #KPP_source_code>`__ |
| sol                  | user-configurable    |                      |
| ver <https://github. | input files to       |                      |
| com/geoschem/kpp>`__ | Fortran-90 source    |                      |
|                      | code.                |                      |
+----------------------+----------------------+----------------------+
| `flex                | Flex is the Fast     | `Install with        |
| <https://github      | Lexical Analyzer.    | Spack <http://github |
| .com/westes/flex>`__ | This is a required   | .com/spack/spack>`__ |
|                      | library for the  KPP |                      |
|                      | chemical solver.     |                      |
+----------------------+----------------------+----------------------+


Required source code and data for HEMCO
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+----------------------+----------------------+----------------------+
| Item                 | Description          | Best way to install  |
+======================+======================+======================+
| `A clone of the      | The HEMCO            | Git clone from       |
| geoschem/HEMCO       | (Harmonized          | geoschem/HEMCO       |
| repository           | Emissions Component  |                      |
| repository <h        | codebase.            |                      |
| ttps://github.com/ge |                      |                      |
| oschem/HEMCO>`__     |                      |                      |
+----------------------+----------------------+----------------------+
| The GEOS-Chem shared | This is the          | `Perform a HEMCO     |
| data directories     | directory structure  | dry                  |
|                      | containing the       | run <Downloading_da  |
|                      | meteorology and      | ta_with_the_GEOS-Che |
|                      | emissions data that  | m_dry-run_option>`__ |
|                      | GEOS-Chem reads as   |                      |
|                      | input. For more      |                      |
|                      | information, please  |                      |
|                      | see our `Downloading |                      |
|                      | GEOS-Chem data       |                      |
|                      | directories <Do      |                      |
|                      | wnloading_GEOS-Chem_ |                      |
|                      | data_directories>`__ |                      |
|                      | wiki page.           |                      |
+----------------------+----------------------+----------------------+
