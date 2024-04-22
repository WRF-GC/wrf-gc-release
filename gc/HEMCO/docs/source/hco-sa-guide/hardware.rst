.. _hco-sa-hard:

############################
Obtain the required hardware
############################

In this chapter, we provide information about the computer equipment
that you will need in order to run :program:`HEMCO` in standalone
mode (aka the :program:`HEMCO standalone`).

.. _hco-sa-hard-computer:

============================
Computer system requirements
============================

Before you can run HEMCO standalone, you will need to have
one the following items.

#. A Unix/Linux based computer system, OR:
#. An account on the `Amazon Web Services cloud computing platform
   <http://geos-chem-cloud.readthedocs.io/>`_.

If your institution has computational resources (e.g. a shared
computer cluster with many cores, sufficient disk storage and memory),
then you can run HEMCO standalone there.  Contact your IT
staff for assistance.

If your institution lacks computational resources (or if you need
additional computational resources beyond what is available), then you
should consider signing up for access to the Amazon Web Services
cloud. Using the cloud has the following advantages:

#. You can run HEMCO standalone without having to invest in
   local hardware and maintenance personnel.
#. You won't have to download any meteorological fields or emissions
   data. All of the necessary data input for HEMCO standalone
   will be available on the cloud.
#. You can initialize your computational environment with all of the
   required software (e.g. compilers,libraries, utilities) that you
   need for HEMCO standalone.
#. Your runs will be 100% reproducible, because you will initialize
   your computational environment the same way every time.
#. You will avoid compilation errors due to library incompatibilities.
#. You will be charged for the computational time that you use, and if
   you download data off the cloud.

.. _hco-sa-hard-mem-disk:

============================
Memory and disk requirements
============================

If you plan to run HEMCO standalone on a local computer
system, please make sure that your system has sufficient memory and
disk space.

We would recommend at least 4 GB of RAM to run HEMCO standalone.
However, if you will be reading data sets at very fine horizional
resolution, you will want to increase the memory to perhaps 20-30
GB/RAM.

Also make sure that you have enough disk space to store the amount of
input data for your HEMCO standalone simulations.
