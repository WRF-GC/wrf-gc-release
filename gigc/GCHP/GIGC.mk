#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: GIGC.mk
#
# !DESCRIPTION: Makefile fragment that specifies include and library
#  paths for the various phases of the GIGC build sequence.
#\\
#\\
# !REMARKS:
#  Updated from Mike Long
#
# !REVISION HISTORY:
#  18 Nov 2014 - M. Long     - Initial version
#  18 Nov 2014 - R. Yantosca - Now use env vars to specify MPI inc & lib dirs
#  01 Dec 2014 - R. Yantosca - Now put FV_LIB before MPI_LIB in link command
#EOP
#------------------------------------------------------------------------------
#BOC

#==============================================================================
# (1) Root directories
#==============================================================================

# %%%%% Root dir for MAPL etc %%%%%
ifndef ESMADIR
export ESMADIR=$(PWD)/GCHP/Shared
endif

# %%%%% Root dir for ESMF %%%%%
ifndef ESMF_DIR
export ESMF_DIR=$(PWD)/GCHP/ESMF
endif

# %%%%% Root dir for FVdycore %%%%%
# This directory is a dummy as we do not use fv3 in wrf-gc
ifndef FVDIR
 export FVDIR=$(PWD)/GCHP/FVdycoreCubed_GridComp
endif

#==============================================================================
# (2) MPI settings
#
# The following are setting for various versions of MPI. Users will need to
# uncomment the settings for MPI_LIB based on the appropriate version. 
# Hopefully this can be automated in future versions.
#==============================================================================
ifeq ($(ESMF_COMM),openmpi)
   # %%%%% OpenMPI settings %%%%%
   MPI_LIB       := $(shell mpif90 --showme:link)
   MPI_LIB       += $(shell mpicxx --showme:link)
else ifeq ($(ESMF_COMM),mvapich2)
   # %%%%% MVAPICH %%%%% 
   MPI_LIB       := -L$(dir $(shell which mpif90))../lib64 -lmpich -lmpichf90
else ifeq ($(ESMF_COMM),mpi)
   # %%%%% Generic MPI (works for SGI) %%%%%
   MPI_LIB       := -L$(dir $(shell which mpif90))../lib -lmpi -lmpi++
else
   $(error ESMF_COMM not defined or not valid at GIGC.mk)
endif

MPI_INC       := $(dir $(shell which mpif90))../include

#==============================================================================
# (3) GIGC/GEOS-Chem general settings
#
# The following are environment settings for GIGC to compile within the
# GEOS-Chem framework. They are dependent upon the settings in Section (1).
#==============================================================================

# %%%%% Architecture %%%%%
ifndef ARCH
  ARCH := $(shell uname -s)
endif
# %%%%% ESMF settings %%%%%
ESMF_MOD      := 
ESMF_INC      := 
ESMF_LIB      := 

# %%%%% MAPL settings %%%%%
MAPL_INC      := 
MAPL_LIB      := 

# %%%%% Link command %%%%%
LINK          := -lGIGC $(MPI_LIB) $(LINK) -lGIGC

# %%%%% Fortran flags %%%%%
FFLAGS        := -double-size 32 -real-size 32 -r4
USER_FFLAGS   += -DMODEL_ -DMODEL_WRF -DSPMD -DGLOBAL_GRID_CREATE
USER_DEFS     += -DMODEL_ -DMODEL_WRF

# %%%%% SDE 2013-03-26: Let HEMCO standalone see MAPL
# Note this is a fix for hemco_standalone.x, NOT for the final "geos" EXE binary.
# Don't get confused and remove this - I sure did (2018-03-21 hplin)
LINK_HCO      += $(MPI_LIB)