#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile
#
# !DESCRIPTION: 
#\\
#\\
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#                                                                             .
#   make -jN TARGET REQUIRED-FLAGS [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#                                                                             .
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% NOTE: Normally you will not have to call this Makefile directly,     %%%
# %%% it will be called automatically from the Makefile in the directory   %%%
# %%% just above this one!                                                 %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                             .
# Makefile uses the following variables:
#                                                                             .
# Variable   Description
# --------   -----------
# SHELL      Specifies the shell for "make" to use (usually SHELL=/bin/sh)
# ROOTDIR    Specifies the root-level directory of the GEOS-Chem code
# HDR        Specifies the directory where GEOS-Chem include files are found
# LIB        Specifies the directory where library files (*.a) are stored
# MOD        Specifies the directory where module files (*.mod) are stored
# AR         Sys var w/ name of library creator program (i.e., "ar", "ranlib")
# MAKE       Sys var w/ name of Make command (i.e, "make" or "gmake")
# R8         Specifies the c
#
# !REVISION HISTORY: 
#  18 Sep 2013 - M. Long     - Initial version
#  17 Oct 2014 - R. Yantosca - Added "gigc_debug" target to print debug output
#  17 Oct 2014 - R. Yantosca - Added "the_nuclear_option" target to totally
#                              clean the ESMF, MAPL, and FVDYCORE directories
#  17 Oct 2014 - R. Yantosca - Cosmetic changes
#  18 Nov 2014 - R. Yantosca - Remove the fvdycore.install file again after
#                              making distclean, for safety's sake.
#  12 Dec 2014 - R. Yantosca - Define defaults for  ESMF, MAPL, FVdycore
#                              that can be overridden by cmd line or env var
#  21 Mar 2018 - H.P. Lin    - Remove FV3Dycore / Cubed-Sphere & "geos" EXE Deps.
#                              To restore a 1-D version of GIGC.
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Use the bash shell
SHELL=/bin/bash

# ROOTDIR is the directory just above this one.
ROOTDIR=..

# Directories where GEOS-Chem code live
HDR=$(ROOTDIR)/Headers
HELP=$(ROOTDIR)/help
LIB=$(ROOTDIR)/lib
MOD=$(ROOTDIR)/mod

# This code here should get the exact number of the intel version (MDY)
ifeq ($(ESMF_COMPILER),intel)
  INTELVERSIONTEXT :=$(shell ifort --version))
  INTELVERSIONTEXT :=$(sort $(INTELVERSIONTEXT))
  LONGVERSION      :=$(word 3, $(INTELVERSIONTEXT))
  MAJORVERSION     :=$(subst ., ,$(LONGVERSION))
  MAJORVERSION     :=$(firstword $(MAJORVERSION))
  MAJORVERSION     :=$(strip $(MAJORVERSION))
else ifeq ($(ESMF_COMPILER),gfortran)
  GNUVERSIONTEXT   :=$(shell gfortran -dumpversion))
  LONGVERSION      :=$(GNUVERSIONTEXT))
  MAJORVERSION     :=$(subst ., ,$(LONGVERSION))
  MAJORVERSION     :=$(firstword $(MAJORVERSION))
  MAJORVERSION     :=$(strip $(MAJORVERSION))
else
  MAJORVERSION     :="0"
endif

# Include header file.  This returns variables CC, F90, FREEFORM, LD, R8,
# as well as the default Makefile compilation rules for source code files.
include $(ROOTDIR)/Makefile_header.mk

# BASEDIR is a synonym for ROOTDIR
ifndef BASEDIR
 export BASEDIR=$(realpath $(ROOTDIR))
endif 

###############################################################################
###                                                                         ###
###  List of files to compile                                               ###
###                                                                         ###
###############################################################################

# List of source files
SRC = $(wildcard *.F) $(wildcard *.F90)

# Replace .F and .F90 extensions with *.o
TMP = $(SRC:.F=.o)
OBJ = $(TMP:.F90=.o)

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: clean help baselibs baselibs_mapl

baselibs:
	@echo "GCHP Shared libraries do not need to be compiled if no ESMF"

baselibs_esmf:
	@echo "No ESMF Installation"

baselibs_mapl:
	@echo "No MAPL Installation"

lib: $(OBJ)
	$(AR) crs libGIGC.a $(OBJ)
	mv libGIGC.a $(LIB)

libesmf:
	@echo "No ESMF Installation"

libmapl:
	@echo "No MAPL Installation"

clean:
	rm -f *.o *.mod *___.h *___.rc

help:
	@$(MAKE) -C $(HELP)

gigc_debug gigc_help:
	@echo "Directories:"
	@echo "----------------------------------------------------------"
	@echo "Current working dir    : $(CURDIR)"
	@echo "ROOTDIR                : $(ROOTDIR)"
	@echo "HDR                    : $(HDR)"
	@echo "HELP                   : $(HELP)"
	@echo "LIB                    : $(LIB)"
	@echo "MOD                    : $(MOD)"
	@echo "BASEDIR                : $(BASEDIR)"

###############################################################################
###                                                                         ###
###  Targets to remove ESMF, MAPL, even EVERYTHING GEOS-CHEM!               ###
###  USE WITH EXTREME CAUTION!!!                                            ###
###                                                                         ###
###############################################################################

.PHONY: the_nuclear_option the_thermonuclear_option

the_nuclear_option:
	@$(MAKE) clean

the_thermonuclear_option:
	cd $ROOTDIR; @$(MAKE) HPC=y realclean
	@$(MAKE) clean


###############################################################################
###                                                                         ###
###  Dependencies listing                                                   ###
###  (grep "USE " to get the list of module references!)                    ###
###                                                                         ###
###  From this list of dependencies, the "make" utility will figure out     ###
###  correct order of compilation (so we don't have to do that ourselves).  ###
###  This also allows us to compile on multiple processors with "make -j".  ###
###                                                                         ###
###  NOTES:                                                                 ###
###  (1) Only specify object-file dependencies that are within this         ###
###       directory.  Object files in other directories will be referenced  ### 
###       at link-time.                                                     ###
###  (2) For "make -jN" (i.e. compile N files simultaneously), all files    ###
###       in this directory must have a listed dependency.                  ###
###                                                                         ###
###############################################################################

gchp_utils.o                : gchp_utils.F90

gigc_stateful_mod.o         : gigc_stateful_mod.F90

gigc_chunk_mod.o            : gigc_chunk_mod.F90 gchp_utils.o gigc_stateful_mod.o

