#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (Main-level)
#
# !DESCRIPTION: This is a "router" makefile.  It calls the main GEOS-Chem 
# Makefile (in the GeosCore subdirectory) to direct the Unix "make" utility 
# how to build the GEOS-Chem source code.
#\\
#\\
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#                                                                             .
#   make -jN TARGET REQUIRED-FLAGS [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#                                                                             .
# Makefile uses the following variables:
#                                                                             .
# Variable   Description
# --------   -----------
# GEOSDIR    Specifies the directory where GEOS-Chem "core" routines are found
# GTMM       Specifies the directory where the GTMM routines are found
#
# !REVISION HISTORY: 
#  16 Sep 2009 - R. Yantosca - Initial version
#  24 Nov 2009 - R. Yantosca - Now call libbpch and libcore targets in
#                              the Makefile in the GeosCore sub-directory
#  11 Dec 2009 - R. Yantosca - Now get SHELL from Makefile_header.mk
#  25 Jan 2010 - R. Yantosca - Added Makefile targets for TOMAS microphysics
#  16 Feb 2011 - R. Yantosca - Added Makefile targets for APM microphysics
#  04 Nov 2011 - R. Yantosca - Remove ESMF targets, those are not needed
#  24 Jan 2012 - R. Yantosca - Also add libnc target to build netCDF utils
#  11 May 2012 - R. Yantosca - Now make sure that all targets of the 
#                              GeosCore/Makefile are pointed to properly
#  20 Aug 2013 - R. Yantosca - Make sure Makefile names are consistent
#  18 Sep 2013 - R. Yantosca - Remove GeosTomas, that is now gone
#  18 Sep 2013 - M. Long     - Add Makefile target "hpc"
#  15 Jan 2014 - R. Yantosca - Updated comments
#  19 Mar 2014 - R. Yantosca - Add more visible comment section dividers
#  04 Jun 2015 - R. Yantosca - Now add "wipeout" and "debug" targets
#  23 Jun 2016 - R. Yantosca - Remove references to GeosApm
#EOP
#------------------------------------------------------------------------------
#BOC

# Directories
GEOSDIR :=GeosCore
GTMM    :=GTMM

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete list!                ###
###                                                                         ###
###############################################################################

.PHONY: all lib libcore libheaders libkpp libiso libnc librad libutil
.PHONY: exe clean realclean doc docclean tauclean help wipeout debug
.PHONY: realclean_except_rrtmg libhemcosa

all:
	@$(MAKE) -C $(GEOSDIR) all

hpc:
	@$(MAKE) -C $(GEOSDIR) hpc

lib:
	@$(MAKE) -C $(GEOSDIR) lib

libcore:
	@$(MAKE) -C $(GEOSDIR) libcore

libheaders:
	@$(MAKE) -C $(GEOSDIR) libheaders

libiso:
	@$(MAKE) -C $(GEOSDIR) libiso

libkpp:
	@$(MAKE) -C $(GEOSDIR) libkpp

libnc:
	@$(MAKE) -C $(GEOSDIR) libnc	

ncdfcheck:
	@$(MAKE) -C $(GEOSDIR) ncdfcheck

librad:
	@$(MAKE) -C $(GEOSDIR) librad

libhemcosa:
	@$(MAKE) -C $(GEOSDIR) libhemcosa

libutil:
	@$(MAKE) -C $(GEOSDIR) libutil

exe:
	@$(MAKE) -C $(GEOSDIR) exe

clean:
	@$(MAKE) -C $(GEOSDIR) clean

distclean:
	@$(MAKE) -C $(GEOSDIR) distclean

realclean:
	@$(MAKE) -C $(GEOSDIR) realclean

realclean_except_rrtmg:
	@$(MAKE) -C $(GEOSDIR) realclean_except_rrtmg

doc:
	@$(MAKE) -C $(GEOSDIR) doc

docclean: 
	@$(MAKE) -C $(GEOSDIR) docclean

tauclean:
	find . -name '*.pdb' -o -name '*.inst.*' -o -name '*.pp.*' -o -name '*.continue.*' -o -name '*.tau.inc~' | xargs rm -f

debug:
	@$(MAKE) -C $(GEOSDIR) debug

wipeout:
	@$(MAKE) -C $(GEOSDIR) wipeout

help:
	@$(MAKE) -C $(GEOSDIR) help

headerinfo:
	@$(MAKE) -C $(GEOSDIR) headerinfo

###############################################################################
###                                                                         ###
###  Targets for Hg simulation w/ Global Terrestrial Mercury Model (GTMM)   ###
###                                                                         ###
###############################################################################

.PHONY: allhg libhg libgtmm exehg

allhg:
	@$(MAKE) -C $(GEOSDIR) GTMM_Hg=yes allhg

libhg:
	@$(MAKE) -C $(GEOSDIR) GTMM_Hg=yes libhg

ligbtmm:
	@$(MAKE) -C $(GEOSDIR) GTMM_Hg=yes libgtmm

exehg:
	@$(MAKE) -C $(GEOSDIR) GTMM_Hg=yes exehg

#EOC



