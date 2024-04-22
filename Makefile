################################################################################
#
#   WRF-GCHP
#   GEOS-Chem High Performance-powered Chemistry Add-On for WRF Model
#
#   WRF & GCHP are (c) their original authors.
#   WRF-GCHP coupling layer (WGCL) is (c) Atmospheric Chemistry and Climate Group, Peking University
#
#   Developed by Haipeng Lin <linhaipeng@pku.edu.cn>, Xu Feng, 2018-01
#   Peking University, School of Physics
#
################################################################################
#
#   Codename Pumpkin: Abstracted Bindings for Chemistry-to-WRF
#
#   This Chemical Interface (chem/) is written after comprehensive study of
#   the original chem_driver.f from WRF-Chem v3.6.1
#   which is (c) their respective authors.
#
################################################################################
#
#  Makefile
#
#  Compiles WRF-GCHP-Pumpkin (Chemistry-to-WRF bindings)
#  and the chemistry component (GEOS-Chem), bindings based on GEOS-Chem HP
#
#  Author: Haipeng Lin <linhaipeng@pku.edu.cn>, April 2nd, 2018
#
################################################################################
SHELL := /bin/bash

LN      =       ln -sf
MAKE    =       make -i -r
RM      =       rm -f
MODULES =                                 \
        module_data_radm2.o \
        module_data_sorgam.o \
        module_data_mosaic_asect.o \
        module_chem_utilities.o \
        module_convection_prep.o \
        module_interpolate.o \
        module_tropopause.o \
        module_upper_bc_driver.o \
        module_mosaic_driver.o \
        module_aerosols_sorgam.o \
        module_gocart_aerosols.o \
        module_aerosols_soa_vbs.o \
        module_input_chem_data.o \
        module_input_tracer.o \
        module_input_chem_bioemiss.o \
        module_optical_averaging.o \
        module_mixactivate_wrappers.o \
        module_diag_aero_size_info.o \
        module_data_rrtmgaeropt.o \
        module_data_mosaic_other.o \
        module_data_mosaic_therm.o \
        module_data_gigc_asect.o \
        module_peg_util.o \
        wrfgc_convert_state_mod.o

ifdef PNETCDF
MODULES +=	wrfgc_history_mod.o\
			wrfgc_io_pnetcdf.o
endif

OBJS    =                           \
        chemics_init.o              \
        optical_driver.o            \
        mixactivate_driver.o        \
        chem_driver.o

TARGETDIR    =  ./

ifndef ROOTDIR
  export ROOTDIR=$(realpath ..)
endif

# For some reason this gets lost...
export WRF_SRC_ROOT_DIR=$(ROOTDIR)

# Makefile Targets (.PHONY)
# Phony targets as defined by make will ALWAYS run when you ask to run them.
# We have added chemics to this list of targets - this is to support "repeat" compiles of WRF,
# as libwrflib.a will need to regen'd by the ./compile script sometimes, but chemics does not
# always run in this instance.
#
# If your make does not support phony targets... upgrade or check man make for FORCE
.PHONY: chemics clean devclean install_common install_registry install_registry_ch4 install_registry_co2 compile_chem troubleshooting test_registry_dummy

# This is the main target (chemics)
# It compiles modules, drivers, then links to create a libwrflib.a
# with command "ar ru" (Replace, u - newer/modified only)
#
# hplin 3/31/2018: Added "make compile_chem" to default chemics target.
# This compiles GEOS-Chem so we can link it later in MODULE/DRIVERS
#
# Also, "ar ru" seems to be problematic on some systems (such as mine)
# If this is the case, replace $(ARFLAGS) with "r". This may be related to an issue with the
# D modifier.
chemics: compile_chem MODULE DRIVERS
	$(AR) $(ARFLAGS) $(ROOTDIR)/main/$(LIBWRFLIB) $(MODULES) $(OBJS)

MODULE: $(MODULES)

DRIVERS: $(OBJS)

# Include WRF Environmental Settings...
# THIS FILE, CONTRARY TO COMMON SENSE, CONTAINS MAKEFILE TARGETS. THIS MEANS THAT IT MUST GO AFTER CHEMICS.
# IF YOU PUT THIS BEFORE CHEMICS, THE WORLD WILL BLOW UP.
# A TOTAL OF 23 HOURS OF DEBUGGING WAS WASTED HERE.
include ../configure.wrf

# MUST use DMPARALLEL=1
ifndef DMPARALLEL
  $(error WRF-GCHP only supports the dmpar parallel option as it requires mvapich2.)
endif

ifneq ($(DMPARALLEL), 1)
  $(error WRF-GCHP only supports the dmpar parallel option as it requires mvapich2.)
endif

# Requirements for GIGC Compilation
# Compiler for ESMF - eg intel, intelgcc, gfortran
ifndef ESMF_COMPILER
  $(error ESMF_COMPILER is not defined - use gfortran or intel)
endif

# MPI type for ESMF
ifndef ESMF_COMM
  $(error ESMF_COMM is not defined - use mvapich2 or openmpi)
endif

# Configuration variables for GCHP interop.
# To ease the burden on the user, we port the options from WRF to here
export GC_BIN=$(NETCDFPATH)/bin
export GC_INCLUDE=$(NETCDFPATH)/include
export GC_LIB=$(NETCDFPATH)/lib

export MPI_ROOT=$(NETCDFPATH)/

export NETCDF_FORTRAN_HOME=$(NETCDFPATH)/
export NETCDF_FORTRAN_INCLUDE=$(NETCDFPATH)/include
export NETCDF_FORTRAN_LIB=$(NETCDFPATH)/lib

export GC_F_BIN=$(GC_BIN)
export GC_F_INCLUDE=$(GC_INCLUDE)
export GC_F_LIB=$(GC_LIB)

export OMPI_FC=$(SFC)
export OMPI_CC=$(SCC)

# Specify compiler-specific options (hplin, 6/23/19)
ifeq ($(ESMF_COMPILER),intel)
	export OMPI_CXX=icpc
	TRACEBACK_OPT := -traceback
else ifeq ($(ESMF_COMPILER),gfortran)
	export OMPI_CXX=g++
	TRACEBACK_OPT := -fbacktrace
endif

# Specify MPI-specific options (hplin, 6/23/19)
ifeq ($(ESMF_COMM),openmpi)
	MPI_OPT := $(shell mpif90 --showme:link)
	MPI_OPT += $(shell mpicxx --showme:link)
else ifeq ($(ESMF_COMM),mvapich2)
	MPI_OPT := -lmpich -lmpichf90
else ifeq ($(ESMF_COMM),intelmpi)
	MPI_OPT := -lmpi
else
	$(error Unknown MPI communicator ESMF_COMM, valid are openmpi or mvapich2)
endif


export COMPILER=$(SFC)
export HDF5DIR=$(HDF5PATH)

# Do not declare / override compilers here - instead pass them to compile_chem ONLY
# otherwise you will be breaking the includes in ./compile em_real (up in WRF)

# .ONESHELL has proven to be problematic...
# DO NOT USE "export" either, it will completely crash GMAO_gfio's compile in MAPL/ESMA

clean:
	rm -f *.o
	rm -f *.f90
	rm -f *.mod
	rm -f gc/GCHP/*.o
	cd gc; FC=$(SFC) CC=$(SCC) F77=$(SFC) F90=$(SFC) CXX=$(OMPI_CXX) make HPC=y realclean
	cd gc/GCHP; FC=$(SFC) CC=$(SCC) F77=$(SFC) F90=$(SFC) CXX=$(OMPI_CXX) make USE_EXTERNAL_GRID=y MET=geos-fp the_nuclear_option
	cd gc/GCHP; FC=$(SFC) CC=$(SCC) F77=$(SFC) F90=$(SFC) CXX=$(OMPI_CXX) make clean
	rm -f ../dyn_em/module_convtrans_prep.f90
	rm -f ../dyn_em/module_convtrans_prep.o
	@echo "Cleaning chem may not be enough - check subdirectories."

softclean:
	@echo "Softclean only cleans the chemistry bindings."
	rm -f *.o
	rm -f *.f90
	rm -f *.mod
	rm -f *.G
	rm -f ../dyn_em/module_convtrans_prep.f90
	rm -f ../dyn_em/module_convtrans_prep.o
	@echo "To run a full clean, run make clean and check subdirectories."

devclean:
	@echo "This is for development purposes only.\n"
	git pull origin master
	rm -f *.f90
	rm -f *.o
	rm -f *.mod
	git checkout -- .
	@echo "Done"

install_common:
	@echo "*****************************************************************"
	@echo "       __          _______  ______       _____  _____            "
	@echo "       \ \        / /  __ \|  ____|     / ____|/ ____|           "
	@echo "        \ \  /\  / /| |__) | |__ ______| |  __| |                "
	@echo "         \ \/  \/ / |  _  /|  __|______| | |_ | |                "
	@echo "          \  /\  /  | | \ \| |         | |__| | |____            "
	@echo "           \/  \/   |_|  \_\_|          \_____|\_____|           "
	@echo "*****************************************************************"
	@echo "   THIS IS THE WRF-GC 'PUMPKIN' CHEMISTRY ABSTRACTION LAYER      "
	@echo "                  FOR THE WRF MODEL VERSION 3+                   "
	@echo "*****************************************************************"
	@echo "THIS WILL INSTALL AND REPLACE THE WRF-CHEM STANDARD REGISTRY.    "
	mv -n ../Registry/registry.chem ../Registry/registry.chem.bak
	cp ./registry.chem ../Registry/registry.chem
	@echo "Your original registry is now in registry.chem.bak.              "
	@echo "To compile, return to root directory and run ./compile em_real   "
	@echo "*****************************************************************"
	@echo "THIS WILL REPLACE PHYS/ WITH TWO-WAY COUPLING DATA PARAM FILES (for WRFv3) "
	mv ../phys/module_data_gocart_dust.F ../phys/module_data_gocart_dust.F.bak
	mv ../share/module_chem_share.F ../share/module_chem_share.F.bak
	cp ./module_data_gocart_dust.cpy ../phys/module_data_gocart_dust.F
	cp ./module_chem_share.cpy ../share/module_chem_share.F
	sed -i -e "s@  integer                , parameter      :: MaxVars          = 3000@  integer                , parameter      :: MaxVars          = 10000@" "$(ROOTDIR)/external/io_netcdf/wrf_io.F90";
	# source common bash functions from scripts 
	source ./newUserRegistration.sh && registerNewUserwrapper

install_registry: install_common
	rm ./wrfgc_convert_state_mod.F
	cp -v ./config/clim_p_trop.nc ../run/
	cp -v ./config/species_database.yml ../run/
	cp -v ./config/fullchem/* ../run/
	ln -sf ./config/couplers/wrfgc_convert_state_mod.F.fullchem ./wrfgc_convert_state_mod.F
	@echo "Configuration files copied successfully to your default WRF rundir"
	@echo "If you are using a custom run directory, copy the files above for"
	@echo "correct GEOS-Chem operation."

install_registry_ch4: install_common
	rm ./wrfgc_convert_state_mod.F
	cp -v ./config/clim_p_trop.nc ../run/
	cp -v ./config/species_database.yml ../run/
	cp -v ./config/ch4/* ../run/
	cp -v ./config/couplers/wrfgc_convert_state_mod.F.ch4 ./wrfgc_convert_state_mod.F
	@echo "INSTALLED WRF-GC, CH4 specialty simulation"
	@echo "Configuration files copied successfully to your default WRF rundir"
	@echo "If you are using a custom run directory, copy the files above for"
	@echo "correct GEOS-Chem operation."

install_registry_co2: install_common
	rm ./wrfgc_convert_state_mod.F
	cp -v ./config/clim_p_trop.nc ../run/
	cp -v ./config/species_database.yml ../run/
	cp -v ./config/co2/* ../run/
	cp -v ./config/couplers/wrfgc_convert_state_mod.F.co2 ./wrfgc_convert_state_mod.F
	@echo "INSTALLED WRF-GC, CO2 specialty simulation"
	@echo "Configuration files copied successfully to your default WRF rundir"
	@echo "If you are using a custom run directory, copy the files above for"
	@echo "correct GEOS-Chem operation."

test_registry_dummy:
	test -s ../Registry/registry.chem.bak || { echo "***** WRF-GC REGISTRY NOT FOUND: DO ./clean -a, cd chem/, make install_registry, THEN COME BACK AGAIN. *****"; false; }

compile_chem: test_registry_dummy
	@echo "*****************************************************************"
	@echo "       __          _______  ______       _____  _____            "
	@echo "       \ \        / /  __ \|  ____|     / ____|/ ____|           "
	@echo "        \ \  /\  / /| |__) | |__ ______| |  __| |                "
	@echo "         \ \/  \/ / |  _  /|  __|______| | |_ | |                "
	@echo "          \  /\  /  | | \ \| |         | |__| | |____            "
	@echo "           \/  \/   |_|  \_\_|          \_____|\_____|           "
	@echo "*****************************************************************"
	@echo "            THIS IS THE WRF-GC PROJECT FOR WRFV4                 "
	@echo "                  BASED ON GEOS-CHEM HP v12.0.0                  "
	@echo "      GEOS-Chem Version: 14.1.1      WRF-GC 2023/03/08           "
	@echo "*****************************************************************"
	@echo " (c) 2018-2023 Haipeng Lin, Xu Feng, Tzung-May Fu*               "
	@echo "*****************************************************************"
	@echo "THIS COMMAND WILL COMPILE GEOS-CHEM IN THE gc SUBDIRECTORY       "

	# if you want to compile with LUO_WETDEP, add LUO_WETDEP=y here
	cd $(ROOTDIR)/chem/gc && FC=$(SFC) CC=$(SCC) F77=$(SFC) F90=$(SFC) CXX=$(OMPI_CXX) $(MAKE) $(J) NC_DIAG=y CHEM=fullchem EXTERNAL_GRID=y DEBUG=n TRACEBACK=y BOUNDS=y MET=geos-fp GRID=4x5 NO_REDUCED=y UCX=yes NO_EXE=y hpc

	@echo "*****************************************************************"
	@echo " WE WILL NOW UPDATE WRF/main EM_REAL COMPILE RULES               "
	@echo " TO LINK AGAINST GEOS-CHEM LIBRARIES. "
	@echo "        /everybody stand back/ I know regular expressions.       "
	@echo "       - The original Makefile is at main/Makefile.bak -         "
	if   grep -q DESMF_ "$(ROOTDIR)/main/Makefile"; then cp $(ROOTDIR)/main/Makefile.bak $(ROOTDIR)/main/Makefile; fi
	if   grep -q chem/gigc "$(ROOTDIR)/main/Makefile"; then cp $(ROOTDIR)/main/Makefile.bak $(ROOTDIR)/main/Makefile; fi
	if   grep -q Isoropia "$(ROOTDIR)/main/Makefile"; then cp $(ROOTDIR)/main/Makefile.bak $(ROOTDIR)/main/Makefile; fi
	if ! grep -q DMODEL_WRF "$(ROOTDIR)/main/Makefile"; then cp $(ROOTDIR)/main/Makefile $(ROOTDIR)/main/Makefile.bak; sed -i -e "s@-o wrf\.exe .*@& -DLINUX_IFORT -DEXTERNAL_GRID -DNC_DIAG -DUCX -DGEOS_FP -DNC_HAS_COMPRESSION -DMODEL_ -DMODEL_WRF -DUSE_REAL8 -lGIGC $(MPI_OPT) -L$(ROOTDIR)/chem/gc/lib -lHistory -lGeosCore -lHistory -lKpp -lGeosCore -lIsorropia -lHCOI -lHCOX -lHCO -lGeosUtil -lKpp -lHeaders -lNcUtils -L$(NETCDFPATH) -lnetcdff -lnetcdf -L$(HDF5PATH)/lib -I$(ROOTDIR)/chem/gc/mod -lnetcdf@" "$(ROOTDIR)/main/Makefile"; sed -i -e "s@-o real\.exe .*@& -DLINUX_IFORT -DEXTERNAL_GRID -DNC_DIAG -DUCX -DGEOS_FP -DNC_HAS_COMPRESSION -DMODEL_ -DMODEL_WRF -DUSE_REAL8 -lGIGC $(MPI_OPT) -L$(ROOTDIR)/chem/gc/lib -lHistory -lGeosCore -lHistory -lKpp -lGeosCore -lIsorropia -lHCOI -lHCOX -lHCO -lGeosUtil -lKpp -lHeaders -lNcUtils -L$(NETCDFPATH) -lnetcdff -lnetcdf -L$(HDF5PATH)/lib -I$(ROOTDIR)/chem/gc/mod -lnetcdf@" "$(ROOTDIR)/main/Makefile"; sed -i -e "s@-o nup\.exe .*@& -DLINUX_IFORT -DEXTERNAL_GRID -DNC_DIAG -DUCX -DGEOS_FP -DNC_HAS_COMPRESSION -DMODEL_ -DMODEL_WRF -DUSE_REAL8 -lGIGC $(MPI_OPT) -L$(ROOTDIR)/chem/gc/lib -lHistory -lGeosCore -lHistory -lKpp -lGeosCore -lIsorropia -lHCOI -lHCOX -lHCO -lGeosUtil -lKpp -lHeaders -lNcUtils -L$(NETCDFPATH) -lnetcdff -lnetcdf -L$(HDF5PATH)/lib -I$(ROOTDIR)/chem/gc/mod -lnetcdf@" "$(ROOTDIR)/main/Makefile"; sed -i -e "s@-o ndown\.exe .*@& -DLINUX_IFORT -DEXTERNAL_GRID -DNC_DIAG -DUCX -DGEOS_FP -DNC_HAS_COMPRESSION -DMODEL_ -DMODEL_WRF -DUSE_REAL8 -lGIGC $(MPI_OPT) -L$(ROOTDIR)/chem/gc/lib -lHistory -lGeosCore -lHistory -lKpp -lGeosCore -lIsorropia -lHCOI -lHCOX -lHCO -lGeosUtil -lKpp -lHeaders -lNcUtils -L$(NETCDFPATH) -lnetcdff -lnetcdf -L$(HDF5PATH)/lib -I$(ROOTDIR)/chem/gc/mod -lnetcdf@" "$(ROOTDIR)/main/Makefile"; sed -i -e "s@-o tc\.exe .*@& -DLINUX_IFORT -DEXTERNAL_GRID -DNC_DIAG -DUCX -DGEOS_FP -DNC_HAS_COMPRESSION -DMODEL_ -DMODEL_WRF -DUSE_REAL8 -lGIGC $(MPI_OPT) -L$(ROOTDIR)/chem/gc/lib -lHistory -lGeosCore -lHistory -lKpp -lGeosCore -lIsorropia -lHCOI -lHCOX -lHCO -lGeosUtil -lKpp -lHeaders -lNcUtils -L$(NETCDFPATH) -lnetcdff -lnetcdf -L$(HDF5PATH)/lib -I$(ROOTDIR)/chem/gc/mod -lnetcdf@" "$(ROOTDIR)/main/Makefile"; fi

	@echo "   *********** GEOS-CHEM HAS BEEN INSTALLED IN WRF ***********   "
	@echo "   Please remember that:                                         "
	@echo "   - Documentation is at wrfgc.readthedocs.io.                   "
	@echo "   - GEOS-Chem output species are different than in WRF-Chem.    "
	@echo "   - You can configure GEOS-Chem in the geoschem_config.yml file,"
	@echo "     but:                                                        "
	@echo "     (1) individual process switches are in namelist.input       "
	@echo "         gc_do_convection, _hemco, _pblmix, _chemistry, _drydep, _wetdep"
	@echo "     (2) simulation types can only be switched with commands     "
	@echo "         make_install_registry{_ch4/_co2} and a quick recompile  "
	@echo "   - Not all HEMCO emissions are supported in all areas.         "
	@echo "   With great power comes great responsibility.                  "
	@echo "   Please refer to 'make about' or contact WRF-GC team.          "


troubleshooting:
	@echo "*****************************************************************"
	@echo "              THIS IS THE WRF-GC PROJECT FOR WRFV3               "
	@echo "                 BASED ON GEOS-CHEM HP v12.0.0                   "
	@echo "*****************************************************************"
	@echo " TROUBLESHOOTING INFORMATION:                                    "
	@echo "   - CC: $(SCC), FC: $(SFC) (WRF)"
	@echo "   - CC: $(CC), FC: $(FC) (GEOS-Chem)"
	@echo "   - Parallel Make (-j N): $(J)"
	@echo "   - MPIROOT: $(MPI_ROOT)"
	@echo "   - OMPI_CC: $(OMPI_CC), _FC: $(OMPI_FC), _CXX: $(OMPI_CXX)"
	@echo "   - ROOTDIR: $(ROOTDIR)"
	@echo "   - WRF_SRC_ROOT_DIR: $(WRF_SRC_ROOT_DIR)"

about:
	@echo "*****************************************************************"
	@echo "       __          _______  ______       _____  _____            "
	@echo "       \ \        / /  __ \|  ____|     / ____|/ ____|           "
	@echo "        \ \  /\  / /| |__) | |__ ______| |  __| |                "
	@echo "         \ \/  \/ / |  _  /|  __|______| | |_ | |                "
	@echo "          \  /\  /  | | \ \| |         | |__| | |____            "
	@echo "           \/  \/   |_|  \_\_|          \_____|\_____|           "
	@echo "*****************************************************************"
	@echo "   THIS IS THE WRF-GCHP 'PUMPKIN' CHEMISTRY ABSTRACTION LAYER    "
	@echo "                    FOR THE WRF MODEL VERSION 4                  "
	@echo "*****************************************************************"
	@echo " FOR ERRORS, SUGGESTIONS AND FEEDBACK, CONTACT HAIPENG LIN AT    "
	@echo "           HPLIN@SEAS.HARVARD.EDU | JIMMIE.LIN@GMAIL.COM         "
	@echo "*****************************************************************"
	@echo " (c) 2018-2023 Haipeng Lin, Xu Feng, Tzung-May Fu*               "
	@echo " Atmospheric Chemistry and Climate Group, SUSTech                "
	@echo "*****************************************************************"
	@echo "Commands:                                                        "
	@echo "    make about - Show this about (help) screen                   "
	@echo "    make clean - Clean chemistry (full, use with caution)        "
	@echo "    make softclean - Soft clean 'Pumpkin' bindings only          "
	@echo "    make install_registry - Install chemistry species into WRF   "
	@echo "       (replacing existing Registry.chem) and configs into run   "
	@echo "    make compile_chem - Compile target chemistry                 "
	@echo "    make troubleshooting - Show debugging configuration info     "
	@echo "    make devclean - Dev purposes only, get from git origin/master"
	@echo "                                                                 "
	@echo "  * Not all commands, especially install_registry, compile_chem  "
	@echo "    are always available. Check with your chemistry maintainer.  "
	@echo "    These commands are mostly developed for WRF-GCHP's GC part.  "
	@echo "  * If you compiled WRF without install_registry, you need to    "
	@echo "    clean all (./clean -a) to rebuild module_state_description.  "
	@echo "    A non-clean compile cannot refresh the registry.             "
	@echo "*****************************************************************"

# DEPENDENCIES : only dependencies after this line
module_data_radm2.o:

module_data_sorgam.o: module_data_radm2.o

# Required by module_cu_kfcup
module_data_mosaic_asect.o:

module_chem_utilities.o:

module_convection_prep.o:

module_interpolate.o:

module_vertmx_wrf.o:

module_aer_opt_out.o:

module_tropopause.o: module_interpolate.o

module_upper_bc_driver.o: module_tropopause.o

module_input_chem_bioemiss.o:

module_mosaic_driver.o:

module_aerosols_sorgam.o:

module_aerosols_soa_vbs.o:

module_input_chem_data.o: module_data_sorgam.o

# Required by optical_driver
module_data_rrtmgaeropt.o:

module_data_mosaic_other.o:

module_data_mosaic_therm.o:

module_peg_util.o:

module_optical_averaging.o: ../phys/module_data_gocart_dust.o ../phys/module_cam_infnan.o module_data_rrtmgaeropt.o module_peg_util.o module_data_sorgam.o module_data_mosaic_asect.o module_data_mosaic_therm.o module_data_mosaic_other.o

optical_driver.o: module_optical_averaging.o module_data_rrtmgaeropt.o module_peg_util.o

# Required by module_diag_aero_size_info
module_data_gigc_asect.o:

module_diag_aero_size_info.o: ../phys/module_data_gocart_dust.o module_data_gigc_asect.o module_data_sorgam.o

module_mixactivate_wrappers.o: ../phys/module_mixactivate.o module_data_gigc_asect.o

mixactivate_driver.o: module_mixactivate_wrappers.o

wrfgc_history_mod.o: compile_chem 
	$(RM) $@
	$(CPP) -I$(WRF_SRC_ROOT_DIR)/inc $(CPPFLAGS) $(OMPCPP) $*.F > $*.G
	$(SED_FTN) $*.G | $(CPP) $(TRADFLAG) > $*.f90
	$(RM) $*.G
	@ if echo $(ARCHFLAGS) | $(FGREP) 'DVAR4D'; then \
          echo COMPILING $*.F for 4DVAR ; \
          $(WRF_SRC_ROOT_DIR)/var/build/da_name_space.pl $*.f90 > $*.f90.tmp ; \
          mv $*.f90.tmp $*.f90 ; \
        fi
	$(FC) -o $@ -c $(FCFLAGS) -check bounds $(TRACEBACK_OPT) $(OMP) $(MODULE_DIRS) -DLINUX_IFORT -DEXTERNAL_GRID -DNC_DIAG -DUCX -DGEOS_FP -DNC_HAS_COMPRESSION -DMODEL_ -DMODEL_WRF -DUSE_REAL8 -L$(PNETCDF)/lib -I$(PNETCDF)/include -lpnetcdf -I$(ROOTDIR)/chem/gc/mod $(PROMOTION) $(FCSUFFIX) $*.f90

wrfgc_io_pnetcdf.o: compile_chem wrfgc_history_mod.o
	$(RM) $@
	$(CPP) -I$(WRF_SRC_ROOT_DIR)/inc $(CPPFLAGS) $(OMPCPP) $*.F > $*.G
	$(SED_FTN) $*.G | $(CPP) $(TRADFLAG) > $*.f90
	$(RM) $*.G
	@ if echo $(ARCHFLAGS) | $(FGREP) 'DVAR4D'; then \
          echo COMPILING $*.F for 4DVAR ; \
          $(WRF_SRC_ROOT_DIR)/var/build/da_name_space.pl $*.f90 > $*.f90.tmp ; \
          mv $*.f90.tmp $*.f90 ; \
        fi
	$(FC) -o $@ -c $(FCFLAGS) -check bounds $(TRACEBACK_OPT) $(OMP) $(MODULE_DIRS) -DLINUX_IFORT -DEXTERNAL_GRID -DNC_DIAG -DUCX -DGEOS_FP -DNC_HAS_COMPRESSION -DMODEL_ -DMODEL_WRF -DUSE_REAL8 -L$(PNETCDF)/lib -I$(PNETCDF)/include -lpnetcdf -I$(ROOTDIR)/chem/gc/mod $(PROMOTION) $(FCSUFFIX) $*.f90

CONVERT_STATE_RELA = 	compile_chem \
						module_input_chem_data.o \
						module_tropopause.o \
						module_upper_bc_driver.o 
CHEMICS_INIT_RELA = 	compile_chem \
						module_input_chem_data.o \
						module_mixactivate_wrappers.o \
						module_diag_aero_size_info.o \
						module_tropopause.o \
						module_upper_bc_driver.o \
						wrfgc_convert_state_mod.o 

ifdef PNETCDF
CHEMICS_INIT_RELA += wrfgc_history_mod.o
CONVERT_STATE_RELA += wrfgc_io_pnetcdf.o
PNETCDFFLAG = -L$(PNETCDF)/lib -I$(PNETCDF)/include -lpnetcdf
WRFGC_HISTORY_FLAG = -Duse_wrfgc_history_output
endif

wrfgc_convert_state_mod.o: $(CONVERT_STATE_RELA)
	$(RM) $@
	$(CPP) -I$(WRF_SRC_ROOT_DIR)/inc $(CPPFLAGS) $(OMPCPP) $(WRFGC_HISTORY_FLAG) $*.F > $*.G
	$(SED_FTN) $*.G | $(CPP) $(TRADFLAG) > $*.f90
	$(RM) $*.G
	@ if echo $(ARCHFLAGS) | $(FGREP) 'DVAR4D'; then \
          echo COMPILING $*.F for 4DVAR ; \
          $(WRF_SRC_ROOT_DIR)/var/build/da_name_space.pl $*.f90 > $*.f90.tmp ; \
          mv $*.f90.tmp $*.f90 ; \
        fi
	$(FC) -o $@ -c $(FCFLAGS) -check bounds $(TRACEBACK_OPT) $(OMP) $(MODULE_DIRS) -DLINUX_IFORT -DEXTERNAL_GRID -DNC_DIAG -DUCX -DGEOS_FP -DNC_HAS_COMPRESSION -DMODEL_ -DMODEL_WRF -DUSE_REAL8 -L$(PNETCDF)/lib -I$(PNETCDF)/include -lpnetcdf -I$(ROOTDIR)/chem/gc/mod $(PROMOTION) $(FCSUFFIX) $*.f90

chemics_init.o: $(CHEMICS_INIT_RELA)
	$(RM) $@
	sed -e "s/grid%mu/gridmu/g" -e "s/grid%Mu/gridMu/g" -e "s/^\!.*'.*//" -e "s/^ *\!.*'.*//" $*.F > $*.G
	$(CPP) -I$(WRF_SRC_ROOT_DIR)/inc $(CPPFLAGS) $(OMPCPP)  $(WRFGC_HISTORY_FLAG) $*.G  > $*.H
	sed -e "s/gridmu/grid%mu/g" -e "s/gridMu/grid%Mu/g" $*.H > $*.bb
	$(SED_FTN) $*.bb | $(CPP) $(TRADFLAG) > $*.f90
	$(RM) $*.G $*.H $*.bb
	@ if echo $(ARCHFLAGS) | $(FGREP) 'DVAR4D'; then \
          echo COMPILING $*.F for 4DVAR ; \
          $(WRF_SRC_ROOT_DIR)/var/build/da_name_space.pl $*.f90 > $*.f90.tmp ; \
          mv $*.f90.tmp $*.f90 ; \
        fi
	$(FC) -o $@ -c $(FCFLAGS) -check bounds $(TRACEBACK_OPT) $(OMP) $(MODULE_DIRS) -DLINUX_IFORT -DEXTERNAL_GRID -DNC_DIAG -DUCX -DGEOS_FP -DNC_HAS_COMPRESSION -DMODEL_ -DMODEL_WRF -DUSE_REAL8 -I$(ROOTDIR)/chem/gc/mod $(PROMOTION) $(FCSUFFIX) $*.f90

chem_driver.o: compile_chem ../dyn_em/module_convtrans_prep.o module_input_chem_data.o module_chem_utilities.o module_tropopause.o  module_upper_bc_driver.o wrfgc_convert_state_mod.o optical_driver.o module_diag_aero_size_info.o mixactivate_driver.o
	$(RM) $@
	sed -e "s/grid%mu/gridmu/g" -e "s/grid%Mu/gridMu/g" -e "s/^\!.*'.*//" -e "s/^ *\!.*'.*//" $*.F > $*.G
	$(CPP) -I$(WRF_SRC_ROOT_DIR)/inc $(CPPFLAGS) $(OMPCPP) $*.G  > $*.H
	sed -e "s/gridmu/grid%mu/g" -e "s/gridMu/grid%Mu/g" $*.H > $*.bb
	$(SED_FTN) $*.bb | $(CPP) $(TRADFLAG) > $*.f90
	$(RM) $*.G $*.H $*.bb
	@ if echo $(ARCHFLAGS) | $(FGREP) 'DVAR4D'; then \
          echo COMPILING $*.F for 4DVAR ; \
          $(WRF_SRC_ROOT_DIR)/var/build/da_name_space.pl $*.f90 > $*.f90.tmp ; \
          mv $*.f90.tmp $*.f90 ; \
        fi
	$(FC) -o $@ -c $(FCFLAGS) -check bounds $(TRACEBACK_OPT) $(OMP) $(MODULE_DIRS) -DLINUX_IFORT -DEXTERNAL_GRID -DNC_DIAG -DUCX -DGEOS_FP -DNC_HAS_COMPRESSION -DMODEL_ -DMODEL_WRF -DUSE_REAL8 -I$(ROOTDIR)/chem/gc/mod $(PROMOTION) $(FCSUFFIX) $*.f90

# Note that linking is suppressed for chemics & chem_driver, so you have to do this linking in wrf.exe, ndown/nup.exe, ...
