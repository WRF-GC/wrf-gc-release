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
        gigc_convert_state_mod.o

OBJS    =                           \
        chemics_init.o              \
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
.PHONY: chemics clean devclean install_registry install_configs compile_chem troubleshooting

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
else ifeq($(ESMF_COMPILER),gfortran)
	export OMPI_CXX=g++
	TRACEBACK_OPT := -fbacktrace
endif

# Specify MPI-specific options (hplin, 6/23/19)
ifeq ($(ESMF_COMM),openmpi)
	MPI_OPT := $(shell mpif90 --showme:link)
	MPI_OPT += $(shell mpicxx --showme:link)
else ifeq ($(ESMF_COMM),mvapich2)
	MPI_OPT := -lmpich -lmpichf90
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
	cd gigc; make HPC=y realclean
	cd gigc/GCHP; make USE_EXTERNAL_GRID=y MET=geos-fp the_nuclear_option
	cd gigc/GCHP; make clean
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

install_registry:
	@echo "*****************************************************************"
	@echo "  __          _______  ______       _____  _____ _    _ _____    "
	@echo "  \ \        / /  __ \|  ____|     / ____|/ ____| |  | |  __ \   "
	@echo "   \ \  /\  / /| |__) | |__ ______| |  __| |    | |__| | |__) |  "
	@echo "    \ \/  \/ / |  _  /|  __|______| | |_ | |    |  __  |  ___/   "
	@echo "     \  /\  /  | | \ \| |         | |__| | |____| |  | | |       "
	@echo "      \/  \/   |_|  \_\_|          \_____|\_____|_|  |_|_|       "
	@echo "*****************************************************************"
	@echo "   THIS IS THE WRF-GC 'PUMPKIN' CHEMISTRY ABSTRACTION LAYER      "
	@echo "                  FOR THE WRF MODEL VERSION 3+                   "
	@echo "*****************************************************************"
	@echo " (c) 2018 Haipeng Lin                                            "
	@echo " Peking University, Atmospheric Chemistry and Climate Group      "
	@echo "*****************************************************************"
	@echo "THIS WILL INSTALL AND REPLACE THE WRF-CHEM STANDARD REGISTRY.    "
	mv ../Registry/registry.chem ../Registry/registry.chem.bak
	cp ./registry.chem ../Registry/registry.chem
	@echo "Your original registry is now in registry.chem.bak.              "
	@echo "To compile, return to root directory and run ./compile em_real   "

install_configs:
	cp -v ./config/* ../run/
	@echo "Configuration files copied successfully to your default WRF rundir"
	@echo "If you are using a custom run directory, copy the files above for"
	@echo "correct GEOS-Chem operation."

compile_chem: install_registry install_configs
	@echo "*****************************************************************"
	@echo "  __          _______  ______       _____  _____ _    _ _____    "
	@echo "  \ \        / /  __ \|  ____|     / ____|/ ____| |  | |  __ \   "
	@echo "   \ \  /\  / /| |__) | |__ ______| |  __| |    | |__| | |__) |  "
	@echo "    \ \/  \/ / |  _  /|  __|______| | |_ | |    |  __  |  ___/   "
	@echo "     \  /\  /  | | \ \| |         | |__| | |____| |  | | |       "
	@echo "      \/  \/   |_|  \_\_|          \_____|\_____|_|  |_|_|       "
	@echo "*****************************************************************"
	@echo "            THIS IS THE WRF-GC PROJECT FOR WRFV3+                "
	@echo "                  BASED ON GEOS-CHEM HP v12.0.0                  "
	@echo "*****************************************************************"
	@echo " (c) 2018 Haipeng Lin, Xu Feng, Tzung-May Fu*                    "
	@echo " Peking University, Atmospheric Chemistry and Climate Group      "
	@echo "*****************************************************************"
	@echo "THIS COMMAND WILL PRE-COMPILE GEOS-CHEM IN THE gigc SUBDIRECTORY "
	cd $(ROOTDIR)/chem/gigc && FC=$(SFC) CC=$(SCC) F77=$(SFC) F90=$(SFC) CXX=$(OMPI_CXX) $(MAKE) $(J) NC_DIAG=y CHEM=Standard EXTERNAL_GRID=y DEBUG=n TRACEBACK=y MET=geos-fp GRID=4x5 NO_REDUCED=y UCX=yes NO_EXE=y hpc

	@echo "*****************************************************************"
	@echo " WE WILL NOW UPDATE WRFV3/main EM_REAL COMPILE RULES             "
	@echo " TO LINK AGAINST GEOS-CHEM LIBRARIES. "
	@echo "        /everybody stand back/ I know regular expressions.       "
	@echo "       - The original Makefile is at main/Makefile.bak -         "
	if   grep -q DESMF_ "$(ROOTDIR)/main/Makefile"; then cp $(ROOTDIR)/main/Makefile.bak $(ROOTDIR)/main/Makefile; fi
	if ! grep -q DMODEL_WRF "$(ROOTDIR)/main/Makefile"; then cp $(ROOTDIR)/main/Makefile $(ROOTDIR)/main/Makefile.bak; sed -i -e "s@-o wrf\.exe .*@& -DLINUX_IFORT -DEXTERNAL_GRID -DNC_DIAG -DUCX -DGEOS_FP -DNC_HAS_COMPRESSION -DMODEL_ -DMODEL_WRF -DUSE_REAL8 -lGIGC $(MPI_OPT) -L$(ROOTDIR)/chem/gigc/lib -lHistory -lGeosCore -lHistory -lKpp -lGeosCore -lIsoropia -lHCOI -lHCOX -lHCO -lGeosUtil -lKpp -lHeaders -lNcUtils -L$(NETCDFPATH) -lnetcdff -lnetcdf -L$(HDF5PATH)/lib -I$(ROOTDIR)/chem/gigc/mod -lnetcdf@" "$(ROOTDIR)/main/Makefile"; sed -i -e "s@-o real\.exe .*@& -DLINUX_IFORT -DEXTERNAL_GRID -DNC_DIAG -DUCX -DGEOS_FP -DNC_HAS_COMPRESSION -DMODEL_ -DMODEL_WRF -DUSE_REAL8 -lGIGC $(MPI_OPT) -L$(ROOTDIR)/chem/gigc/lib -lHistory -lGeosCore -lHistory -lKpp -lGeosCore -lIsoropia -lHCOI -lHCOX -lHCO -lGeosUtil -lKpp -lHeaders -lNcUtils -L$(NETCDFPATH) -lnetcdff -lnetcdf -L$(HDF5PATH)/lib -I$(ROOTDIR)/chem/gigc/mod -lnetcdf@" "$(ROOTDIR)/main/Makefile"; sed -i -e "s@-o nup\.exe .*@& -DLINUX_IFORT -DEXTERNAL_GRID -DNC_DIAG -DUCX -DGEOS_FP -DNC_HAS_COMPRESSION -DMODEL_ -DMODEL_WRF -DUSE_REAL8 -lGIGC $(MPI_OPT) -L$(ROOTDIR)/chem/gigc/lib -lHistory -lGeosCore -lHistory -lKpp -lGeosCore -lIsoropia -lHCOI -lHCOX -lHCO -lGeosUtil -lKpp -lHeaders -lNcUtils -L$(NETCDFPATH) -lnetcdff -lnetcdf -L$(HDF5PATH)/lib -I$(ROOTDIR)/chem/gigc/mod -lnetcdf@" "$(ROOTDIR)/main/Makefile"; sed -i -e "s@-o ndown\.exe .*@& -DLINUX_IFORT -DEXTERNAL_GRID -DNC_DIAG -DUCX -DGEOS_FP -DNC_HAS_COMPRESSION -DMODEL_ -DMODEL_WRF -DUSE_REAL8 -lGIGC $(MPI_OPT) -L$(ROOTDIR)/chem/gigc/lib -lHistory -lGeosCore -lHistory -lKpp -lGeosCore -lIsoropia -lHCOI -lHCOX -lHCO -lGeosUtil -lKpp -lHeaders -lNcUtils -L$(NETCDFPATH) -lnetcdff -lnetcdf -L$(HDF5PATH)/lib -I$(ROOTDIR)/chem/gigc/mod -lnetcdf@" "$(ROOTDIR)/main/Makefile"; sed -i -e "s@-o tc\.exe .*@& -DLINUX_IFORT -DEXTERNAL_GRID -DNC_DIAG -DUCX -DGEOS_FP -DNC_HAS_COMPRESSION -DMODEL_ -DMODEL_WRF -DUSE_REAL8 -lGIGC $(MPI_OPT) -L$(ROOTDIR)/chem/gigc/lib -lHistory -lGeosCore -lHistory -lKpp -lGeosCore -lIsoropia -lHCOI -lHCOX -lHCO -lGeosUtil -lKpp -lHeaders -lNcUtils -L$(NETCDFPATH) -lnetcdff -lnetcdf -L$(HDF5PATH)/lib -I$(ROOTDIR)/chem/gigc/mod -lnetcdf@" "$(ROOTDIR)/main/Makefile"; fi

	@echo "   *********** GEOS-CHEM HAS BEEN INSTALLED IN WRF ***********   "
	@echo "   Please remember that:"
	@echo "   - WRF-GC is a highly experimental project."
	@echo "   - GEOS-Chem output species are different than in WRF-Chem."
	@echo "   - You can configure GEOS-Chem in the input.geos file, except"
	@echo "     individual process switches are in namelist.input as follows"
	@echo "        gc_do_convection, _hemco, _pblmix, _chemistry, _drydep, _wetdep"
	@echo "   - Not all HEMCO emissions are supported in all areas."
	@echo "   - You have to use WRF-Chem style emissions, if you want to"
	@echo "     run in multiple-domains. Turn on wrfgc_legacy_emis, and use"
	@echo "     regular anthropogenic emissions inventory as WRF-Chem."
	@echo "   With great power comes great responsibility."
	@echo "   Please refer to 'make about' or contact WRF-GC team."


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
	@echo "  __          _______  ______       _____  _____ _    _ _____    "
	@echo "  \ \        / /  __ \|  ____|     / ____|/ ____| |  | |  __ \   "
	@echo "   \ \  /\  / /| |__) | |__ ______| |  __| |    | |__| | |__) |  "
	@echo "    \ \/  \/ / |  _  /|  __|______| | |_ | |    |  __  |  ___/   "
	@echo "     \  /\  /  | | \ \| |         | |__| | |____| |  | | |       "
	@echo "      \/  \/   |_|  \_\_|          \_____|\_____|_|  |_|_|       "
	@echo "*****************************************************************"
	@echo "   THIS IS THE WRF-GCHP 'PUMPKIN' CHEMISTRY ABSTRACTION LAYER    "
	@echo "                    FOR THE WRF MODEL VERSION 3                  "
	@echo "*****************************************************************"
	@echo " FOR ERRORS, SUGGESTIONS AND FEEDBACK, CONTACT HAIPENG LIN AT    "
	@echo "           LINHAIPENG@PKU.EDU.CN | JIMMIE.LIN@GMAIL.COM          "
	@echo "*****************************************************************"
	@echo " (c) 2018 Haipeng Lin                                            "
	@echo " Peking University, Atmospheric Chemistry and Climate Group      "
	@echo "*****************************************************************"
	@echo "Commands:                                                        "
	@echo "    make about - Show this about (help) screen                   "
	@echo "    make clean - Clean chemistry (full, use with caution)        "
	@echo "    make softclean - Soft clean 'Pumpkin' bindings only          "
	@echo "    make install_registry - Install chemistry species into WRF   "
	@echo "       (replacing existing Registry.chem)                        "
	@echo "    make install_configs - Copies config files into run directory"
	@echo "    make compile_chem - Compile target chemistry                 "
	@echo "    make troubleshooting - Show debugging configuration info     "
	@echo "    make devclean - Dev purposes only, get from git origin/master"
	@echo "                                                                 "
	@echo "  * Not all commands, especially install_registry, compile_chem  "
	@echo "    are always available. Check with your chemistry maintainer.  "
	@echo "    These commands are mostly developed for WRF-GCHP's GIGC part."
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

gigc_convert_state_mod.o: compile_chem module_input_chem_data.o module_tropopause.o module_upper_bc_driver.o
	$(RM) $@
	$(CPP) -I$(WRF_SRC_ROOT_DIR)/inc $(CPPFLAGS) $(OMPCPP) $*.F > $*.G
	$(SED_FTN) $*.G | $(CPP) $(TRADFLAG) > $*.f90
	$(RM) $*.G
	@ if echo $(ARCHFLAGS) | $(FGREP) 'DVAR4D'; then \
          echo COMPILING $*.F for 4DVAR ; \
          $(WRF_SRC_ROOT_DIR)/var/build/da_name_space.pl $*.f90 > $*.f90.tmp ; \
          mv $*.f90.tmp $*.f90 ; \
        fi
	$(FC) -o $@ -c $(FCFLAGS) $(TRACEBACK_OPT) $(OMP) $(MODULE_DIRS) -DLINUX_IFORT -DEXTERNAL_GRID -DNC_DIAG -DUCX -DGEOS_FP -DNC_HAS_COMPRESSION -DMODEL_ -DMODEL_WRF -DUSE_REAL8 -I$(ROOTDIR)/chem/gigc/mod $(PROMOTION) $(FCSUFFIX) $*.f90

chemics_init.o: compile_chem module_input_chem_data.o module_tropopause.o module_upper_bc_driver.o gigc_convert_state_mod.o
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
	$(FC) -o $@ -c $(FCFLAGS) $(TRACEBACK_OPT) $(OMP) $(MODULE_DIRS) -DLINUX_IFORT -DEXTERNAL_GRID -DNC_DIAG -DUCX -DGEOS_FP -DNC_HAS_COMPRESSION -DMODEL_ -DMODEL_WRF -DUSE_REAL8 -I$(ROOTDIR)/chem/gigc/mod $(PROMOTION) $(FCSUFFIX) $*.f90

chem_driver.o: compile_chem ../dyn_em/module_convtrans_prep.o module_input_chem_data.o module_chem_utilities.o module_tropopause.o  module_upper_bc_driver.o gigc_convert_state_mod.o
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
	$(FC) -o $@ -c $(FCFLAGS) $(TRACEBACK_OPT) $(OMP) $(MODULE_DIRS) -DLINUX_IFORT -DEXTERNAL_GRID -DNC_DIAG -DUCX -DGEOS_FP -DNC_HAS_COMPRESSION -DMODEL_ -DMODEL_WRF -DUSE_REAL8 -I$(ROOTDIR)/chem/gigc/mod $(PROMOTION) $(FCSUFFIX) $*.f90

# Note that linking is suppressed for chemics & chem_driver, so you have to do this linking in wrf.exe, ndown/nup.exe, ...
