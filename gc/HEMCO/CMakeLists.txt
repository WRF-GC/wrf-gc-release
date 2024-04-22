# HEMCO/CMakeLists.txt

cmake_minimum_required(VERSION 3.5)
project(HEMCO VERSION 3.6.2 LANGUAGES Fortran)
# Reminder: Make sure to also update version in src/Core/hco_error_mod.F90

#-----------------------------------------------------------------------------
# Add CMakeScripts/ to the module path
#-----------------------------------------------------------------------------
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/CMakeScripts)
include(HCO-Helpers)

#-----------------------------------------------------------------------------
# Print header
#-----------------------------------------------------------------------------
get_repo_version(HCO_REPO_VERSION ${CMAKE_CURRENT_SOURCE_DIR})
message("=================================================================")
message("HEMCO ${PROJECT_VERSION}")
message("Current status: ${HCO_REPO_VERSION}")
message("=================================================================")

#-----------------------------------------------------------------------------
# Set policies
#-----------------------------------------------------------------------------
cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0057 NEW)
if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()
if(POLICY CMP0079)
    cmake_policy(SET CMP0079 NEW)
endif()

#-----------------------------------------------------------------------------
# Declare the HEMCOBuildProperties library
#-----------------------------------------------------------------------------
add_library(HEMCOBuildProperties INTERFACE)

#---------------------------------------------------------------------
# Default compiler options
#
# HEMCO_Fortran_FLAGS_<COMPILER_ID> : common flags (all build types) 
# HEMCO_Fortran_FLAGS_<CONFIG>_<COMPILER_ID> : types-specific flags
#
# Valid COMPILER_ID: Intel and GNU
# Valid CONFIG: RELEASE, DEBUG, RELWITHDEBINFO
#---------------------------------------------------------------------
set(HEMCO_Fortran_FLAGS_Intel
   -cpp -w -auto -noalign "SHELL:-convert big_endian" "SHELL:-fp-model source" -mcmodel=medium
   -shared-intel -traceback -DLINUX_IFORT
   CACHE STRING "HEMCO compiler flags for all build types with Intel compilers"
)
set(HEMCO_Fortran_FLAGS_RELEASE_Intel
   -O2
   CACHE STRING "HEMCO compiler flags for build type release with Intel compilers"
)
set(HEMCO_Fortran_FLAGS_RELWITHDEBINFO_Intel
   -O2
   CACHE STRING "HEMCO compiler flags for build type relwithdebinfo with Intel compilers"
)
set(HEMCO_Fortran_FLAGS_DEBUG_Intel
   -g -O0 "SHELL:-check arg_temp_created" "SHELL:-debug all" -fpe0 -ftrapuv -check,bounds -DDEBUG
   CACHE STRING "HEMCO compiler flags for build type debug with Intel compilers"
)

set(HEMCO_Fortran_FLAGS_GNU
   -cpp -w -std=legacy -fautomatic -fno-align-commons -fconvert=big-endian
   -fno-range-check -mcmodel=medium -fbacktrace -g -DLINUX_GFORTRAN
   -ffree-line-length-none
   CACHE STRING "HEMCO compiler flags for all build types with GNU compilers"
)
set(HEMCO_Fortran_FLAGS_RELEASE_GNU
   -O3 -funroll-loops
   CACHE STRING "HEMCO compiler flags for build type release with GNU compilers"
)
set(HEMCO_Fortran_FLAGS_RELWITHDEBINFO_GNU
   -O3 -funroll-loops
   CACHE STRING "HEMCO compiler flags for build type relwithdebinfo with GNU compilers"
)
set(HEMCO_Fortran_FLAGS_DEBUG_GNU
   -g -gdwarf-2 -gstrict-dwarf -O0 -Wall -Wextra -Wconversion -Warray-temporaries 
   -fcheck=array-temps -ffpe-trap=invalid,zero,overflow -finit-real=snan -fcheck=bounds -fcheck=pointer
   -fcheck=no-recursion
   CACHE STRING "HEMCO compiler flags for build type debug with GNU compilers"
)

set(HEMCO_SUPPORTED_COMPILER_IDS "Intel" "GNU")
if(NOT CMAKE_Fortran_COMPILER_ID IN_LIST HEMCO_SUPPORTED_COMPILER_IDS)
   message(FATAL_ERROR "HEMCO does not support ${CMAKE_Fortran_COMPILER_ID} compilers")
endif()

#---------------------------------------------------------------------
# Assign comiler options to build properties
#---------------------------------------------------------------------
target_compile_options(HEMCOBuildProperties
   INTERFACE
   $<$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},Intel>:
      ${HEMCO_Fortran_FLAGS_Intel}
      $<$<CONFIG:Debug>:${HEMCO_Fortran_FLAGS_DEBUG_Intel}>
      $<$<CONFIG:RelWithDebInfo>:${HEMCO_Fortran_FLAGS_RELWITHDEBINFO_Intel}>
      $<$<CONFIG:Release>:${HEMCO_Fortran_FLAGS_RELEASE_Intel}>
   >
   $<$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},GNU>:
      ${HEMCO_Fortran_FLAGS_GNU}
      $<$<CONFIG:Debug>:${HEMCO_Fortran_FLAGS_DEBUG_GNU}>
      $<$<CONFIG:RelWithDebInfo>:${HEMCO_Fortran_FLAGS_RELWITHDEBINFO_GNU}>
      $<$<CONFIG:Release>:${HEMCO_Fortran_FLAGS_RELEASE_GNU}>
   >
)

#-----------------------------------------------------------------------------
# Set USE_REAL8 as cache variable to not override the existing definition
# See https://github.com/geoschem/geos-chem/issues/43.
#-----------------------------------------------------------------------------
set(USE_REAL8 ON CACHE BOOL
  "Switch to set HEMCO precision (hp) to 8-byte floating point real"
)
target_compile_definitions(HEMCOBuildProperties INTERFACE
  $<$<BOOL:${USE_REAL8}>:USE_REAL8>
)

#-----------------------------------------------------------------------------
# Add OpenMP compile command to the compilation settings
# Define OMP as a cached variable to prevent overriding existing definition
#-----------------------------------------------------------------------------
set(OMP "TRUE" CACHE STRING "Enables/disables HEMCO built with OpenMP")
if(${OMP})
    find_package(OpenMP REQUIRED)
    ###########################################################################
    # NOTE: Newer versions of CMake (maybe > 3.8) prefer OpenMP::OpenMP
    # rather than ${OpenMP_Fortran_FLAGS} to specify compilation options
    # for OpenMP.  However, this is not supported in older versions.
    # For backwards compatibility, especially with Azure DevOps, we will
    # leave the new syntax commented out.  It can be restored later.
    #  -- Bob Yantosca (28 Jul 2020)
    #
    #    target_compile_options(HEMCOBuildProperties
    #      INTERFACE OpenMP::OpenMP_Fortran
    #    )
    #    target_link_libraries(HEMCOBuildProperties
    #      INTERFACE OpenMP::OpenMP_Fortran
    #    )
    ###########################################################################
    target_compile_options(HEMCOBuildProperties
      INTERFACE ${OpenMP_Fortran_FLAGS}
    )
    target_link_libraries(HEMCOBuildProperties
      INTERFACE ${OpenMP_Fortran_FLAGS}
    )
else()
    target_compile_definitions(HEMCOBuildProperties
        INTERFACE "NO_OMP"
    )
endif()


#-----------------------------------------------------------------------------
# Put all of HEMCO's mod files in build subdir called mod
#-----------------------------------------------------------------------------
set(CMAKE_Fortran_MODULE_DIRECTORY ${HEMCO_BINARY_DIR}/mod)
target_include_directories(HEMCOBuildProperties
    INTERFACE ${HEMCO_BINARY_DIR}/mod
)

#-----------------------------------------------------------------------------
# Append HEMCO's environment variables to CMAKE_PREFIX_PATH
#-----------------------------------------------------------------------------
list(APPEND CMAKE_PREFIX_PATH
    # Possible NetCDF environment variables
    $ENV{NetCDF_F_ROOT}     $ENV{NetCDF_C_ROOT}     $ENV{NetCDF_ROOT}
    $ENV{NETCDF_F_ROOT}     $ENV{NETCDF_C_ROOT}     $ENV{NETCDF_ROOT}
    $ENV{NetCDF_Fortran_ROOT}
    $ENV{NETCDF_FORTRAN_ROOT}

    # Possible HEMCO environmnet variables (still needed? -ewl)
    $ENV{GC_F_BIN} 	    $ENV{GC_BIN}
    $ENV{GC_F_INCLUDE} 	    $ENV{GC_INCLUDE}
    $ENV{GC_F_LIB} 	    $ENV{GC_LIB}
)

#-----------------------------------------------------------------------------
# Link NetCDF-F to HEMCOBuildProperties
#-----------------------------------------------------------------------------
find_package(NetCDF REQUIRED)
target_include_directories(HEMCOBuildProperties
    INTERFACE ${NETCDF_INCLUDE_DIRS}
)
target_link_libraries(HEMCOBuildProperties
    INTERFACE ${NETCDF_LIBRARIES}
)

#-----------------------------------------------------------------------------
# Use the NC_HAS_COMPRESSION def if nf_def_var_deflate is in netcdf.inc
#-----------------------------------------------------------------------------
if(EXISTS ${NETCDF_F77_INCLUDE_DIR}/netcdf.inc)
    file(READ ${NETCDF_F77_INCLUDE_DIR}/netcdf.inc NCINC)
    if("${NCINC}" MATCHES ".*nf_def_var_deflate.*")
        target_compile_definitions(HEMCOBuildProperties
            INTERFACE "NC_HAS_COMPRESSION"
        )
    endif()
endif()

#-----------------------------------------------------------------------------
# Print a description of the source code repo's version
# (Needs more work to enable)
#-----------------------------------------------------------------------------
#get_repo_version(HEMCO_REPO_VERSION ${CMAKE_CURRENT_SOURCE_DIR})
#message(STATUS "HEMCO @ ${HEMCO_REPO_VERSION})

#-----------------------------------------------------------------------------
# For HEMCO Standalone
#
# This conditional block configures the HEMCO build for HEMCO standalone.
# It sets HEMCO_EXE_TARGETS and it configures the HEMCOBuildProperties.
#-----------------------------------------------------------------------------
if(NOT HEMCO_EXTERNAL_CONFIG)

    # Set HEMCO compile flag to standalone
    target_compile_definitions(HEMCOBuildProperties
        INTERFACE "HEMCO_STANDALONE"
    )

    # Specify the HEMCO interface being used
    set(HEMCO_INTERFACE
        standalone
        CACHE STRING "HEMCO interface type (used to choose appropriate I/O option)"
    )

    # Set CMAKE_BUILD_TYPE to Release by default
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release"
            CACHE STRING
            "Set the build type"
            FORCE
    	)
    endif()

    # Set the run directory to one directory up
    set(RUNDIR ".." CACHE PATH "Path to your run directory")
    message(STATUS "Bootstrapping  ${RUNDIR}")

    # Make RUNDIR an absolute path
    get_filename_component(RUNDIR "${RUNDIR}"
        ABSOLUTE BASE_DIR "${CMAKE_BINARY_DIR}"
    )

    # Configure the run directory, if any.
    set(BUILD_WITHOUT_RUNDIR FALSE)
    if(EXISTS ${RUNDIR}/HEMCO_sa_Spec.rc)
        file(STRINGS ${RUNDIR}/HEMCO_sa_Spec.rc HEMCO_SA)
    elseif(NOT EXISTS RUNDIR)
        set(RUNDIR ${CMAKE_BINARY_DIR})
        set(BUILD_WITHOUT_RUNDIR TRUE)
    else()
        message(FATAL_ERROR "Your run directory doesn't
                have HEMCO_sa_Spec.rc! Set RUNDIR to a
                valid run directory.")
    endif()

    # Determine which executables should be built
    set(HEMCO_EXE_TARGETS "hemco_standalone"
        CACHE STRING "Executable targets that get built as a part of \"all\""
    )

    # Add code sanitization options (GNU Fortran only) to HEMCO standalone
    # We need to add these options to both compiler & linker.
    set(SANITIZE OFF CACHE BOOL
        "Switch to turn on code sanitation (i.e. identify memory leaks and similar conditions)"
    )
    if(${SANITIZE})
      if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
	target_compile_options(HEMCOBuildProperties
          INTERFACE "-fsanitize=address"
	)
	target_link_libraries(HEMCOBuildProperties
          INTERFACE "-fsanitize=address"
	)
	target_compile_options(HEMCOBuildProperties
          INTERFACE "-fsanitize=leak"
	)
	target_link_libraries(HEMCOBuildProperties
          INTERFACE "-fsanitize=leak"
	)
	target_compile_options(HEMCOBuildProperties
          INTERFACE "-fsanitize=undefined"
	)
	target_link_libraries(HEMCOBuildProperties
          INTERFACE "-fsanitize=undefined"
	)
      else()
        message( FATAL_ERROR "The SANITIZE option is only defined for GNU Fortran.")
      endif()
    endif()

endif()

#-----------------------------------------------------------------------------
# For HEMCO as submodule in GEOS-Chem "Classic" (GCClassic) superproject
#-----------------------------------------------------------------------------
if(HEMCO_EXTERNAL_CONFIG AND GCCLASSIC_WRAPPER)
    target_compile_definitions(HEMCOBuildProperties
        INTERFACE MODEL_GCCLASSIC
    )

    # Specify the HEMCO interface being used
    set(HEMCO_INTERFACE
        gcclassic
        CACHE STRING "HEMCO interface type (used to choose appropriate I/O option)"
    )
endif()

#-----------------------------------------------------------------------------
# For HEMCO as submodule in MAPL/ESMF superproject
#-----------------------------------------------------------------------------
if(HEMCO_EXTERNAL_CONFIG AND MAPL_ESMF)
    target_link_libraries(HEMCOBuildProperties
        INTERFACE $<LINK_ONLY:MAPL.base>
        INTERFACE $<LINK_ONLY:MAPL.generic>
    )
    target_include_directories(HEMCOBuildProperties 
        INTERFACE $<TARGET_PROPERTY:MAPL.base,INTERFACE_INCLUDE_DIRECTORIES>
        INTERFACE $<TARGET_PROPERTY:MAPL.generic,INTERFACE_INCLUDE_DIRECTORIES>
    )
    target_compile_definitions(HEMCOBuildProperties
        INTERFACE ESMF_ MAPL_ESMF
    )

    # Specify the HEMCO interface being used
    set(HEMCO_INTERFACE
        mapl
        CACHE STRING "HEMCO interface type (used to choose appropriate I/O option)"
    )
endif()

#-----------------------------------------------------------------------------
# Write HEMCOBuildProperties's configuration to a file
#-----------------------------------------------------------------------------
get_target_property(BT_DEFINITIONS  HEMCOBuildProperties
    INTERFACE_COMPILE_DEFINITIONS
)
get_target_property(BT_OPTIONS      HEMCOBuildProperties
    INTERFACE_COMPILE_OPTIONS
)
get_target_property(BT_LIBRARIES    HEMCOBuildProperties
    INTERFACE_LINK_LIBRARIES
)
get_target_property(BT_INCLUDES     HEMCOBuildProperties
    INTERFACE_INCLUDE_DIRECTORIES
)
file(WRITE ${CMAKE_BINARY_DIR}/HEMCOBuildProperties.txt
    "# This file shows the HEMCOBuildProperties's configuration.\n"
    "\n"
    "HEMCOBuildProperties::INTERFACE_COMPILE_DEFINITIONS:${BT_DEFINITIONS}\n"
    "HEMCOBuildProperties::INTERFACE_COMPILE_OPTIONS:${BT_OPTIONS}\n"
    "HEMCOBuildProperties::INTERFACE_LINK_LIBRARIES:${BT_LIBRARIES}\n"
    "HEMCOBuildProperties::INTERFACE_INCLUDE_DIRECTORIES:${BT_INCLUDES}\n"
)

#-----------------------------------------------------------------------------
# Try to compile a simple program that uses NetCDF-Fortran and OpenMP
#-----------------------------------------------------------------------------
if(NOT HEMCO_EXTERNAL_CONFIG AND NOT HEMCO_TRY_RUN_PASSED)

    # Format definitions with -D prefix
    set(TRY_RUN_LIBRARIES ${NETCDF_LIBRARIES})
    set(TRY_RUN_INCLUDE_DIRS ${NETCDF_INCLUDE_DIRS})
    set(TRY_RUN_COMPILER_FLAGS "")
    set(TRY_RUN_DEFINITIONS "")

    if(OMP)
        list(APPEND TRY_RUN_LIBRARIES ${OpenMP_Fortran_LIBRARIES})
        list(APPEND TRY_RUN_COMPILER_FLAGS ${OpenMP_Fortran_FLAGS})
    else()
        list(APPEND TRY_RUN_DEFINITIONS -DNO_OMP)
    endif()

    # Try to compile and run try_compile.F90
    try_run(RUN_FAILED COMPILED_OK
        ${CMAKE_CURRENT_BINARY_DIR}/try_compile                  # binary dir
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeScripts/try_compile.F90 # test src file
        LINK_LIBRARIES ${TRY_RUN_LIBRARIES}
        CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${TRY_RUN_INCLUDE_DIRS}" "-DCMAKE_Fortran_FLAGS=${TRY_RUN_COMPILER_FLAGS}"
        COMPILE_DEFINITIONS ${TRY_RUN_DEFINITIONS}
        COMPILE_OUTPUT_VARIABLE COMPILE_OUTPUT
        RUN_OUTPUT_VARIABLE RUN_OUTPUT
    )

    # Display a warning if its compilation failed
    if(NOT COMPILED_OK)
        if(OMP)
            set(CONDITIONAL_AND_OMP " and OpenMP")
        endif()
        message(WARNING
            "Failed to compile a simple program that uses "
            "NetCDF-Fortran ${CONDITIONAL_AND_OMP}! Could "
            "your NetCDF installation be broken?\nSee "
            "\"FailedCompile.txt\" for more info."
        )
        file(WRITE ${CMAKE_BINARY_DIR}/FailedCompile.txt
            "${COMPILE_OUTPUT}"
        )
    else()
        file(REMOVE ${CMAKE_BINARY_DIR}/FailedCompile.txt)
    endif()

    # Display a warning if its execution failed
    if(RUN_FAILED)
        if(OMP)
            set(CONDITIONAL_AND_OMP "and OpenMP ")
        endif()
        message(WARNING
            "A simple program that uses NetCDF-Fortran "
            "${CONDITIONAL_AND_OMP}compiled successfully, "
            "but its execution failed!\n\nSee "
            "\"FailedExecution.txt\" for more info."
        )
        file(WRITE ${CMAKE_BINARY_DIR}/FailedEasyRun.txt
            "${COMPILE_OUTPUT}\n${RUN_OUTPUT}"
        )
    else()
        file(REMOVE ${CMAKE_BINARY_DIR}/FailedEasyRun.txt
            ${CMAKE_BINARY_DIR}/simple_xy.nc
        )
        set(GC_TRY_RUN_PASSED TRUE CACHE INTERNAL
            "try_run passed" FORCE
        )
    endif()
endif()

#-----------------------------------------------------------------------------
# Print output of environment variable settings
#-----------------------------------------------------------------------------
if(NOT HEMCO_EXTERNAL_CONFIG)
    hco_pretty_print(SECTION "Settings")
    hco_pretty_print(VARIABLE OMP IS_BOOLEAN)
    hco_pretty_print(VARIABLE USE_REAL8 IS_BOOLEAN)
    hco_pretty_print(VARIABLE SANITIZE IS_BOOLEAN)
endif()

#-----------------------------------------------------------------------------
# Add source code directory
#-----------------------------------------------------------------------------
add_subdirectory(src)
