set(src_HCOCore
	hco_arr_mod.F90
	hco_calc_mod.F90
	hco_chartools_mod.F90
	hco_clock_mod.F90
	hco_config_mod.F90
	hco_datacont_mod.F90
	hco_diagn_mod.F90
	hco_driver_mod.F90
	hco_emislist_mod.F90
	hco_error_mod.F90
	hco_extlist_mod.F90
	hco_filedata_mod.F90
	hco_fluxarr_mod.F90
	hco_geotools_mod.F90
	hco_interp_mod.F90
	hco_logfile_mod.F90
	hco_readlist_mod.F90
	hco_restart_mod.F90
	hco_scale_mod.F90
	hco_state_mod.F90
	hco_tidx_mod.F90
	hco_timeshift_mod.F90
	hco_types_mod.F90
	hco_unit_mod.F90
	hco_vertgrid_mod.F90
	hcoio_dataread_mod.F90
	hcoio_diagn_mod.F90
	hcoio_messy_mod.F90
	hcoio_util_mod.F90
	interpreter.F90
	messy_ncregrid_base.F90
	messy_ncregrid_mpi.F90
)

# Choose the appropriate HEMCO I/O library for the given interface.
# Currently available: standard ("HEMCO default") and mapl (MAPL)
#
# To add a future option:
# elseif(HEMCO_INTERFACE STREQUAL future)
# 	list(APPEND src_HCOCore hcoio_read_future_mod.F90 hcoio_write_future_mod.F90)
#
# Note that models that do not use CMake, i.e. WRF-GC (std); CESM2-GC (std) will
# completely bypass CMake and include their own build option, hence all the HCOIO_
# modules still have to include C-Preprocessor variables.
#
# For discussion, refer to: https://github.com/geoschem/HEMCO/issues/87
if(HEMCO_INTERFACE STREQUAL standalone)
	list(APPEND src_HCOCore hcoio_read_std_mod.F90 hcoio_write_std_mod.F90)
elseif(HEMCO_INTERFACE STREQUAL gcclassic)
	list(APPEND src_HCOCore hcoio_read_std_mod.F90 hcoio_write_std_mod.F90)
elseif(HEMCO_INTERFACE STREQUAL mapl)
	list(APPEND src_HCOCore hcoio_read_mapl_mod.F90 hcoio_write_mapl_mod.F90)
endif()

add_library(HCO STATIC EXCLUDE_FROM_ALL
	${src_HCOCore}
)

target_link_libraries(HCO
	PUBLIC GeosUtilHco NcdfUtilHco
)
target_compile_options(HCO
	PRIVATE
		""
		$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-ffree-line-length-none>
)
target_include_directories(HCO
	INTERFACE ${HEMCO_BINARY_DIR}/mod
)
