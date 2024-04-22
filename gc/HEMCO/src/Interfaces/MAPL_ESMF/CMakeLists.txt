# HEMCO/src/Interfaces/MAPL_ESMF/CMakeLists.txt

add_library(HCOI_MAPL_ESMF STATIC EXCLUDE_FROM_ALL
        hcoi_esmf_mod.F90
)
target_link_libraries(HCOI_MAPL_ESMF
    PUBLIC HCOI_Shared
)
target_include_directories(HCOI_MAPL_ESMF
	INTERFACE ${HEMCO_BINARY_DIR}/mod
)
target_compile_options(HCOI_MAPL_ESMF
	PRIVATE
		""
		$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-ffree-line-length-none>
)
