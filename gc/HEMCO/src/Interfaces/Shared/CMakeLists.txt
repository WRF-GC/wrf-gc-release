# HEMCO/src/Interfaces/Shared/CMakeLists.txt

add_library(HCOI_Shared STATIC EXCLUDE_FROM_ALL
        hco_interface_common.F90
)
target_link_libraries(HCOI_Shared
	PUBLIC HCOX
)
target_include_directories(HCOI_Shared
	INTERFACE ${HEMCO_BINARY_DIR}/mod
)
target_compile_options(HCOI_Shared
	PRIVATE
		""
		$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-ffree-line-length-none>
)
