# HEMCO/src/Interfaces/Standalone/CMakeLists.txt

add_library(HCOI_Standalone STATIC EXCLUDE_FROM_ALL
	hcoi_standalone_mod.F90
)
target_link_libraries(HCOI_Standalone
	PUBLIC HCOI_Shared
)
target_include_directories(HCOI_Standalone
	INTERFACE ${HEMCO_BINARY_DIR}/mod
)
target_compile_options(HCOI_Standalone
	PRIVATE
		""
		$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-ffree-line-length-none>
)

if("hemco_standalone" IN_LIST HEMCO_EXE_TARGETS)
	add_executable(hemco_standalone
		hemco_standalone.F90
	)
	target_link_libraries(hemco_standalone
		PUBLIC HCOI_Standalone
	)
	set_target_properties(hemco_standalone PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
	install(TARGETS hemco_standalone
		RUNTIME DESTINATION ${RUNDIR}
	)
endif()
