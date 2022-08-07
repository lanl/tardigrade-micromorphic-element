# Project

# Build the micromorphic balance equations
add_library (${BALANCE_EQUATION_LIBRARY} SHARED "${BALANCE_EQUATION_LIBRARY_FILENAME}.cpp" "${BALANCE_EQUATION_LIBRARY_FILENAME}.h")
set_target_properties(${BALANCE_EQUATION_LIBRARY} PROPERTIES CXX_STANDARD 11
                                                 PUBLIC_HEADER ${BALANCE_EQUATION_LIBRARY_FILENAME}.h)
target_link_libraries(${BALANCE_EQUATION_LIBRARY} error_tools)
target_compile_options(${BALANCE_EQUATION_LIBRARY} PUBLIC)

# Local builds of upstream projects require local include paths
if(NOT cmake_build_type_lower STREQUAL "release")
    target_include_directories(${BALANCE_EQUATION_LIBRARY} PUBLIC
                               "${vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
                               "${error_tools_SOURCE_DIR}/${CPP_SRC_PATH}")
endif()

foreach(package ${BALANCE_EQUATION_LIBRARY} ${MATERIAL_MODEL_LIBRARY})
    install(TARGETS ${package}
            EXPORT ${package}_Targets
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endforeach(package)
