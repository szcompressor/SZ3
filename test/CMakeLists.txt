file(GLOB cpu_source_files "${CMAKE_SOURCE_DIR}/test/*.cpp")


foreach (SRCFILE IN LISTS cpu_source_files)
    get_filename_component(EXE ${SRCFILE} NAME_WLE)
    message(STATUS "Target = ${EXE}")
    add_executable(${EXE} ${EXE}.cpp)

    set(deps "")
    list(APPEND deps SZ3)
    list(APPEND deps ${ZSTD_dep})
    if (GSL_FOUND)
        list(APPEND deps GSL::gsl)
    endif ()
    if (OpenMP_CXX_FOUND)
        list(APPEND deps OpenMP::OpenMP_CXX)
    endif ()
    target_link_libraries(${EXE} ${deps})

    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${EXE} TYPE BIN)
endforeach ()
install(FILES ${PROJECT_SOURCE_DIR}/test/testfloat_8_8_128.dat DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})
