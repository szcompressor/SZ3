file(GLOB cpu_source_files "${CMAKE_SOURCE_DIR}/tools/sz3/*.cpp")

foreach (SRCFILE IN LISTS cpu_source_files)
    get_filename_component(EXE ${SRCFILE} NAME_WLE)
    message(STATUS "Target = ${EXE}")
    add_executable(${EXE} ${EXE}.cpp)
    target_link_libraries(${EXE} PRIVATE SZ3)
    # Ask for the binary rather than spelling out where it lands: a multi-config generator puts it
    # under a per-config subdirectory, and the install then cannot find it.
    install(PROGRAMS $<TARGET_FILE:${EXE}> TYPE BIN)

endforeach ()
install(FILES testfloat_8_8_128.dat DESTINATION ${CMAKE_INSTALL_DATADIR}/SZ3)

# Demo executables
file(GLOB demo_source_files "${CMAKE_SOURCE_DIR}/tools/sz3/demo/*.cpp")
foreach (SRCFILE IN LISTS demo_source_files)
    get_filename_component(EXE ${SRCFILE} NAME_WLE)
    message(STATUS "Target = ${EXE}")
    add_executable(${EXE} ${SRCFILE})
    target_link_libraries(${EXE} PRIVATE SZ3)
    install(PROGRAMS $<TARGET_FILE:${EXE}> TYPE BIN)
endforeach ()
