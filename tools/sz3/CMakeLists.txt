file(GLOB cpu_source_files "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

foreach (SRCFILE IN LISTS cpu_source_files)
    get_filename_component(EXE ${SRCFILE} NAME_WLE)
    message(STATUS "Target = ${EXE}")
    add_executable(${EXE} ${EXE}.cpp)
    target_link_libraries(${EXE} PRIVATE SZ3)
    # Ask for the binary rather than spelling out where it lands: a multi-config generator puts it
    # under a per-config subdirectory, and the install then cannot find it.
    if (SZ3_INSTALL AND (EXE STREQUAL "sz3" OR EXE STREQUAL "sz3_smoke_test"))
        install(PROGRAMS $<TARGET_FILE:${EXE}> TYPE BIN)
    endif ()

endforeach ()

#add_subdirectory(demo)
