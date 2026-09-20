install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MDZ)

set(exes mdz
        mdz_smoke_test)
foreach (EXE IN LISTS exes)
    message(STATUS "Target = ${EXE}")
    add_executable(${EXE} ${EXE}.cpp)
    target_include_directories(${EXE} PUBLIC
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/MDZ>"
            )
    target_link_libraries(${EXE} SZ3)
    # As in tools/sz3: the spelled-out path misses both the per-config subdirectory and the .exe.
    install(PROGRAMS $<TARGET_FILE:${EXE}> TYPE BIN)
endforeach ()

