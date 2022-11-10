set(exes mdz)

foreach (EXE IN LISTS exes)
    message(STATUS "Target = ${EXE}")
    add_executable(${EXE} ${EXE}.cpp)
    target_include_directories(${EXE} PRIVATE ${PROJECT_SOURCE_DIR}/include $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
    target_link_libraries(${EXE} SZ3 ${ZSTD_dep})
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${EXE} TYPE BIN)
endforeach ()
