set(exes sz_demo
        sz_fast
        sz_truncate
        sz_pw
        test_huffman_coding
        test_arithmetic_coding
        sz_interp
        sz_progressive_interp
        sz_progressive_interp_v2
        sz_progressive_interp_v3
        sz_progressive_interp_v4
        sz_rtm
        sz_interp_double
        )
foreach (EXE IN LISTS exes)
    message(STATUS "Target = ${EXE}")
    add_executable(${EXE} ${EXE}.cpp)
    target_link_libraries(${EXE} SZ3)
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${EXE} TYPE BIN)
endforeach ()
install(FILES ${PROJECT_SOURCE_DIR}/test/testfloat_8_8_128.dat DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})
