add_executable(sz_huffman_test sz_huffman_test.cpp)
target_include_directories(sz_huffman_test PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(sz_huffman_test SZ3 ${ZSTD_LIBS})
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/sz_demo TYPE BIN)

add_executable(sz_demo sz_demo.cpp)
target_include_directories(sz_demo PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(sz_demo SZ3 ${ZSTD_LIBS})
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/sz_huffman_test TYPE BIN)


add_executable(test_arithmetic test_arithmetic_coding.cpp)
target_include_directories(test_arithmetic PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(test_arithmetic SZ3 ${ZSTD_LIBS})
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/test_arithmetic TYPE BIN)


install(FILES ${PROJECT_SOURCE_DIR}/test/testfloat_8_8_128.dat  DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})
