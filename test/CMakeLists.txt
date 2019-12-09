add_executable (sz_test sz_compress_test.cpp)
add_executable (sz_test2 sz_compress_test2.cpp)

target_include_directories (sz_test PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories (sz_test2 PRIVATE ${ZSTD_INCLUDES} ${PROJECT_SOURCE_DIR}/include)

target_link_libraries (sz_test SZ3)
target_link_libraries (sz_test2 SZ3 ${ZSTD_LIBS})

add_executable (sz_huffman_test sz_huffman_test.cpp)
# target_include_directories (sz_test PRIVATE ${ZSTD_INCLUDES})
target_include_directories (sz_huffman_test PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries (sz_huffman_test SZ3)
