add_library (SZ3 SHARED 
	dummy.cpp
	# encoder/Encoder.cpp
	# predictor/Predictor.cpp
	# quantizer/Quantizer.cpp
	# compressor/BlockCompressor.cpp
	# compressor/Compressor.cpp
	# utils/BlockIterator.cpp
	)
# target_include_directories(SZ3 PRIVATE ${ZSTD_INCLUDES})
target_include_directories(SZ3 PRIVATE ${PROJECT_SOURCE_DIR}/include)
# target_link_libraries(SZ3 ${ZSTD_LIBS})
# file(GLOB header_files ${PROJECT_SOURCE_DIR}/include/*.hpp)
# install(FILES ${header_files} DESTINATION ${PROJECT_BINARY_DIR}/include)

install (TARGETS SZ3
  LIBRARY DESTINATION ${PROJECT_BINARY_DIR}/lib
  ARCHIVE DESTINATION ${PROJECT_BINARY_DIR}/lib
  INCLUDES DESTINATION ${PROJECT_BINARY_DIR}/include
)
