add_executable(test_pattern_compressor test_pattern_compressor.cpp)
target_include_directories(test_pattern_compressor PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(test_pattern_compressor SZ3_EXAALT)

add_executable(test_persist test_persist.cpp)
target_include_directories(test_persist PRIVATE ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/exaalt/include)
target_link_libraries(test_persist SZ3_EXAALT)

add_executable(test_instant test_instant.cpp)
target_include_directories(test_instant PRIVATE ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/exaalt/include)
target_link_libraries(test_instant SZ3_EXAALT)

add_executable(test_multilevel test_multilevel.cpp)
target_include_directories(test_multilevel PRIVATE ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/exaalt/include)
target_link_libraries(test_multilevel SZ3_EXAALT)
