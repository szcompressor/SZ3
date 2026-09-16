add_executable(module_bench module_bench.cpp)
target_link_libraries(module_bench PRIVATE SZ3)
target_compile_features(module_bench PUBLIC cxx_std_17)
install(TARGETS module_bench EXPORT SZ3Targets)
