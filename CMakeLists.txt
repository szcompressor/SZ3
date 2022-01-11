cmake_minimum_required(VERSION 3.13)
project(SZ3 VERSION 3.0.0.0)
include(GNUInstallDirs)
include(CTest)

find_package(PkgConfig)
find_package(OpenMP)

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(
        ${PROJECT_NAME} INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        )
target_compile_features(${PROJECT_NAME}
    INTERFACE cxx_std_17
  )

target_link_libraries(${PROJECT_NAME} INTERFACE OpenMP::OpenMP_CXX)

pkg_search_module(ZSTD IMPORTED_TARGET libzstd)
if(ZSTD_FOUND)
  target_link_libraries(${PROJECT_NAME} INTERFACE PkgConfig::ZSTD)
else()
  add_subdirectory(zstd)
  target_link_libraries(${PROJECT_NAME} INTERFACE zstd)
endif()

find_package(GSL)
if (GSL_FOUND)
  target_compile_definitions(${PROJECT_NAME} INTERFACE -DENABLE_GSL="TRUE")
endif()

add_subdirectory(test)

install(TARGETS ${PROJECT_NAME} 
  EXPORT SZ3Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT SZ3Targets NAMESPACE SZ3:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZ3)
include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/SZ3Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/SZ3Config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZ3
  )
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/SZ3ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
  )
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/SZ3Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/SZ3ConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZ3
  )
