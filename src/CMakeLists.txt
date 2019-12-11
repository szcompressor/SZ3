configure_file(${PROJECT_SOURCE_DIR}/include/def.hpp.in ${PROJECT_SOURCE_DIR}/include/def.hpp)

add_library(SZ3 SHARED
        dummy.cpp
        )
target_include_directories(SZ3 PUBLIC ${PROJECT_SOURCE_DIR}/include)
# install(FILES ${header_files} DESTINATION ${PROJECT_BINARY_DIR}/include)

install(TARGETS SZ3
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/lib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/lib
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/include
        )
