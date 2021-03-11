add_library(SZ3_EXAALT
        dummy.cpp
        )

target_link_libraries(SZ3_EXAALT ${ZSTD_dep})
target_include_directories(SZ3_EXAALT PUBLIC
        ${PROJECT_SOURCE_DIR}/exaalt/include
        ${PROJECT_SOURCE_DIR}/include
        ${ZSTD_INCLUDES})

install(TARGETS SZ3_EXAALT
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/lib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/lib
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/include
        )
