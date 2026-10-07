add_library(
        SZ3c
        src/sz3c.cpp
)

target_link_libraries(
        SZ3c
        PRIVATE SZ3
)

target_compile_definitions(SZ3c PRIVATE SZ3C_EXPORTS)
if (NOT BUILD_SHARED_LIBS)
    # Otherwise sz3c.h marks the functions of a static SZ3c dllimport on Windows.
    target_compile_definitions(SZ3c PUBLIC SZ3C_STATIC)
endif ()
set_target_properties(SZ3c PROPERTIES CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
set_target_properties(SZ3c PROPERTIES VERSION ${SZ3_VERSION} SOVERSION ${SZ3_VERSION_MAJOR}.${SZ3_VERSION_MINOR})

target_include_directories(
        SZ3c
        PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SZ3c>"
)

if (SZ3_INSTALL)
    install(TARGETS SZ3c EXPORT SZ3c
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            )
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SZ3c)
    install(EXPORT SZ3c NAMESPACE SZ3:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZ3)
endif ()
