add_library(
        sz3c
        src/sz3c.cpp
)

target_link_libraries(
        sz3c
        PUBLIC SZ3
)

target_include_directories(
        sz3c
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:sz3c>
)

install(TARGETS sz3c EXPORT SZ3C
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        )
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/sz3c)
export(TARGETS sz3c FILE sz3c.cmake)
