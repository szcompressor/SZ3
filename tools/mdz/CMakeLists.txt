add_library(MDZ INTERFACE)
target_include_directories(
        MDZ INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(MDZ INTERFACE SZ3 ${ZSTD_dep})

install(TARGETS MDZ EXPORT MDZTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MDZ)
export(TARGETS MDZ FILE MDZ.cmake)

set(exes mdz)
foreach (EXE IN LISTS exes)
    message(STATUS "Target = ${EXE}")
    add_executable(${EXE} ${EXE}.cpp)
    target_include_directories(${EXE} PRIVATE ${PROJECT_SOURCE_DIR}/include $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
    target_link_libraries(${EXE} MDZ)
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${EXE} TYPE BIN)
endforeach ()