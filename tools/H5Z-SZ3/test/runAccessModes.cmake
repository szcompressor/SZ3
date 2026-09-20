# The suite tests an install tree, not a build tree.
file(REMOVE_RECURSE "${SZ3_PREFIX}")
execute_process(
        COMMAND ${CMAKE_COMMAND} --install "${SZ3_BUILD_DIR}" --prefix "${SZ3_PREFIX}"
        RESULT_VARIABLE install_result OUTPUT_QUIET)
if (NOT install_result EQUAL 0)
    message(FATAL_ERROR "could not install into ${SZ3_PREFIX}")
endif ()
execute_process(COMMAND bash "${SZ3_SCRIPT}" "${SZ3_PREFIX}" "${SZ3_H5BIN}" RESULT_VARIABLE suite_result)
if (NOT suite_result EQUAL 0)
    message(FATAL_ERROR "filterAccessModes.sh reported failures")
endif ()
