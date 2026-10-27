file(MAKE_DIRECTORY "${WORK}")

function(run_sz3)
    execute_process(COMMAND "${SZ3_EXECUTABLE}" ${ARGN}
                    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
    if (NOT status STREQUAL "0")
        message(FATAL_ERROR "sz3 ${ARGN} failed (${status}): ${output}${error}")
    endif ()
endfunction()

# Make a real stream and a reference reconstruction using explicit dimensions.
run_sz3(-f -i "${INPUT}" -3 128 8 8 -M ABS 0.001
        -z "${WORK}/input.sz3" -o "${WORK}/reference.dat")

# Both CLI syntaxes must infer dimensions from the stream, with and without
# the original input supplied for error reporting (-a).
foreach(mode IN ITEMS modern modern_analysis sz2 sz2_analysis)
    set(args -f)
    if (mode MATCHES "^modern")
        list(APPEND args -z "${WORK}/input.sz3" -o "${WORK}/${mode}.dat")
    else()
        list(APPEND args -s "${WORK}/input.sz3" -x "${WORK}/${mode}.dat")
    endif()
    if (mode MATCHES "analysis$")
        list(APPEND args -i "${INPUT}" -a)
    endif()
    run_sz3(${args})
    execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
                    "${WORK}/reference.dat" "${WORK}/${mode}.dat"
                    RESULT_VARIABLE different)
    if (NOT different STREQUAL "0")
        message(FATAL_ERROR "${mode} reconstruction differs from explicit-dimension decoding")
    endif()
endforeach()

# Missing compression dimensions and explicit zero dimensions remain errors.
foreach(mode IN ITEMS missing zero)
    set(args -f -i "${INPUT}" -M ABS 0.001 -z "${WORK}/invalid.sz3")
    if (mode STREQUAL "zero")
        list(APPEND args -1 0)
    endif()
    execute_process(COMMAND "${SZ3_EXECUTABLE}" ${args}
                    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
    if (status STREQUAL "0" OR NOT error MATCHES "dimension")
        message(FATAL_ERROR "${mode} compression dimensions were not rejected: ${output}${error}")
    endif()
endforeach()
