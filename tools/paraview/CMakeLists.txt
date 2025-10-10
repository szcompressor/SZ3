find_package(ParaView REQUIRED)

# Build the plugin by default (for ParaView v6.0).
set(PARAVIEW_PLUGIN_ENABLE_SZ3Reader ON)

# Scan the plugin.
paraview_plugin_scan(
  PLUGIN_FILES "${CMAKE_CURRENT_SOURCE_DIR}/SZ3Reader/paraview.plugin"
  PROVIDES_PLUGINS plugins
  REQUIRES_MODULES modules)

# Ensure we have required modules.
foreach (module IN LISTS REQUIRES_MODULES)
  if (NOT TARGET "${module}")
    message(FATAL_ERROR
      "The required ${module} module is not available in the "
      "ParaView build.")
  endif ()
endforeach ()

# Build the plugin.
paraview_plugin_build(
  RUNTIME_DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY_DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY_SUBDIRECTORY "${PARAVIEW_PLUGIN_SUBDIR}"
  TARGET paraview_plugins
  PLUGINS ${plugins})