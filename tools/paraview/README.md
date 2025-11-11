# SZ3Reader Plugin
ParaView/VTK reader for direct visualization of data files compressed using the SZ3 compressor.

Version: 0.1

General information
--------------------
For more information on SZ3, see <https://szcompressor.org/>.

For more information on ParaView, see <https://www.paraview.org/>.

The SZ3Reader plugin is built on the SZ3 C++ API. It reads SZ3-compressed files (e.g., `.sz3`).

The current version supports regular data with up to three dimensions in single- or double-precision floating-point format.


Installation
-------------

Before the plugin is included in the official ParaView release, you will need to manually build the plugin from sources and load it in ParaView.

#### Prerequisites

- Download the ParaView sources from <https://www.paraview.org/download/>. If you have ParaView already installed, please make sure you select the same version from the drop-down list.
- Follow the [ParaView building process](https://www.paraview.org/paraview-docs/latest/cxx/md__builds_gitlab-kitware-sciviz-ci_Documentation_dev_build.html).

#### Building the Plugin

1. Configure the SZ3 project with CMake and enable the ParaView plugin:

   ```bash
   cmake -DBUILD_PARAVIEW_PLUGIN=ON [other options] ..
   ```

2. If CMake does not find ParaView automatically, set `ParaView_DIR` as `[-DParaView_DIR=/path/to/paraview/lib/cmake/paraview-<version>]`. The option `PARAVIEW_PLUGIN_ENABLE_SZ3Reader` will be set to `ON` by default.

3. Build the project (optional: use parallel jobs with `-j`):

   ```bash
   cmake --build . -j [num_threads]
   ```

#### Loading the Plugin in ParaView

1. After the build completes, open ParaView.
2. Navigate to `Tools -> Manage Plugins -> Load New`.
3. Browse to and select the plugin file:
   - **Linux/macOS**: `{SZ3_BUILD_DIR}/lib/paraview-{VERSION}/plugins/SZ3Reader/SZ3Reader.so`
   - **Windows**: `{SZ3_BUILD_DIR}/lib/paraview-{VERSION}/plugins/SZ3Reader/SZ3Reader.dll`
4. Check `Auto Load` if you want the plugin to be loaded automatically every time you run ParaView.


Contact Info
-------------

This plugin can certainly be improved and any suggestions or contributions are welcome.

Please report bugs, problems, and suggestions to:

- Guoxi Liu (email: [liu.12722@osu.edu](mailto:liu.12722@osu.edu))
- Siheng Zhang (email: [zhang.16458@buckeyemail.osu.edu](mailto:zhang.16458@buckeyemail.osu.edu))
