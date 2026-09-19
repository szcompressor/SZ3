"""
Setup script for pysz - Python bindings for SZ3
Automatically downloads and builds SZ3 with bundled zstd.
"""

import sys
import shutil
import subprocess
from pathlib import Path
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from Cython.Build import cythonize
import numpy as np



# A released tag, never a branch or a bare commit: a published wheel has to be buildable from a
# published source, and a branch is gone once it merges. pysz is tagged on its own schedule, so
# this moves when SZ3 releases, not when SZ3 changes.
SZ3_VERSION = "3.3.2"

# Which of those releases this is decides where the bundled Zstd lives and what it is called:
# v3.3.2 fetches it into build/_deps and builds `libzstd`, later trees vendor it under tools/zstd
# and build `libsz3_zstd`. Both are found in the tree that was just built, so bumping the tag
# above is the only edit the next release needs.
ZSTD_HEADER_DIRS = (("tools", "zstd", "lib"), ("build", "_deps", "zstdfetched-src", "lib"))
ZSTD_LIBRARY_DIRS = (("build", "tools", "zstd"),
                     ("build", "tools", "zstd", "Release"),
                     ("build", "tools", "zstd", "Debug"))
ZSTD_LIBRARY_NAMES = ("sz3_zstd", "zstd")


def find_zstd_header_dir(sz3_dir):
    for parts in ZSTD_HEADER_DIRS:
        candidate = sz3_dir.joinpath(*parts)
        if (candidate / "zstd.h").is_file():
            return candidate
    raise RuntimeError(f"no bundled zstd.h under {sz3_dir}; is SZ3_VERSION a tag that bundles Zstd?")


def find_zstd_library(sz3_dir):
    """The link name of the bundled Zstd this SZ3 build produced, and the directory holding it."""
    for name in ZSTD_LIBRARY_NAMES:
        for parts in ZSTD_LIBRARY_DIRS:
            directory = sz3_dir.joinpath(*parts)
            for pattern in (f"lib{name}.*", f"{name}.lib"):
                if any(directory.glob(pattern)):
                    return name, directory
    raise RuntimeError(f"no bundled Zstd library under {sz3_dir / 'build' / 'tools' / 'zstd'}")


class BuildSZ3Extension(_build_ext):

    def run(self):
        sz3_dir = self.download_and_build_sz3()
        zstd_name, zstd_dir = find_zstd_library(sz3_dir)
        print(f"Linking bundled Zstd: {zstd_name} from {zstd_dir}")

        for ext in self.extensions:
            ext.include_dirs.insert(0, str(sz3_dir / "include"))
            ext.include_dirs.insert(0, str(sz3_dir / "build" / "include"))
            ext.include_dirs.append(str(find_zstd_header_dir(sz3_dir)))
            ext.libraries.append(zstd_name)
            ext.library_dirs.append(str(zstd_dir))

        super().run()

        # A shared bundled Zstd has to ride along in the wheel next to the extension that loads it
        # -- the rpath below points there. A static one is already inside the extension.
        package_dir = Path(self.build_lib) / "pysz"
        if package_dir.exists():
            for pattern in (f"lib{zstd_name}.dylib", f"lib{zstd_name}.so", f"{zstd_name}.dll"):
                for shared in sorted(zstd_dir.glob(pattern)):
                    shutil.copy2(shared, package_dir / shared.name)
                    print(f"Copied {shared.name} to package")
                    return

    def download_and_build_sz3(self):
        build_temp = Path(self.build_temp).absolute()
        build_temp.mkdir(parents=True, exist_ok=True)
        sz3_dir = build_temp / "SZ3"
        
        if (sz3_dir / "build" / "include" / "SZ3" / "version.hpp").exists():
            print(f"SZ3 already built at: {sz3_dir}")
            return sz3_dir
        
        if not sz3_dir.exists():
            print(f"Cloning SZ3 v{SZ3_VERSION}...")
            subprocess.run([
                "git", "clone", "--depth", "1",
                "--branch", f"v{SZ3_VERSION}",
                "--single-branch",
                "https://github.com/szcompressor/SZ3.git",
                str(sz3_dir)
            ], check=True)

        build_dir = sz3_dir / "build"
        build_dir.mkdir(exist_ok=True)
        
        cmake_args = ["cmake"]
        cmake_args.extend([
            "-DCMAKE_BUILD_TYPE=Release",
            "-DBUILD_TESTING=OFF",
            "-DBUILD_SZ3_BINARY=OFF",
            "-DSZ3_USE_BUNDLED_ZSTD=ON",
            ".."
        ])
        subprocess.run(cmake_args, cwd=build_dir, check=True)
        subprocess.run(["cmake", "--build", ".", "-j"], cwd=build_dir, check=True)
        print(f"Built SZ3 v{SZ3_VERSION}")
        return sz3_dir




def create_extensions():
    include_dirs = [np.get_include()]
    library_dirs = []
    libraries = []  # the bundled Zstd is appended once the build knows what it is called
    extra_compile_args = []
    extra_link_args = []
    
    if sys.platform == 'win32':
        extra_compile_args.extend(['/std:c++17', '/O2'])
    elif sys.platform == 'darwin':
        extra_compile_args.extend(['-std=c++17', '-O3', '-stdlib=libc++', '-mmacosx-version-min=10.9'])
        extra_link_args.extend(['-stdlib=libc++', '-Wl,-rpath,@loader_path'])
    elif sys.platform == 'linux':
        extra_compile_args.extend(['-std=c++17', '-O3'])
        extra_link_args.extend(['-Wl,-rpath,$ORIGIN'])
    
    extensions = [
        Extension(
            "pysz.sz",
            sources=["src/pysz/sz.pyx"],
            include_dirs=include_dirs,
            libraries=libraries,
            library_dirs=library_dirs,
            language='c++',
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args,
        ),
    ]
    
    return cythonize(extensions, compiler_directives={'language_level': '3', 'embedsignature': True})


if __name__ == "__main__":
    setup(
        name="pysz",
        version="1.0.3",
        packages=["pysz"],
        package_dir={"": "src"},
        ext_modules=create_extensions(),
        cmdclass={'build_ext': BuildSZ3Extension},
        test_suite="tests",
        tests_require=["pytest"],
    )
