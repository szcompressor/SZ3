"""
Setup script for pysz - Python bindings for SZ3

This package downloads and builds SZ3 automatically.
SZ3 is downloaded from GitHub and built with its own CMake configuration,
which automatically handles zstd dependency.
"""

import sys
import os
import subprocess
import shutil
from pathlib import Path
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from Cython.Build import cythonize
import numpy as np



def get_sz3_version():
    """Get SZ3 version from GitHub."""
    return "3.3.1"

class BuildSZ3Extension(_build_ext):
    """Custom build command that downloads and builds SZ3 first."""
    
    def run(self):
        # Download and build SZ3
        sz3_dir = self.download_and_build_sz3()
        
        # Update include directories with SZ3 paths
        for ext in self.extensions:
            # Add SZ3 include directories
            ext.include_dirs.insert(0, str(sz3_dir / "include"))
            ext.include_dirs.insert(0, str(sz3_dir / "build" / "include"))
            
            # Add bundled zstd from SZ3's build
            zstd_build_lib = sz3_dir / "build" / "tools" / "zstd"
            zstd_include = sz3_dir / "build" / "_deps" / "zstdfetched-src" / "lib"
            if zstd_include.exists():
                ext.include_dirs.append(str(zstd_include))
            if zstd_build_lib.exists():
                ext.library_dirs.append(str(zstd_build_lib))
        
        # Continue with normal build
        super().run()
        
        # Copy zstd library to package directory after build
        zstd_lib = sz3_dir / "build" / "tools" / "zstd" / ("libzstd.dylib" if sys.platform == "darwin" else "libzstd.so")
        if zstd_lib.exists():
            package_dir = Path(self.build_lib) / "pysz"
            if package_dir.exists():
                dest_lib = package_dir / zstd_lib.name
                shutil.copy2(zstd_lib, dest_lib)
                print(f"✓ Copied {zstd_lib.name} to package")
                
                # Fix rpath on macOS to load from package directory
                if sys.platform == "darwin":
                    for so_file in package_dir.glob("*.so"):
                        subprocess.run([
                            "install_name_tool", "-change",
                            "@rpath/libzstd.dylib",
                            "@loader_path/libzstd.dylib",
                            str(so_file)
                        ], check=False)  # Don't fail if already correct
                    print(f"✓ Fixed rpath for .so files")
    
    def download_and_build_sz3(self):
        """Download SZ3 from GitHub and build it."""
        version = get_sz3_version()
        build_temp = Path(self.build_temp).absolute()
        build_temp.mkdir(parents=True, exist_ok=True)
        
        sz3_dir = build_temp / "SZ3"
        
        # Check if already downloaded and built
        if (sz3_dir / "build" / "include" / "SZ3" / "version.hpp").exists():
            print(f"SZ3 already built at: {sz3_dir}")
            return sz3_dir
        
        print("=" * 70)
        print(f"Downloading and building SZ3 v{version}...")
        print("=" * 70)
        
        # Clone SZ3 if not exists
        if not sz3_dir.exists():
            print(f"\n[1/3] Cloning SZ3 v{version} from GitHub...")
            subprocess.run([
                "git", "clone", "--depth", "1",
                "--branch", f"v{version}",
                "--single-branch",
                "https://github.com/szcompressor/SZ3.git",
                str(sz3_dir)
            ], check=True)
            print(f"✓ Clone complete (tag: v{version})")
        
        # Build SZ3
        print(f"\n[2/3] Building SZ3...")
        build_dir = sz3_dir / "build"
        build_dir.mkdir(exist_ok=True)
        
        # Configure with CMake - force bundled zstd for portability
        subprocess.run([
            "cmake",
            "-DCMAKE_BUILD_TYPE=Release",
            "-DBUILD_TESTING=OFF",
            "-DSZ3_USE_BUNDLED_ZSTD=ON",
            ".."
        ], cwd=build_dir, check=True)
        
        # Build
        subprocess.run([
            "cmake", "--build", ".", "-j"
        ], cwd=build_dir, check=True)
        
        print("✓ SZ3 build complete")
        print(f"\n[3/3] Using SZ3 v{version} from: {sz3_dir}")
        print("=" * 70 + "\n")
        
        return sz3_dir




def create_extensions():
    """Create Cython extension modules."""
    
    # Basic setup - paths will be updated during build
    include_dirs = [np.get_include()]
    library_dirs = []
    libraries = ['zstd']
    
    # Compiler flags
    extra_compile_args = ['-std=c++17', '-O3']
    extra_link_args = []
    
    if sys.platform == 'darwin':
        extra_compile_args.extend(['-stdlib=libc++', '-mmacosx-version-min=10.9'])
        extra_link_args.extend(['-stdlib=libc++'])
    elif sys.platform == 'linux':
        extra_link_args.extend(['-Wl,-rpath,$ORIGIN'])
    
    # Create extensions
    extensions = [
        Extension(
            "pysz.pyConfig",
            sources=["src/pysz/pyConfig.pyx"],
            include_dirs=include_dirs,
            libraries=libraries,
            library_dirs=library_dirs,
            language='c++',
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args,
        ),
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
    
    return cythonize(
        extensions,
        compiler_directives={
            'language_level': '3',
            'embedsignature': True,
        }
    )


# Run setup
if __name__ == "__main__":
    sz3_ver = get_sz3_version()
    setup(
        name="pysz",
        version=f"1.0.0+sz3.{sz3_ver}",  # Version includes SZ3 version
        packages=["pysz"],
        package_dir={"": "src"},
        ext_modules=create_extensions(),
        cmdclass={'build_ext': BuildSZ3Extension},
    )
