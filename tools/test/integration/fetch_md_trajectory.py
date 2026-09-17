#!/usr/bin/env python3
"""Turn a published MD trajectory into the raw float32 arrays the integration tests read.

The SDRBench fields the suite already covers are 1D and 2D, so the trajectory layout
ALGO_BIOMDXTC was written for -- {frames, atoms, xyz} -- is never exercised on molecular
dynamics data. These are the MDAnalysisData benchmark trajectories, downloaded from the
figshare files that package points at and converted here, which keeps the test data to
stable checksummed URLs rather than a Python package's release cadence.

Positions come out in whatever units MDAnalysis reports, which is Angstrom for every reader
used here, so an absolute error bound means the same thing across all of them. The XTC sets
arrive already quantized to GROMACS's own 0.001 nm grid, so below an error bound of 0.01 the
compressor reproduces them nearly exactly; the DCD and NetCDF sets carry full float precision.

Usage:
    python3 fetch_md_trajectory.py --list
    python3 fetch_md_trajectory.py <dataset> <output_dir> [--fields xyz,x,y,z] [--max-frames N]
"""

import argparse
import hashlib
import os
import random
import time

import numpy as np
import requests

FIGSHARE = "https://ndownloader.figshare.com/files/"

# topology and trajectory for each dataset, as (filename, figshare file id, sha256).
DATASETS = {
    "adk-equilibrium": {
        "topology": ("adk4AKE.psf", "8672230",
                     "1aa947d58fb41b6805dc1e7be4dbe65c6a8f4690f0bd7fc2ae03e7bd437085f4"),
        "trajectory": ("1ake_007-nowater-core-dt240ps.dcd", "8672074",
                       "598fcbcfcc425f6eafbe9997238320fcacc6a4613ecce061e1521732bab734bf"),
        "shape": (4187, 3341),
    },
    "ifabp-water": {
        "topology": ("ifabp_water.psf", "12980639",
                     "ba40714318aabec537015dc550fe5bd5ac1ac0b853f5abdd2f0ae63af9cfcafa"),
        "trajectory": ("rmsfit_ifabp_water_1.dcd", "12980642",
                       "cebb48e58015abc8ff2f5bb7ba3eb7a289047f256351a8252bf1f29f9aaacf0e"),
        "shape": (500, 12445),
    },
    "membrane-peptide": {
        "topology": ("memb_pept.tpr", "14993171",
                     "677a3ae55e35c24f37f2610eafa92d19285d1774731d6ffb9a99dfde39b8c437"),
        "trajectory": ("memb_pept.xtc", "14993174",
                       "f9bdfee4e1aa69ccfeef21cb74703202f6728f514543c4125382bd5250773eb7"),
        "shape": (1001, 18727),
    },
    "cg-fiber": {
        "topology": ("126chains.psf", "13374146",
                     "3ddb654b68549ac2ad5107a4282899f41fad233d09ea572446031711af4e57da"),
        "trajectory": ("126chains.dcd", "13375838",
                       "e0b47d422f31ec209ea810edcf6cf3830da04bb2e1540f520477c27f4433d849"),
        "shape": (2221, 8364),
    },
    "peg-1chain": {
        "topology": ("PEG.prmtop", "13532462",
                     "2d7955b9a8cb6e008171e0c5a1c31e3e458246ea3ee7302281eafefafa7cede9"),
        "trajectory": ("PEG_03_prod.nc", "13532465",
                       "b978714ec2f93d1cbe99564cb257959f0cb38872359aa745c8eba720a7d85225"),
        "shape": (50, 13406),
    },
    "nhaa-equilibrium": {
        "topology": ("NhaA_non_water.gro", "13222709",
                     "ae42f4cfcfe312476f9e5121fe47764a11aff962197799671c0c5a8f83637420"),
        "trajectory": ("NhaA_non_water.xtc", "13222712",
                       "c9ab7ba8c9c271d535cfadebc33da1d90fbf00d9a01f48afedd0f7a703128eaf"),
        "shape": (5001, 60702),
    },
    "yiip-equilibrium": {
        "topology": ("YiiP_system.pdb", "15286808",
                     "3c2b96bbd2f95105e1a4f37140132ee073a947df8fe209a8170f09ca5b73e6cf"),
        "trajectory": ("YiiP_system_9ns_center.xtc", "15285461",
                       "97f6a93acc1e330915338b290625d81d84928a7e4c1e5aed63d209881fbe268b"),
        "shape": (901, 111815),
    },
}


def download(file_id, filename, sha256, into):
    """Fetch one figshare file, reusing it when it is already there and intact."""
    target = os.path.join(into, filename)
    if os.path.exists(target) and file_digest(target) == sha256:
        print(f"{filename} already present")
        return target

    url = FIGSHARE + file_id
    for attempt in range(5):
        try:
            print(f"Downloading {filename} from {url}")
            with requests.get(url, stream=True, timeout=60) as response:
                response.raise_for_status()
                with open(target, "wb") as out:
                    for chunk in response.iter_content(chunk_size=1 << 20):
                        out.write(chunk)
            break
        except (requests.exceptions.RequestException, OSError) as error:
            if attempt == 4:
                raise
            print(f"Download failed ({error}), retrying")
            time.sleep(random.uniform(5, 20))

    digest = file_digest(target)
    if digest != sha256:
        raise SystemExit(f"{filename}: expected sha256 {sha256}, got {digest}")
    return target


def file_digest(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def convert(topology, trajectory, output_dir, expected_shape, fields, max_frames=0):
    """Write xyz.f32 as {frames, atoms, 3}, and one 2D {frames, atoms} file per named coordinate.

    Each coordinate file is a third the size of xyz.f32, so a large trajectory is usually asked
    for one field at a time to stay inside a runner's disk.
    """
    import MDAnalysis  # imported here so --list works without it

    universe = MDAnalysis.Universe(topology, trajectory)
    atoms = len(universe.atoms)
    frames = len(universe.trajectory)
    print(f"{frames} frames, {atoms} atoms")
    if expected_shape and (frames, atoms) != tuple(expected_shape):
        raise SystemExit(f"expected {expected_shape} frames x atoms, found {(frames, atoms)}")
    if max_frames and max_frames < frames:
        print(f"writing the leading {max_frames} frames")
        frames = max_frames

    paths = {name: os.path.join(output_dir, f"{name}.f32") for name in fields}
    handles = {name: open(path, "wb") for name, path in paths.items()}
    try:
        for step, _ in enumerate(universe.trajectory):
            if step == frames:
                break
            positions = universe.atoms.positions.astype(np.float32, copy=False)
            if "xyz" in handles:
                handles["xyz"].write(positions.tobytes())
            for axis, name in enumerate(("x", "y", "z")):
                if name in handles:
                    handles[name].write(np.ascontiguousarray(positions[:, axis]).tobytes())
            if step % 500 == 0:
                print(f"  frame {step}/{frames}")
    finally:
        for handle in handles.values():
            handle.close()

    for name, path in paths.items():
        dims = [frames, atoms, 3] if name == "xyz" else [frames, atoms]
        print(f"  {os.path.basename(path)}: dims {dims}, {os.path.getsize(path)} bytes")
    return frames, atoms


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", nargs="?", help="one of: " + ", ".join(sorted(DATASETS)))
    parser.add_argument("output_dir", nargs="?")
    parser.add_argument("--keep-downloads", action="store_true",
                        help="keep the topology and trajectory files after converting")
    parser.add_argument("--fields", default="xyz,x,y,z",
                        help="which of xyz,x,y,z to write (default all)")
    parser.add_argument("--max-frames", type=int, default=0,
                        help="write only the leading N frames (default: all)")
    parser.add_argument("--list", action="store_true", help="print the dataset names and exit")
    args = parser.parse_args()

    if args.list:
        for name, spec in sorted(DATASETS.items()):
            frames, atoms = spec["shape"]
            print(f"{name:<20} {frames} frames x {atoms} atoms")
        return

    if not args.dataset or not args.output_dir:
        parser.error("dataset and output_dir are required")
    if args.dataset not in DATASETS:
        parser.error(f"unknown dataset {args.dataset}; try --list")

    spec = DATASETS[args.dataset]
    os.makedirs(args.output_dir, exist_ok=True)

    sources = [download(file_id, filename, checksum, args.output_dir)
               for filename, file_id, checksum in (spec["topology"], spec["trajectory"])]
    fields = [f for f in args.fields.split(",") if f]
    unknown = [f for f in fields if f not in ("xyz", "x", "y", "z")]
    if unknown:
        parser.error(f"unknown fields {unknown}")
    convert(sources[0], sources[1], args.output_dir, spec["shape"], fields, args.max_frames)

    if not args.keep_downloads:
        for path in sources:
            os.remove(path)
            print(f"Removed {os.path.basename(path)}")


if __name__ == "__main__":
    main()
