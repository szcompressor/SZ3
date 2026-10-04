# Compressing molecular-dynamics data with SZ3

## Which compressor for which data

| Data | Use | Notes |
|---|---|---|
| Biomolecular trajectories: coordinates `{frames, atoms, 3}` (GROMACS, H5MD) | `ALGO_BIOMD` | Absolute bound, strict. Uses water and bond geometry. |
| The same coordinates, when you want xtc's behaviour | `ALGO_BIOMDXTC` | Follows GROMACS's xtc; a coordinate can come back up to 10% past the bound. |
| Velocities and forces | `ALGO_NOPRED` | Absolute bound. Neighbours and earlier frames do not predict them at usual output intervals. |
| Box, step, time, energies | `ALGO_LOSSLESS`, or HDF5 shuffle + gzip | A few values per frame, needed exactly. |
| MD of solid materials (crystals, metals) | MDZ, [tools/mdz](../tools/mdz/README.md) | Atoms vibrating around lattice sites. |
| Fields on a grid (densities, potentials) | the default `ALGO_INTERP_LORENZO` | |

## How well ALGO_BIOMD does

At XTC's default precision, 5e-4 nm, on 101 trajectories (geometric means; 1 / 20 frames per chunk):

- **Compression ratio.** All-atom with water: 3.81 / 5.05, against 3.31 for XTC; the best other HDF5 filter within
  the bound gives 2.33 / 3.18, lossless filters 1.14–1.28. All-atom without water: 3.27 / 5.36, against 3.12 for
  XTC. Coarse-grained: 2.32 / 3.40, against 2.39 for XTC.
- **Speed.** Compression takes 1.74 / 1.36 times XTC's time, decompression 1.34 / 1.44 times. Through HDF5 it writes
  439 / 498 MB/s and reads 574 / 572 MB/s on one core.
- **In GROMACS.** H5MD files of all-atom systems are 0.86–0.87 / 0.62–0.66 times the size of XTC files, and writing
  them takes 0.09%–0.62% / 0.09%–0.60% of a GPU run (XTC 0.06%–0.43%).
- **Error.** Every coordinate comes back within the bound.

Every run, the other bounds, the setup and the comparison with other HDF5 filters are in
[#169](https://github.com/szcompressor/SZ3/pull/169).

## ALGO_BIOMD

- **Shape.** `{frames, atoms, 3}`; one frame (`{atoms, 3}`) also works. Other shapes throw `std::invalid_argument`. A
  single call or HDF5 chunk holds at most 2^31 - 1 values.
- **Bound.** Absolute (`EB_ABS`), in the units of the coordinates. MD trajectories are usually kept at 5e-4 nm (XTC's
  default precision) or 5e-3 nm for coarse-grained runs. float and double.
- **Values stored losslessly.** Coordinates go on an integer lattice that covers `|x| < B`, B the smallest power of
  two above 2^21 times the bound for float (2^27 for double). A chunk holding a coordinate with `|x| >= B`, or NaN or
  Inf outside trailing fill frames (frames at the end of the chunk in which every value is the same), is stored
  losslessly, at a ratio of about 1.1 to 1.5 instead of 3 to 6:

  | Bound (nm) | 1e-5 | 5e-5 | 1e-4 | 5e-4 | 1e-3 | 5e-3 |
  |---|---|---|---|---|---|---|
  | B for float (nm) | 32 | 128 | 256 | 2048 | 4096 | 16384 |

- **Recompression.** Decompressing and compressing again with the same bound gives the same values, so a chunk HDF5
  writes again stays within the bound.
- **Frames per chunk.** Each chunk is compressed on its own, and BIOMD finds water and bonds on its first frame. More
  frames per chunk give a higher ratio (see above). Coarse-grained runs at one frame per chunk come out below XTC:
  most beads have no bonds or water to predict from.

## The HDF5 filter for trajectories

[tools/H5Z-SZ3/README.md](../tools/H5Z-SZ3/README.md) covers installing and using the filter. For trajectories:

- **Chunk shape.** An `ALGO_BIOMD` chunk must have 3 as its last dimension and should hold whole frames,
  `(k, atoms, 3)`. Set `chunks` yourself; h5py's automatic chunking splits the atoms and xyz, and creating such a
  dataset fails.
- **Registering the filter.** An application that links the filter calls `H5Zregister(H5PLget_plugin_info())` when it
  opens any file, so that a program that only reads can read SZ3 data.
- **Appending frames.** Set `H5D_CHUNK_DONT_FILTER_PARTIAL_CHUNKS` and give the dataset a chunk cache that holds a
  whole chunk (`H5Pset_chunk_cache`), also when the file is reopened. Otherwise HDF5 compresses a chunk again each
  time frames are added to it. The last, partial chunk is stored uncompressed at full chunk size.
- **Other filters.** Add `H5Pset_fletcher32` after SZ3 to detect damaged chunks; never put shuffle before SZ3.

## Which setting for which GROMACS output

| Output | SZ3 setting | Compression ratio, 1 / 20 frames per chunk |
|---|---|---|
| positions | `ALGO_BIOMD`, absolute bound (5e-4 nm is XTC's default), chunks of `(frames, atoms, 3)` | 3.81 / 5.05 (all-atom with water) |
| velocities | `ALGO_NOPRED`, absolute bound (for example 1e-3 nm/ps) | 2.90–2.96 / 2.94–2.99 (Zstd 1.07) |
| forces | `ALGO_NOPRED`, absolute bound (for example 1 kJ/mol/nm) | 3.08–3.68 / 3.10–3.70 (Zstd 1.08–1.37) |
| box, step, time, energies | `ALGO_LOSSLESS`, or shuffle + gzip | – |
