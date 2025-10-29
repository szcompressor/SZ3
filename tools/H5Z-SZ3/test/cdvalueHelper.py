#!/usr/bin/env python3
"""
Convert SZ3 parameters to HDF5 cd_values
Note: many configuration parameters are set to placeholder values
      and will be overwritten when the HDF5 filter calls the set_local function.
"""

import sys
import struct

class SZ3:
    def __init__(self, absolute=None, relative=None, psnr=None, l2norm=None):
        if sum(x is not None for x in [absolute, relative, psnr, l2norm]) != 1:
            raise ValueError("Specify exactly one error bound: absolute, relative, psnr, or l2norm.")

        # placeholder values will be overwritten when the HDF5 filter calls the set_local function,
        serialized = bytearray()

        serialized.extend(struct.pack('<B', 40)) # size of the config structure
        serialized.extend(struct.pack('<I', 0))  # placeholder for SZ3_MAGIC_NUMBER
        serialized.extend(struct.pack('<I', 0))  # placeholder for SZ3_DATA_VER
        serialized.extend(struct.pack('<B', 0))  # placeholder for dimension count
        serialized.extend(struct.pack('<B', 0))  # placeholder for dimension value bit width
        serialized.extend(struct.pack('<Q', 0))  # placeholder for total data elements
        serialized.extend(struct.pack('<B', 1))  # compression algorithm, 1: ALGO_INTERP_LORENZO

        if absolute is not None:
            serialized.extend(struct.pack('<B', 0)) # Error bound mode for EB_ABS
            serialized.extend(struct.pack('<d', absolute))
        elif relative is not None:
            serialized.extend(struct.pack('<B', 1)) # Error bound mode for EB_REL
            serialized.extend(struct.pack('<d', relative))
        elif psnr is not None:
            serialized.extend(struct.pack('<B', 2)) # Error bound mode for EB_PSNR
            serialized.extend(struct.pack('<d', psnr))
        elif l2norm is not None:
            serialized.extend(struct.pack('<B', 3)) # Error bound mode for EB_L2NORM
            serialized.extend(struct.pack('<d', l2norm))


        serialized.extend(struct.pack('<B', 0)) # boolean flags
        serialized.extend(struct.pack('<B', 0))  # placeholder for data type
        serialized.extend(struct.pack('<I', 65536))  # quantization bin count
        serialized.extend(struct.pack('<I', 0))  # placeholder for block size
        serialized.extend(struct.pack('<B', 0))  # prediction dimension

        serialized_data = bytes(serialized)
        self.cd_values = [int.from_bytes(serialized_data[i:i + 4], 'little')
                     for i in range(0, len(serialized_data), 4)]


    def print_h5repack_args(self):
        cd_nelmts = len(self.cd_values)
        args = f"-f UD=32024,0,{cd_nelmts},"
        args += ",".join(map(str, self.cd_values))
        print(args)


def main():
    """Command-line interface."""
    import argparse

    parser = argparse.ArgumentParser(
        description='Convert SZ3 error bounds to HDF5 cd_values',
        epilog="""
Examples:
  %(prog)s --abs 1e-3
  %(prog)s --rel 0.001
  %(prog)s --psnr 80
  %(prog)s --l2norm 0.01
  
Use with h5repack:
  h5repack $(%(prog)s --abs 1e-3) input.h5 output.h5
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--abs', '--absolute', type=float, dest='absolute',
                      help='Absolute error bound')
    group.add_argument('--rel', '--relative', type=float, dest='relative',
                      help='Relative error bound (value-range based)')
    group.add_argument('--psnr', type=float,
                      help='PSNR error bound')
    group.add_argument('--l2norm', '--norm', type=float, dest='l2norm',
                      help='L2 norm error bound')

    args = parser.parse_args()

    try:
        if args.absolute is not None:
            config = SZ3(absolute=args.absolute)
        elif args.relative is not None:
            config = SZ3(relative=args.relative)
        elif args.psnr is not None:
            config = SZ3(psnr=args.psnr)
        else:
            config = SZ3(l2norm=args.l2norm)

        config.print_h5repack_args()
        return 0

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())

