import numpy as np
import subprocess
import os
import sys
from shutil import rmtree
import tempfile


def create_sz3_config(algo, path):
    """
    Creates the sz3.config file.
    """
    with open(os.path.join(path, "sz3.config"), "w") as f:
        f.write(f"[GlobalSettings]\n")
        f.write(f"CmprAlgo = {algo}\n")


def run_sz3_compress(sz3_executable, input_file, output_file, bound, dims, cwd, dtype_flag):
    """
    Runs the sz3 executable for compression.
    """
    try:
        cmd = [sz3_executable, dtype_flag, '-i', input_file, '-z', os.path.basename(output_file), '-c', 'sz3.config', '-M',
               'ABS', str(bound), f'-{len(dims)}'] + dims
        print(f"Running sz3 compression with command: {' '.join(cmd)}")
        subprocess.run(cmd, check=True, cwd=cwd)
        print(f"sz3 compression completed. Output file: '{output_file}'")
    except subprocess.CalledProcessError as e:
        print(f"Error running sz3 compression: {e}")
        print(f"stderr: {e.stderr}")
        sys.exit(1)


def run_sz3_decompress(sz3_executable, compressed_file, decompressed_file, dims, cwd, dtype_flag):
    """
    Runs the sz3 executable for decompression.
    """
    try:
        cmd = [sz3_executable, dtype_flag, '-s', os.path.basename(compressed_file), '-o',
               os.path.basename(decompressed_file), f'-{len(dims)}'] + dims
        print(f"Running sz3 decompression with command: {' '.join(cmd)}")
        subprocess.run(cmd, check=True, cwd=cwd)
        print(f"sz3 decompression completed. Output file: '{decompressed_file}'")
    except subprocess.CalledProcessError as e:
        print(f"Error running sz3 decompression: {e}")
        print(f"stderr: {e.stderr}")
        sys.exit(1)


def read_raw_data(raw_file_path, shape, dtype):
    """
    Reads raw data from a file and reshapes it.
    """
    try:
        # Mapping from dtype string to numpy dtype
        np_dtypes = {
            'float32': np.float32,
            'float64': np.float64,
            'int32': np.int32,
            'int64': np.int64,
            # Add more as needed
        }
        np_dtype = np_dtypes.get(dtype, np.float32)  # Default to float32
        data = np.fromfile(raw_file_path, dtype=np_dtype)
        expected_size = np.prod(shape)
        if data.size != expected_size:
            raise ValueError(f"Data size {data.size} does not match expected shape {shape}")
        data = data.reshape(shape)
        print(f"Successfully read and reshaped raw data from '{raw_file_path}' with shape {shape}")
        return data
    except Exception as e:
        print(f"Error reading raw data: {e}")
        sys.exit(1)


def compare_data(original_data, decompressed_data):
    """
    Compares two numpy arrays and prints the maximum error.
    """
    max_error = np.max(np.abs(original_data - decompressed_data))
    print(f"Max error between original and decompressed data: {max_error}")
    return max_error


def parse_shape(args):
    """
    Parses shape dimensions from command-line arguments.
    """
    try:
        shape = tuple(int(dim) for dim in args)
        return shape
    except ValueError:
        print("Error: All shape dimensions must be integers.")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error parsing shape: {e}")
        sys.exit(1)


def main():
    if len(sys.argv) < 7:
        print(
            "Usage: python test_sz3_executable.py <sz3_exe_folder> <cmpr_algo> <error_bound> <raw_file_path> <dtype> <dim1> [<dim2> ... <dimN>]")
        sys.exit(1)

    sz3_exe_folder = sys.argv[1]
    cmpr_algo = sys.argv[2]
    bound = float(sys.argv[3])
    raw_file = sys.argv[4]
    dtype = sys.argv[5]
    shape_args = sys.argv[6:]

    # Mapping from dtype string to SZ3 flag
    dtype_flags = {
        'float32': '-f',
        'float64': '-d',
        'int32': '-i',
        'int64': '-l',
        # Add more as needed
    }
    dtype_flag = dtype_flags.get(dtype, '-f')  # Default to -f if unknown

    sz3_executable = os.path.join(sz3_exe_folder, "sz3") if os.path.isdir(sz3_exe_folder) else sz3_exe_folder

    output_dir = tempfile.mkdtemp(prefix='output_sz3_')

    shape = parse_shape(shape_args)

    base_name = os.path.splitext(os.path.basename(raw_file))[0]
    compressed_file = os.path.join(output_dir, f"{base_name}_compressed.sz3")
    decompressed_file = os.path.join(output_dir, f"{base_name}_decompressed.dat")

    create_sz3_config(cmpr_algo, output_dir)

    original_wd = os.getcwd()
    os.chdir(output_dir)

    run_sz3_compress(sz3_executable, raw_file, compressed_file, bound, [str(d) for d in shape[::-1]], output_dir, dtype_flag)

    run_sz3_decompress(sz3_executable, compressed_file, decompressed_file, [str(d) for d in shape[::-1]], output_dir, dtype_flag)

    os.chdir(original_wd)

    original_data = read_raw_data(raw_file, shape, dtype)
    decompressed_data = read_raw_data(decompressed_file, shape, dtype)

    max_error = compare_data(original_data, decompressed_data)

    if max_error <= bound * 1.5:
        result = "PASS"
    else:
        result = "FAIL"

    print(f"Test Result for AbsErrorBound = {bound}: {result}")

    rmtree(output_dir)

    if result == "FAIL":
        sys.exit(1)


if __name__ == '__main__':
    main()
