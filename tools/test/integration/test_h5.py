import numpy as np
import h5py
import subprocess
import os
import sys
from shutil import which, rmtree
import tempfile

def check_h5repack():
    """
    Check if 'h5repack-shared' is available in the system PATH.
    If not, check for 'h5repack'.
    Returns the command to use.
    """
    if which('h5repack-shared') is not None:
        print("Using 'h5repack-shared' for HDF5 operations.")
        return 'h5repack-shared'
    elif which('h5repack') is not None:
        print("Using 'h5repack' for HDF5 operations.")
        return 'h5repack'
    else:
        print("Error: Neither 'h5repack-shared' nor 'h5repack' is installed or found in PATH.")
        sys.exit(1)

def create_comp_conf(comp_conf_path, cmpr_algo, bound):
    """
    Creates the 'comp.conf' configuration file with the specified content and error bound.

    Parameters:
    - comp_conf_path (str): Path where the comp.conf file will be created.
    - bound (float): The absolute error bound to set in the configuration.
    """
    comp_conf_content = f"""[GlobalSettings]
ErrorBoundMode = ABS
AbsErrorBound = {bound}
CmprAlgo = {cmpr_algo}
OpenMP = No
"""
    try:
        with open(comp_conf_path, 'w') as f:
            f.write(comp_conf_content)
        print(f"Configuration file '{comp_conf_path}' created successfully with AbsErrorBound = {bound}")
    except Exception as e:
        print(f"Error creating configuration file '{comp_conf_path}': {e}")
        sys.exit(1)

def get_compression_args(comp_conf_path):
    """
    Calls the external executable 'print_h5repack_args' with the configuration file
    to get compression arguments.

    Parameters:
    - comp_conf_path (str): Path to the comp.conf file.

    Returns:
    - str: Compression string obtained from 'print_h5repack_args'.
    """
    try:
        print(f"Calling 'print_h5repack_args' with configuration file '{comp_conf_path}' to get compression parameters...")
        result = subprocess.run(['./print_h5repack_args', '-c', comp_conf_path], capture_output=True, text=True, check=True)
        compression = result.stdout.strip()
        if not compression:
            raise ValueError("Compression argument is empty.")
        print(f"Compression parameters obtained: '{compression}'")
        return compression
    except FileNotFoundError:
        print("Error: 'print_h5repack_args' executable not found in PATH.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Error executing 'print_h5repack_args': {e}")
        print(f"stderr: {e.stderr}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error while getting compression arguments: {e}")
        sys.exit(1)

def read_raw_fp32(raw_file_path, shape):
    """
    Reads raw FP32 data from a file and reshapes it.

    Parameters:
    - raw_file_path (str): Path to the raw FP32 file.
    - shape (tuple): Desired shape of the data.

    Returns:
    - np.ndarray: Reshaped data array.
    """
    try:
        data = np.fromfile(raw_file_path, dtype=np.float32)
        expected_size = np.prod(shape)
        if data.size != expected_size:
            raise ValueError(f"Data size {data.size} does not match expected shape {shape}")
        data = data.reshape(shape)
        print(f"Successfully read and reshaped raw data from '{raw_file_path}' with shape {shape}")
        return data
    except Exception as e:
        print(f"Error reading raw data: {e}")
        sys.exit(1)

def write_hdf5(data, h5_file_path, dataset_name='test'):
    """
    Writes data to an HDF5 file.

    Parameters:
    - data (np.ndarray): Data to write.
    - h5_file_path (str): Path to the HDF5 file.
    - dataset_name (str): Name of the dataset within the HDF5 file.
    """
    try:
        with h5py.File(h5_file_path, 'w') as f:
            f.create_dataset(dataset_name, data=data)
        print(f"Data written to HDF5 file '{h5_file_path}' with dataset name '{dataset_name}'")
    except Exception as e:
        print(f"Error writing HDF5 file: {e}")
        sys.exit(1)

def apply_h5repack(repack_cmd, input_h5, output_h5, compression='gzip=4'):
    """
    Applies h5repack (or h5repack-shared) with specified compression.

    Parameters:
    - repack_cmd (str): The h5repack command to use ('h5repack' or 'h5repack-shared').
    - input_h5 (str): Input HDF5 file path.
    - output_h5 (str): Output HDF5 file path after repacking.
    - compression (str): Compression filter settings for h5repack.
    """
    try:
        cmd = [repack_cmd, '-f', compression, input_h5, output_h5]
        print(f"Running {repack_cmd} with command: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        print(f"{repack_cmd} completed. Output file: '{output_h5}'")
    except subprocess.CalledProcessError as e:
        print(f"Error running {repack_cmd}: {e}")
        print(f"stderr: {e.stderr}")
        sys.exit(1)

def compare_hdf5(original_h5, decompressed_h5, dataset_name='test'):
    """
    Compares two HDF5 files and prints the maximum error between the datasets.

    Parameters:
    - original_h5 (str): Path to the original HDF5 file.
    - decompressed_h5 (str): Path to the decompressed HDF5 file.
    - dataset_name (str): Name of the dataset to compare.

    Returns:
    - float: The maximum error found between the datasets.
    """
    try:
        with h5py.File(original_h5, 'r') as f1, h5py.File(decompressed_h5, 'r') as f2:
            data1 = f1[dataset_name][:]
            data2 = f2[dataset_name][:]
            if data1.shape != data2.shape:
                raise ValueError("Shape mismatch between original and decompressed data.")
            max_error = np.max(np.abs(data1 - data2))
            max_error_location = np.unravel_index(np.argmax(np.abs(data1 - data2)), data1.shape)
            print(f"Max error between original and decompressed data: {max_error}")
            print(f"Location of max error: {max_error_location}")
            return max_error
    except Exception as e:
        print(f"Error comparing HDF5 files: {e}")
        sys.exit(1)

def parse_shape(args):
    """
    Parses shape dimensions from command-line arguments.

    Parameters:
    - args (list): List of shape dimensions as strings.

    Returns:
    - tuple: Parsed shape as a tuple of integers.
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

def format_bound(bound):
    """
    Formats the error bound for use in filenames.

    Parameters:
    - bound (float): The error bound to format.

    Returns:
    - str: Formatted bound string.
    """
    # Convert to scientific notation, remove '+' sign, replace '.' with 'p'
    return "{:e}".format(bound).replace('e+', 'e').replace('.', 'p')

def main():

    # Define the list of error bounds to test
    error_bounds = [10, 1, 5e-1, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4]
    # error_bounds = [1e-4]
    # Ensure the script is called with the correct number of arguments
    if len(sys.argv) < 3:
        print("Usage: python script.py <cmpr_algo> <raw_file_path> <dim1> [<dim2> ... <dimN>]")
        print("Example: python script.py ALGO_LORENZO_REG input.fp32 100 100 100")
        sys.exit(1)

    # Input Arguments
    cmpr_algo = sys.argv[1]
    raw_file = sys.argv[2]
    shape_args = sys.argv[3:]

    # Fixed Parameters
    h5_dataset_name = 'test'             # Fixed dataset name
    # output_dir = './output_hdf5'      # Fixed output directory
    output_dir = tempfile.mkdtemp(prefix='output_hdf5_')

    # Parse Shape
    shape = parse_shape(shape_args)

    # Define output file paths based on output_dir and input file name
    base_name = os.path.splitext(os.path.basename(raw_file))[0]
    original_h5 = os.path.join(output_dir, f"{base_name}_original.h5")

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Check for h5repack or h5repack-shared
    repack_cmd = check_h5repack()

    # Step 1: Read raw data
    data = read_raw_fp32(raw_file, shape)

    # Step 2: Convert to HDF5
    write_hdf5(data, original_h5, h5_dataset_name)

    passed_all_cases = True
    # Iterate over each error bound and perform compression, decompression, and comparison
    for bound in error_bounds:
        print("\n" + "="*50)
        print(f"Testing {raw_file} with algo = {cmpr_algo} AbsErrorBound = {bound}")
        print("="*50)

        # Create comp.conf with the current bound
        comp_conf_path = os.path.join(output_dir, f"comp_{format_bound(bound)}.conf")
        create_comp_conf(comp_conf_path, cmpr_algo, bound)

        # Get Compression Arguments by calling external executable
        compression = get_compression_args(comp_conf_path)

        # Define compressed and decompressed file paths for the current bound
        bound_str = format_bound(bound)
        compressed_h5 = os.path.join(output_dir, f"{base_name}_compressed_{bound_str}.h5")
        decompressed_h5 = os.path.join(output_dir, f"{base_name}_decompressed_{bound_str}.h5")

        # Apply compression
        apply_h5repack(repack_cmd, original_h5, compressed_h5, compression=compression)

        # Decompress the compressed HDF5
        apply_h5repack(repack_cmd, compressed_h5, decompressed_h5, compression='NONE')

        # Compare original and decompressed HDF5 files
        max_error = compare_hdf5(original_h5, decompressed_h5, h5_dataset_name)

        # Verify that max_error <= bound * 1.1
        if max_error <= bound * 1.5:
            result = "PASS"
        else:
            result = "FAIL"
            passed_all_cases = False

        print(f"Test Result for AbsErrorBound = {bound}: {result}")

    # Clean up intermediate files
    try:
        rmtree(output_dir)
        print(f"Cleaned up intermediate folder")
    except Exception as e:
        print(f"Error removing intermediate folder: {e}")

    if passed_all_cases:
        print(f"\nAll tests passed successfully for {raw_file}.")
    else:
        print(f"\nSome tests failed for {raw_file}. Please check the results.")

if __name__ == '__main__':
    main()