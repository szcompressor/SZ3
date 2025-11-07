import numpy as np
import os
import sys
from shutil import rmtree
import tempfile

sys.path.append(os.path.join(os.path.dirname(__file__), '../../H5Z-SZ3/test'))
from cdvalueHelper import SZ3

if len(sys.argv) > 1:
    os.environ["HDF5_PLUGIN_PATH"] = sys.argv[1]
import h5py

if not h5py.h5z.filter_avail(32024):
    print("SZ3 filter not available")
    sys.exit(1)

def get_compression_args(cmpr_algo, bound):
    """
    Generates compression arguments using cdvalueHelper.SZ3.

    Parameters:
    - cmpr_algo (str): The compression algorithm.
    - bound (float): The absolute error bound.

    Returns:
    - tuple: Compression filter and options.
    """
    try:
        config = SZ3(cmpr_algo, absolute=bound)
        return (32024, tuple(config.cd_values))
    except Exception as e:
        print(f"Error generating compression arguments: {e}")
        sys.exit(1)

def read_raw_data(raw_file_path, shape, dtype):
    """
    Reads raw data from a file and reshapes it.

    Parameters:
    - raw_file_path (str): Path to the raw file.
    - shape (tuple): Desired shape of the data.
    - dtype (str): Data type, 'float32' or 'float64'.

    Returns:
    - np.ndarray: Reshaped data array.
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

def write_hdf5(data, h5_file_path, dataset_name='test', compression=None, compression_opts=None, chunks=None):
    """
    Writes data to an HDF5 file.

    Parameters:
    - data (np.ndarray): Data to write.
    - h5_file_path (str): Path to the HDF5 file.
    - dataset_name (str): Name of the dataset within the HDF5 file.
    - compression (str): Compression filter to use.
    - compression_opts (tuple): Compression options.
    - chunks (tuple or None): Chunk shape for the dataset.
    """
    try:
        with h5py.File(h5_file_path, 'w') as f:
            f.create_dataset(dataset_name, data=data, compression=compression, compression_opts=compression_opts, chunks=chunks)
        print(f"Data written to HDF5 file '{h5_file_path}' with dataset name '{dataset_name}'")
    except Exception as e:
        print(f"Error writing HDF5 file: {e}")
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


def main():
    if len(sys.argv) < 7:
        print("Usage: python test_h5_filter.py <h5_plugin_path> <cmpr_algo> <error_bound> <raw_file_path> <dtype> <dim1> [<dim2> ... <dimN>]")
        sys.exit(1)

    h5_plugin_path = sys.argv[1]
    cmpr_algo = sys.argv[2]
    bound = float(sys.argv[3])
    raw_file = sys.argv[4]
    dtype = sys.argv[5]
    shape_args = sys.argv[6:]

    h5_dataset_name = 'testdata_compressed'
    output_dir = tempfile.mkdtemp(prefix='output_hdf5_')
    # output_dir='/var/tmp/output_hdf5_'
    # os.makedirs(output_dir, exist_ok=True)

    shape = parse_shape(shape_args)

    base_name = os.path.splitext(os.path.basename(raw_file))[0]
    original_h5 = os.path.join(output_dir, f"{base_name}_original.h5")

    data = read_raw_data(raw_file, shape, dtype)

    write_hdf5(data, original_h5, h5_dataset_name)

    compression, compression_opts = get_compression_args(cmpr_algo, bound)


    all_pass = True
    for chunk in [False, True]:
        print(f"Testing {raw_file} with algo = {cmpr_algo} AbsErrorBound = {bound} Chunk = {chunk}")

        compressed_h5 = os.path.join(output_dir, f"{base_name}_compressed_{chunk}.h5")
        decompressed_h5 = os.path.join(output_dir, f"{base_name}_decompressed_{chunk}.h5")

        if chunk:
            # hd5py will automatically determine chunk sizes if chunks is not set
            write_hdf5(data, compressed_h5, h5_dataset_name, compression=compression, compression_opts=compression_opts)
        else:
            write_hdf5(data, compressed_h5, h5_dataset_name, compression=compression, compression_opts=compression_opts, chunks=shape)
        

        with h5py.File(compressed_h5, 'r') as f_in, h5py.File(decompressed_h5, 'w') as f_out:
            f_out.create_dataset(h5_dataset_name, data=f_in[h5_dataset_name][:])

        max_error = compare_hdf5(original_h5, decompressed_h5, h5_dataset_name)

        if max_error <= bound * 1.5:
            result = "PASS"
        else:
            result = "FAIL"

        print(f"Test Result for AbsErrorBound = {bound} ChunkSize = {chunk}: {result}")

        if result == "FAIL":
            all_pass = False

    rmtree(output_dir)

    if not all_pass:
        sys.exit(1)

if __name__ == '__main__':
    main()
