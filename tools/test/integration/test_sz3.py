import subprocess
import re
import os
from pathlib import Path
import sys

HOME = str(Path.home())


def run_command(cmd):
    """Run a command and return stdout and stderr."""
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout, result.stderr


def parse_errors(output):
    """
    Parse the decompression output for maximum errors.
    Expected output lines:
       Max absolute error = <value>
       Max relative error = <value>
    Returns (abs_error, rel_error) as floats (or None if not found).
    """
    abs_err = None
    rel_err = None
    abs_match = re.search(r"Max absolute error\s*=\s*([\deE\.\+-]+)", output)
    if abs_match:
        abs_err = float(abs_match.group(1))
    rel_match = re.search(r"Max relative error\s*=\s*([\deE\.\+-]+)", output)
    if rel_match:
        rel_err = float(rel_match.group(1))
    return abs_err, rel_err


def test_file(sz3_path, file_path, dims, mode, error_bound):
    """
    Run the SZ3 command with the given file, reversed dimensions, error control mode, and error bound.
    - file_path: path to the original binary file.
    - dims: a list of dimensions (e.g. [100, 500, 500] from the test file).
    - mode: either "ABS" or "REL".
    - error_bound: error bound value to test.

    Returns: (abs_err, rel_err) from the output.
    """
    # Derive filename for the decompressed file
    base = os.path.splitext(os.path.basename(file_path))[0]
    decompressed_file = "{}.{}.{}.sz.out".format(base, mode, error_bound)

    # Reverse the dimensions order for SZ3 (e.g., 100 500 500 becomes 500 500 100)
    dims_str = " ".join(str(d) for d in dims[::-1])

    # Build the single command that compresses and decompresses the file.
    # Note: no -z is provided, only -i (input) and -o (decompressed output).
    cmd = f"{sz3_path} -f -i {file_path} -o {decompressed_file} -3 {dims_str} -M {mode} {error_bound} -a"

    print("Running command:")
    print(cmd)
    stdout, stderr = run_command(cmd)

    if stderr:
        print("stderr:", stderr)

    print("Command output:")
    print(stdout)

    # Parse error values from the output
    abs_err, rel_err = parse_errors(stdout)

    if os.path.exists(decompressed_file):
        os.remove(decompressed_file)

    return abs_err, rel_err

def load_test_files(file_list_path):
    """
    Read the list of test files from the provided text file.
    Each line should have: <file_path> <dim1> <dim2> ... <dimN>
    Returns a list of tuples: (file_path, [dim1, dim2, ...])
    """
    test_files = []
    with open(file_list_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            file_path = parts[0].replace('~', HOME)
            try:
                dims = [int(x) for x in parts[1:]]
            except ValueError:
                print(f"Error parsing dimensions for line: {line}")
                continue
            test_files.append((file_path, dims))
    return test_files


def main():
    if len(sys.argv) > 2:
        sz3_path = sys.argv[1]
        file_list_path = sys.argv[2]
    else:
        print("Usage: python test_sz3.py <path_to_sz3_executable> <path_to_datalist>")
        exit()

    test_files = load_test_files(file_list_path)

    # Error bounds to test for each mode
    error_tests = {
        "ABS": [1, 1e-3, 1e-5, 1e-7],
        "REL": [1e-2, 1e-4, 1e-6]
    }

    all_passed = True
    for file_path, dims in test_files:
        if not os.path.exists(file_path):
            print(f"Test file {file_path} does not exist, skipping.")
            continue
        for mode, bounds in error_tests.items():
            for bound in bounds:
                print("=" * 60)
                print(f"Testing file: {file_path}\nMode: {mode}\nError Bound: {bound}")
                abs_err, rel_err = test_file(sz3_path, file_path, dims, mode, bound)

                # Check results depending on the mode
                if mode == "ABS":
                    if abs_err is not None and abs_err <= bound:
                        print(f"PASS: Max absolute error {abs_err} is within bound {bound}")
                    else:
                        print(f"FAIL: Max absolute error {abs_err} exceeds bound {bound}")
                        all_passed = False
                elif mode == "REL":
                    if rel_err is not None and rel_err <= bound:
                        print(f"PASS: Max relative error {rel_err} is within bound {bound}")
                    else:
                        print(f"FAIL: Max relative error {rel_err} exceeds bound {bound}")
                        all_passed = False
                print("=" * 60 + "\n")

    if all_passed:
        print("All tests passed.")
    else:
        print("Some tests failed.")


if __name__ == "__main__":
    main()
