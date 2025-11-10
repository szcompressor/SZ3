
import subprocess
import os
import sys
import json
import requests
import tarfile
import tempfile
import shutil

def get_tmpdir():
    candidates = ['/scratch', '/var/tmp', '/tmp']
    for cand in candidates:
        if os.path.exists(cand):
            try:
                test_dir = os.path.join(cand, 'sz3_test_write')
                os.makedirs(test_dir, exist_ok=True)
                os.rmdir(test_dir)
                return cand
            except OSError:
                pass
    return tempfile.gettempdir()

def run_test(cmd, description):
    print('\n\n', "="*80)
    print(description)
    print("command:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    
    if result.stderr:
        print("STDERR:", result.stderr)
    if result.returncode == 0:
        print("PASS")
        print("="*80)
        return True
    else:
        print("FAIL")
        print("="*80)
        return False

def prepare_dataset(path, dataset_dir, dataset_info=None):
    """
    Prepares the dataset: if path is URL, download and extract; if local dir, copy it
    Returns the actual directory that directly containing the files.
    """
    if path.startswith('http'):
        # Download and extract
        if not os.path.exists(dataset_dir):
            os.makedirs(dataset_dir)

        tar_filename = os.path.join(dataset_dir, os.path.basename(path))

        if not os.path.exists(tar_filename):
            print(f"Downloading {path} to {tar_filename}")
            with requests.get(path, stream=True) as r:
                r.raise_for_status()
                with open(tar_filename, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
        else:
            print(f"{tar_filename} already exists. Skipping download.")

        print(f"Extracting {tar_filename} to {dataset_dir}")
        with tarfile.open(tar_filename, 'r:gz') as tar_ref:
            tar_ref.extractall(dataset_dir, filter='data')
        
        # Delete the tar file to save space if running in GitHub CI
        if os.getenv('GITHUB_ACTIONS') == 'true':
            os.remove(tar_filename)
            print(f"Deleted {tar_filename} to save space in CI")
        
        # Find the data directory - if there's a single subdirectory, use it
        subdirs = [d for d in os.listdir(dataset_dir) if os.path.isdir(os.path.join(dataset_dir, d))]
        if len(subdirs) == 1:
            return os.path.join(dataset_dir, subdirs[0])
        else:
            return dataset_dir
    else:
        # Local path - assume it's a directory
        if os.path.isdir(path):
            if not os.path.exists(dataset_dir):
                os.makedirs(dataset_dir)
            if dataset_info:
                for field in dataset_info["fields"].keys():
                    src = os.path.join(path, field)
                    dst = os.path.join(dataset_dir, field)
                    if os.path.isfile(src):
                        shutil.copy2(src, dst)
                    else:
                        print(f"Warning: field file {src} not found")
            else:
                shutil.copytree(path, dataset_dir)
            return dataset_dir
    return dataset_dir


def main():
    if len(sys.argv) > 1:
        datasets_json = sys.argv[1]
    else:
        print("No datasets.json provided ")
        sys.exit(1)

    try:
        with open(datasets_json, 'r') as f:
            datasets = json.load(f)
    except FileNotFoundError:
        print(f"Error: Datasets file '{datasets_json}' not found.")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON in '{datasets_json}': {e}")
        sys.exit(1)

    if len(sys.argv) > 2:
        selected_dataset = sys.argv[2]
        if selected_dataset not in datasets:
            print(f"Dataset {selected_dataset} not found in {datasets_json}")
            sys.exit(1)
        datasets = {selected_dataset: datasets[selected_dataset]}

    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_source_dir = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))

    test_dir = os.path.join(get_tmpdir(), "sz3_integration_test")
    data_dir = os.path.join(test_dir, "data")
    build_dir = os.path.join(test_dir, "build")

    for d in [test_dir, data_dir]:
        os.makedirs(d, exist_ok=True)

    original_wd = os.getcwd()
    os.chdir(test_dir)
    
    # Build the project
    os.makedirs(build_dir, exist_ok=True)
    subprocess.run(["cmake", project_source_dir, "-DBUILD_H5Z_FILTER=ON"], cwd=build_dir, check=True)
    subprocess.run(["cmake", "--build", ".", "--", "-j"], cwd=build_dir, check=True)
    sz3_executable_path = os.path.join(build_dir, "tools", "sz3")
    h5_plugin_path = os.path.join(build_dir, "tools", "H5Z-SZ3")

    error_bounds = [1e-1, 1e-2, 1e-3, 1e-4]
    algorithms = ["ALGO_INTERP_LORENZO", "ALGO_LORENZO_REG", "ALGO_BIOMD", "ALGO_BIOMDXTC"]

    results = []

    for dataset_name, dataset_info in datasets.items():
        dataset_dir = os.path.join(data_dir, dataset_name)
        actual_data_dir = prepare_dataset(dataset_info["path"], dataset_dir, dataset_info)

        for field, field_info in dataset_info["fields"].items():
            dims = field_info["dims"]
            dtype = field_info.get("dtype", "float32")
            
            for algo in algorithms:
                for eb in error_bounds:
                    data_file = os.path.join(actual_data_dir, field)
                    if not os.path.isfile(data_file):
                        print(f"Data file {data_file} does not exist. Skipping")
                        continue

                    # Call H5 test
                    cmd = [sys.executable, os.path.join(script_dir, "test_h5_filter.py"), h5_plugin_path, algo, str(eb), data_file, dtype] + [str(d) for d in dims]
                    results.append(run_test(cmd, f"Testing HDF5 {algo} {eb} on {dataset_name}/{field}"))

                    # Call SZ3 test
                    cmd = [sys.executable, os.path.join(script_dir, "test_sz3_executable.py"), sz3_executable_path, algo, str(eb), data_file, dtype] + [str(d) for d in dims]
                    results.append(run_test(cmd, f"Testing SZ3 EXE {algo} {eb} on {dataset_name}/{field}"))
        
        if os.getenv('GITHUB_ACTIONS') == 'true':
            shutil.rmtree(dataset_dir)

    # Summary
    total_tests = len(results)
    passed = sum(results)
    failed = total_tests - passed
    print(f"\nSummary: {passed}/{total_tests} tests passed, {failed} failed.")
    if failed > 0:
        print("Some tests failed. Exiting with error.")
        sys.exit(1)

    os.chdir(original_wd)

if __name__ == "__main__":
    main()
