import csv
import subprocess
import os
import sys
import json
import requests
import tarfile
import tempfile
import shutil
import time
import random

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
        return True, result.stdout
    else:
        print("FAIL")
        print("="*80)
        return False, result.stdout


BASELINE_TOLERANCE = 0.05


def baseline_key(dataset, field, algo, eb, harness):
    return f"{dataset}/{field}|{algo}|{eb}|{harness}"


def load_baseline(path):
    """Compressed sizes a previous run produced. Absent means nothing is checked."""
    try:
        with open(path) as handle:
            return {k: v for k, v in json.load(handle).items() if not k.startswith("_")}
    except FileNotFoundError:
        print(f"No compression baseline at {path}; sizes will be reported, not checked")
        return {}


def check_against_baseline(baseline, row):
    """Compression is deterministic, so a size that moved means the output moved.

    Returns True when the case is within tolerance or has no baseline to compare against.
    """
    key = baseline_key(row["dataset"], row["field"], row["algo"], row["error_bound"], row["harness"])
    expected = baseline.get(key)
    if expected is None or not row["bytes"]:
        if expected is None:
            print(f"BASELINE new {key} = {row['bytes']}")
        return True
    actual = int(row["bytes"])
    drift = actual / expected - 1
    if abs(drift) <= BASELINE_TOLERANCE:
        print(f"BASELINE ok {key}: {actual} vs {expected} ({drift:+.2%})")
        return True
    print(f"BASELINE FAIL {key}: {actual} vs {expected} ({drift:+.2%}), "
          f"outside +/-{BASELINE_TOLERANCE:.0%}")
    return False


def collect_metrics(dataset, field, algo, eb, harness, output):
    """Every METRICS line a test printed. The HDF5 test prints one per chunk mode."""
    row = {"dataset": dataset, "field": field, "algo": algo, "error_bound": eb, "harness": harness,
           "bytes": "", "ratio": "", "compress_s": "", "decompress_s": "", "max_error": ""}
    rows = []
    for line in (output or "").splitlines():
        if not line.startswith("METRICS "):
            continue
        found = dict(row)
        for pair in line[len("METRICS "):].split():
            key, _, value = pair.partition("=")
            if key in found:
                found[key] = value
            elif key == "chunk":
                found["harness"] = f"{harness}:{value}"
        rows.append(found)
    return rows


def write_metrics(rows, path):
    """A ratio and a timing per case, so a regression in either is visible between runs."""
    if not rows:
        return
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows of compression metrics to {path}")


def prepare_dataset(path, dataset_dir, dataset_info=None):
    """
    Prepares the dataset: "mdtraj:<name>" fetches and converts a published MD trajectory,
    an http path is downloaded and extracted, a local path is copied.
    Returns the actual directory that directly containing the files.
    """
    if path.startswith('mdtraj:'):
        # A published trajectory, which fetch_md_trajectory.py downloads and converts to the raw
        # arrays every other dataset here already is. Only the fields this run asks for are
        # written, because a whole trajectory in every layout does not fit a runner's disk.
        os.makedirs(dataset_dir, exist_ok=True)
        wanted = sorted({f.rsplit('.', 1)[0] for f in (dataset_info or {}).get("fields", {})})
        script = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fetch_md_trajectory.py')
        cmd = [sys.executable, script, path.split(':', 1)[1], dataset_dir]
        if wanted:
            cmd += ['--fields', ','.join(wanted)]
        if (dataset_info or {}).get("max_frames"):
            cmd += ['--max-frames', str(dataset_info["max_frames"])]
        subprocess.run(cmd, check=True)
        return dataset_dir

    if path.startswith('http'):
        # Download and extract
        if not os.path.exists(dataset_dir):
            os.makedirs(dataset_dir)

        tar_filename = os.path.join(dataset_dir, os.path.basename(path))

        if not os.path.exists(tar_filename):
            print(f"Downloading {path} to {tar_filename}")
            for attempt in range(5):
                try:
                    with requests.get(path, stream=True) as r:
                        r.raise_for_status()
                        with open(tar_filename, 'wb') as f:
                            for chunk in r.iter_content(chunk_size=8192):
                                f.write(chunk)
                    break  # Success
                except requests.exceptions.ConnectionError as e:
                    if attempt < 4:
                        print(f"Download failed (attempt {attempt+1}/5), retrying in 5 seconds...")
                        time.sleep(random.uniform(30, 60))
                    else:
                        raise
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
        # "name" runs the whole dataset, "name:a.f32,b.f32" runs those fields of it, so a dataset
        # too slow to be one job can be spread over several.
        selected_dataset, _, selected_fields = sys.argv[2].partition(":")
        if selected_dataset not in datasets:
            print(f"Dataset {selected_dataset} not found in {datasets_json}")
            sys.exit(1)
        dataset_info = datasets[selected_dataset]
        if selected_fields:
            wanted = [f for f in selected_fields.split(",") if f]
            missing = [f for f in wanted if f not in dataset_info["fields"]]
            if missing:
                print(f"Fields {missing} not found in dataset {selected_dataset}")
                sys.exit(1)
            dataset_info = dict(dataset_info)
            dataset_info["fields"] = {f: dataset_info["fields"][f] for f in wanted}
        datasets = {selected_dataset: dataset_info}

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
    metrics = []
    baseline = load_baseline(os.path.join(script_dir, "compression_baseline.json"))

    try:
        run_datasets(datasets, data_dir, script_dir, sz3_executable_path, h5_plugin_path,
                     algorithms, error_bounds, baseline, results, metrics)
    finally:
        write_metrics(metrics, os.path.join(original_wd, "integration_metrics.csv"))

    # Summary
    total_tests = len(results)
    passed = sum(results)
    failed = total_tests - passed
    print(f"\nSummary: {passed}/{total_tests} tests passed, {failed} failed.")
    if failed > 0:
        print("Some tests failed. Exiting with error.")
        sys.exit(1)

    os.chdir(original_wd)


def run_datasets(datasets, data_dir, script_dir, sz3_executable_path, h5_plugin_path,
                 algorithms, error_bounds, baseline, results, metrics):
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
                    passed, output = run_test(cmd, f"Testing HDF5 {algo} {eb} on {dataset_name}/{field}")
                    results.append(passed)
                    metrics.extend(collect_metrics(dataset_name, field, algo, eb, "hdf5", output))

                    # Call SZ3 test
                    cmd = [sys.executable, os.path.join(script_dir, "test_sz3_executable.py"), sz3_executable_path, algo, str(eb), data_file, dtype] + [str(d) for d in dims]
                    passed, output = run_test(cmd, f"Testing SZ3 EXE {algo} {eb} on {dataset_name}/{field}")
                    rows = collect_metrics(dataset_name, field, algo, eb, "cli", output)
                    metrics.extend(rows)
                    results.append(passed and all(check_against_baseline(baseline, r) for r in rows))
        
        if os.getenv('GITHUB_ACTIONS') == 'true':
            shutil.rmtree(dataset_dir)


if __name__ == "__main__":
    main()
